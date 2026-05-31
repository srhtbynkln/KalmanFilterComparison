function out = EKF_sagehusa(ds)
    N = length(ds.time);
    dt = mean(diff(ds.time));
    
    lat0 = ds.gps.lat(1); lon0 = ds.gps.lon(1); h0 = ds.gps.h(1);
    Re = 6378137; deg2rad_const = pi/180;
    
    gps_n = (ds.gps.lat - lat0) * (Re * deg2rad_const);
    gps_e = (ds.gps.lon - lon0) * (Re * cos(lat0 * deg2rad_const) * deg2rad_const);
    
    x = [0; 0; 0; ds.gps.heading(1)]; 
    P = eye(4) * 5;
    Q = diag([0.05, 0.05, 0.1, 0.05].^2); 
    
    est_x = zeros(4, N);
    gps_idx = 1;
    prev_gps = [];
    alpha = 0.995; % Unutma faktörü
    R_adaptive = diag([1, 1, deg2rad(5)^2]); % Başlangıç 3x3
    k_sage = 1;

    for k = 1:N
        a_res = ds.ax(k);
        dot_p = ds.ax(k)*cos(x(4)) + ds.ay(k)*sin(x(4));
        if dot_p < 0, a_res = -a_res; end
        
        [x, P] = predict_step(x, P, a_res, ds.gyr_z(k), Q, dt);
        
        if gps_idx <= length(ds.gps.time)
            time_diff = ds.time(k) - ds.gps.time(gps_idx);
            
            if time_diff >= 0 && time_diff < dt
                z_current = [gps_n(gps_idx); gps_e(gps_idx)];
                
                if ~isempty(prev_gps) && norm(z_current - prev_gps) > 0.1
                    % 3 Ölçümlü Durum
                    yaw_gps = atan2(z_current(2)-prev_gps(2), z_current(1)-prev_gps(1));
                    z_meas = [z_current; yaw_gps];
                    H_meas = [1 0 0 0; 0 1 0 0; 0 0 0 1];
                else
                    % 2 Ölçümlü Durum
                    z_meas = z_current;
                    H_meas = [1 0 0 0; 0 1 0 0];
                end
                
                % Sage-Husa Güncellemesi
                [x, P, R_adaptive] = update_step(x, P, z_meas, R_adaptive, H_meas, alpha, k_sage);
                k_sage = k_sage + 1;
                
                prev_gps = z_current;
                gps_idx = gps_idx + 1;
            end
        end
        est_x(:, k) = x;
    end
    
    out = ds;
    out.n = est_x(1,:)'; out.e = est_x(2,:)';
    out.v = est_x(3,:)'; out.yaw_est = rad2deg(est_x(4,:)');
    [out.lat, out.lon] = ned2lla_local(out.n, out.e, lat0, lon0);
end

function [x, P] = predict_step(x, P, a_res, gyr_z, Q, dt)
    F = [ 1   0   cos(x(4))*dt   -x(3)*sin(x(4))*dt ;
          0   1   sin(x(4))*dt    x(3)*cos(x(4))*dt ;
          0   0   1               0 ;
          0   0   0               1 ];
    x(1) = x(1) + x(3) * cos(x(4)) * dt;
    x(2) = x(2) + x(3) * sin(x(4)) * dt;
    x(3) = x(3) + a_res * dt;
    x(4) = x(4) + gyr_z * dt;
    P = F * P * F' + Q;
end

function [x, P, R_adaptive_new] = update_step(x, P, z, R_prev, H, alpha, k)
    % 1. İnovasyon
    y = z - H*x;
    if size(z,1) == 3 
        y(3) = atan2(sin(y(3)), cos(y(3)));
    end
    
    % 2. d katsayısı
    d = (1 - alpha) / (1 - alpha^k);
    
    % 3. Mevcut ölçüm boyutuna göre R güncelleme
    m = size(z, 1); % Ölçüm sayısı (2 veya 3)
    H_sub = H(1:m, :);
    R_sub_prev = R_prev(1:m, 1:m);
    innovation_cov = y * y';
    predicted_cov = H_sub * P * H_sub';
    
    for i = 1:m
        % Sadece diyagonal elemanları adapte et
        R_sub_new(i,i) = (1-d)*R_sub_prev(i,i) + d*(innovation_cov(i,i) - predicted_cov(i,i));
        R_sub_new(i,i) = max(R_sub_new(i,i), 1e-3); % Alt sınır koruması
    end
    % Diyagonal dışı elemanları temizle (opsiyonel ama stabiliteyi artırır)
    R_sub_new = diag(diag(R_sub_new));

    min_R = 1e-4; % 0.01 metrelik minimum standart sapma (varyans için 1e-4)
    for i = 1:size(R_sub_new, 1)
        if R_sub_new(i,i) < min_R
            R_sub_new(i,i) = min_R; 
        end
    end
    
    
    
    % R_adaptive matrisini güncelle (boyuta göre)
    R_adaptive_new = R_prev;
    R_adaptive_new(1:m, 1:m) = R_sub_new;
    
    % 4. Standart Güncelleme
    S = H_sub * P * H_sub' + R_sub_new;
    K = (P * H_sub') / S;
    x = x + K * y;
    P = (eye(size(P,1)) - K * H_sub) * P;
end

function [lat, lon] = ned2lla_local(n, e, lat0, lon0)
    Re = 6378137; deg2rad = pi/180;
    lat = lat0 + (n / Re) / deg2rad;
    lon = lon0 + (e / (Re * cos(lat0 * deg2rad))) / deg2rad;
end