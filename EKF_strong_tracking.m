function out = EKF_strong_tracking(ds)
    N = length(ds.time);
    dt = mean(diff(ds.time));
    
    % --- İlk Değerler ---
    lat0 = ds.gps.lat(1); lon0 = ds.gps.lon(1);
    Re = 6378137; deg2rad_const = pi/180;
    
    gps_n = (ds.gps.lat - lat0) * (Re * deg2rad_const);
    gps_e = (ds.gps.lon - lon0) * (Re * cos(lat0 * deg2rad_const) * deg2rad_const);
    
    x = [0; 0; 0; deg2rad(ds.gps.heading(1))]; % Heading radyana çevrilmeli
    P = eye(4) * 1.0; % Başlangıç belirsizliği çok uçmasın
    Q = diag([0.05, 0.05, 0.1, 0.02].^2); 
    
    est_x = zeros(4, N);
    gps_idx = 1;
    prev_gps = [];
    
    % Strong Tracking için gerekli inovasyon hafızası
    y_last = zeros(2,1); % Başlangıçta 2D (konum) kabul edelim
    
    for k = 1:N
        % İvme işleme
        % Yercekimi sadece ham IMU verisinde (modena) cikarilir.
        if ds.gravity_present, az_lin = ds.az(k) - 9.81; else, az_lin = ds.az(k); end
        a_res = sqrt(ds.ax(k)^2 + ds.ay(k)^2 + az_lin^2);
        dot_p = ds.ax(k)*cos(x(4)) + ds.ay(k)*sin(x(4));
        if dot_p < 0, a_res = -a_res; end
        
        % 1. Tahmin Adımı
        [x, P] = predict_step(x, P, a_res, ds.gyr_z(k), Q, dt);
        
        % 2. GPS Güncelleme Kontrolü
        if gps_idx <= length(ds.gps.time)
            time_diff = ds.time(k) - ds.gps.time(gps_idx);
            
            if time_diff >= 0 && time_diff < dt
                hacc = ds.gps.hacc(gps_idx);
                
                % Kötü GPS Reddi (Hard Gate)
                if hacc > 15 || isnan(hacc)
                    gps_idx = gps_idx + 1;
                    continue;
                end
                
                z_current = [gps_n(gps_idx); gps_e(gps_idx)];
                
                % Ölçüm Modeli Seçimi
                if ~isempty(prev_gps) && norm(z_current - prev_gps) > 0.1
                    yaw_gps = atan2(z_current(2)-prev_gps(2), z_current(1)-prev_gps(1));
                    z_meas = [z_current; yaw_gps];
                    H_meas = [1 0 0 0; 0 1 0 0; 0 0 0 1];
                    R_meas = diag([hacc^2, hacc^2, deg2rad(10)^2]);
                else
                    z_meas = z_current;
                    H_meas = [1 0 0 0; 0 1 0 0];
                    R_meas = eye(2) * hacc^2;
                end
                
                % 3. Güncelleme Adımı (Strong Tracking Mantığı İçeride)
                [x, P, y_last] = update_step_robust(x, P, z_meas, R_meas, H_meas);
                
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
    % Jacobian F
    F = [ 1   0   cos(x(4))*dt   -x(3)*sin(x(4))*dt ;
          0   1   sin(x(4))*dt    x(3)*cos(x(4))*dt ;
          0   0   1               0 ;
          0   0   0               1 ];
          
    % State Prediction
    if abs(a_res) < 0.01 && abs(gyr_z) < deg2rad(0.1)
        x(3) = 0; % Hızı sıfırla
        % P matrisindeki hız ve açı belirsizliğini küçült
        P(3,3) = 0.01; P(4,4) = 0.01; 
    end
    x(1) = x(1) + x(3) * cos(x(4)) * dt;
    x(2) = x(2) + x(3) * sin(x(4)) * dt;
    x(3) = x(3) + a_res * dt;

    x(4) = x(4) + gyr_z * dt;
    x(4) = atan2(sin(x(4)), cos(x(4)));
    
    % Covariance Prediction
    P = F * P * F' + Q;
    
    % --- Kritik: P matrisini simetrik ve pozitif tanımlı tut ---
    P = (P + P') * 0.5; 
end

function [x, P, y] = update_step_robust(x, P, z, R, H)
    y = z - H*x;
    % Açı farkını normalize et
    if size(z,1) == 3, y(3) = atan2(sin(y(3)), cos(y(3))); end
    
    S = H * P * H' + R;
    
    if rcond(S) < 1e-15
        return; 
    end
    
    d2 = y' / S * y;
    gate_limit = 11.3; % 3 DOF için %99 güven aralığı
    
    if d2 > gate_limit
        % ÖNEMLİ: Sadece R'yi şişirmek yerine P'yi de biraz serbest bırak ki 
        % filtre GPS'e geri dönebilsin.
        P = P * 1.1; 
        R = R * (d2 / gate_limit);
        S = H * P * H' + R;
    end
    
    K = (P * H') / S;
    
    % State güncelleme ve hemen ardından açı normalizasyonu
    x = x + K * y;
    x(4) = atan2(sin(x(4)), cos(x(4))); 
    
    % Joseph Formu (Burası iyi, koru)
    I = eye(size(P,1));
    P = (I - K*H) * P * (I - K*H)' + K * R * K';
    
    % Simetri ve stabilite
    P = (P + P') * 0.5 + eye(size(P,1)) * 1e-9;
end
function [lat, lon] = ned2lla_local(n, e, lat0, lon0)
    Re = 6378137; deg2rad = pi/180;
    lat = lat0 + (n / Re) / deg2rad;
    lon = lon0 + (e / (Re * cos(lat0 * deg2rad))) / deg2rad;
end