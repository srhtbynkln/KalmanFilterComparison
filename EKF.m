function out = EKF(ds)
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

    for k = 1:N
        % Yercekimi sadece ham IMU verisinde (modena) cikarilir; telefon verisi
        % zaten lineer (yercekimsiz) oldugundan cikarilmaz.
        if ds.gravity_present, az_lin = ds.az(k) - 9.81; else, az_lin = ds.az(k); end
        a_res = sqrt(ds.ax(k)^2 + ds.ay(k)^2 + az_lin^2);
        dot_p = ds.ax(k)*cos(x(4)) + ds.ay(k)*sin(x(4));
        if dot_p < 0, a_res = -a_res; end
        
        [x, P] = predict_step(x, P, a_res, ds.gyr_z(k), Q, dt);
        
        if gps_idx <= length(ds.gps.time)
            time_diff = ds.time(k) - ds.gps.time(gps_idx);
            
            if time_diff >= 0 && time_diff < dt
                z_current = [gps_n(gps_idx); gps_e(gps_idx)];
                
                if ~isempty(prev_gps) && norm(z_current - prev_gps) > 0.1
                    yaw_gps = atan2(z_current(2)-prev_gps(2), z_current(1)-prev_gps(1));
                    z_meas = [z_current; yaw_gps];
                    H_meas = [1 0 0 0; 0 1 0 0; 0 0 0 1];
                    R_meas = diag([ds.gps.hacc(gps_idx)^2, ds.gps.hacc(gps_idx)^2, deg2rad(5)^2]);
                else
                    z_meas = z_current;
                    H_meas = [1 0 0 0; 0 1 0 0];
                    R_meas = eye(2) * ds.gps.hacc(gps_idx)^2;
                end
                
                [x, P] = update_step(x, P, z_meas, R_meas, H_meas);
                
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

function [x, P] = update_step(x, P, z, R, H)
    y = z - H*x;
    if size(z,1) == 3 
        y(3) = atan2(sin(y(3)), cos(y(3)));
    end
    S = H * P * H' + R;
    K = (P * H') / S;
    x = x + K * y;
    P = (eye(size(P,1)) - K*H) * P;
end

function [lat, lon] = ned2lla_local(n, e, lat0, lon0)
    Re = 6378137; deg2rad = pi/180;
    lat = lat0 + (n / Re) / deg2rad;
    lon = lon0 + (e / (Re * cos(lat0 * deg2rad))) / deg2rad;
end