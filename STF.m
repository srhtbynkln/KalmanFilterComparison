function stf_out = STF(ds)
    % Veri uzunluğu ve ön bellek tahsisi
    N = length(ds.time);
    stf_out.n = zeros(N, 1);
    stf_out.e = zeros(N, 1);
    stf_out.vn = zeros(N, 1);
    stf_out.ve = zeros(N, 1);
    stf_out.lat = zeros(N, 1);
    stf_out.lon = zeros(N, 1);
    stf_out.time = ds.time;
    
    ref_lat = ds.gps.lat(1);
    ref_lon = ds.gps.lon(1);
    
    % Durum Vektörü: x = [Kuzey; Doğu; Hız_N; Hız_E]
    x = [0; 0; 0; 0]; 
    psi = 0;
    if isfield(ds.gps, 'heading') && ~isnan(ds.gps.heading(1))
        psi = ds.gps.heading(1) * pi / 180;
    end
    
    dt = ds.dt;
    
    % --- MATRİSLER ---
    A = [1, 0, dt,  0; 
         0, 1,  0, dt; 
         0, 0,  1,  0; 
         0, 0,  0,  1];
     
    B = [0.5*dt^2, 0; 0, 0.5*dt^2; dt, 0; 0, dt];
    H = [1, 0, 0, 0; 0, 1, 0, 0];
         
    P = diag([10, 10, 1, 1]); 
    Q = diag([0.1, 0.1, 0.5, 0.5]); 
    R = diag([1, 1]); 
    
    % --- STF PARAMETRELERİ [cite: 279, 281] ---
    rho = 0.95; % Unutma faktörü (p)
    V_o = zeros(2, 2); 
    
    gps_idx = 1;
    num_gps = length(ds.gps.time);
    
    % Birim dönüşümleri  (NOT: ivme zaten m/s^2; g->m/s^2 olceklemesi YAPILMAZ)
    scale_acc = 1.0; scale_gyro = 1.0;
    if mean(abs(ds.gz)) > 10, scale_gyro = pi/180; end

    % --- DÖNGÜ ---
    for k = 1:N
        % 1. Tahmin (Prediction)
        gz = ds.gz(k) * scale_gyro;
        psi = psi + gz * dt;
        ax = ds.ax(k) * scale_acc; ay = ds.ay(k) * scale_acc;
        aN = ax * cos(psi) - ay * sin(psi);
        aE = ax * sin(psi) + ay * cos(psi);
        u = [aN; aE];
        
        x_pred = A * x + B * u;
        
        % 2. Güncelleme ve STF Mantığı
        if gps_idx <= num_gps && ds.time(k) >= ds.gps.time(gps_idx)
            % GPS LLA -> NED dönüşümü (Hata burada oluyordu, fonksiyon alta eklendi)
            [z_n, z_e] = lla2ned_local_internal(ds.gps.lat(gps_idx), ds.gps.lon(gps_idx), ref_lat, ref_lon);
            z = [z_n; z_e];
            
            % İnovasyon (Residual)
            v = z - H * x_pred; 
            
            % STF Fading Factor Hesabı 
            if k == 1
                V_o = v * v';
            else
                V_o = (rho * V_o + v * v') / (1 + rho); 
            end
            
            N_mat = V_o - R - H * Q * H'; 
            M_mat = H * (A * P * A') * H'; 
            
            lambda = trace(N_mat) / trace(M_mat);
            if lambda < 1, lambda = 1; end % Abrupt change yoksa lambda=1 [cite: 272]
            
            % Tahmin Kovaryansını Lambda ile genişlet [cite: 263]
            P_pred = lambda * (A * P * A') + Q;
            
            % Standart Update adımları
            S = H * P_pred * H' + R;
            K = P_pred * H' / S;
            x = x_pred + K * v;
            P = (eye(4) - K * H) * P_pred;
            
            gps_idx = gps_idx + 1;
        else
            x = x_pred;
            P = A * P * A' + Q;
        end
        
        % Çıktıları kaydet
        stf_out.n(k) = x(1); stf_out.e(k) = x(2);
        stf_out.vn(k) = x(3); stf_out.ve(k) = x(4);
        [stf_out.lat(k), stf_out.lon(k)] = ned2lla_local_internal(x(1), x(2), ref_lat, ref_lon);
    end
end

% --- YARDIMCI FONKSİYONLAR (Dosya içine gömülü) ---
function [n, e] = lla2ned_local_internal(lat, lon, lat0, lon0)
    Re = 6378137; deg2rad = pi/180;
    n = (lat - lat0) * deg2rad * Re;
    e = (lon - lon0) * deg2rad * Re * cos(lat0 * deg2rad);
end

function [lat, lon] = ned2lla_local_internal(n, e, lat0, lon0)
    Re = 6378137; rad2deg = 180/pi;
    lat = lat0 + (n / Re) * rad2deg;
    lon = lon0 + (e / (Re * cos(lat0 * pi/180))) * rad2deg;
end
