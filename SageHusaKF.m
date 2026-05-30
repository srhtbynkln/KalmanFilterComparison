function sh_out = SageHusaKF(ds)
    % Veri uzunluğu ve ön bellek tahsisi
    N = length(ds.time);
    sh_out.n = zeros(N, 1);
    sh_out.e = zeros(N, 1);
    sh_out.vn = zeros(N, 1);
    sh_out.ve = zeros(N, 1);
    sh_out.lat = zeros(N, 1);
    sh_out.lon = zeros(N, 1);
    sh_out.time = ds.time;
    
    ref_lat = ds.gps.lat(1);
    ref_lon = ds.gps.lon(1);
    
    % Durum Vektörü: x = [Kuzey; Doğu; Hız_N; Hız_E]
    x = [0; 0; 0; 0]; 
    psi = 0;
    if isfield(ds.gps, 'heading') && ~isnan(ds.gps.heading(1))
        psi = ds.gps.heading(1) * pi / 180;
    end
    
    dt = ds.dt;
    
    % --- BAŞLANGIÇ MATRİSLERİ ---
    A = [1, 0, dt,  0; 
         0, 1,  0, dt; 
         0, 0,  1,  0; 
         0, 0,  0,  1];
    B = [0.5*dt^2, 0; 0, 0.5*dt^2; dt, 0; 0, dt];
    H = [1, 0, 0, 0; 0, 1, 0, 0];
    
    P = diag([10, 10, 1, 1]); 
    Q = diag([0.1, 0.1, 0.5, 0.5]); % Başlangıç tahmini
    R = diag([1, 1]);              % Başlangıç tahmini
    
    % --- SAGE-HUSA PARAMETRELERİ [cite: 254, 257] ---
    b = 0.98;      % Unutma faktörü (Forgetting factor, genellikle 0.95-0.99)
    q_hat = zeros(4, 1); % Sistem gürültüsü ortalaması
    r_hat = zeros(2, 1); % Ölçüm gürültüsü ortalaması
    
    gps_idx = 1;
    num_gps = length(ds.gps.time);
    
    % NOT: ivme zaten m/s^2; g->m/s^2 olceklemesi YAPILMAZ.
    scale_acc = 1.0; scale_gyro = 1.0;
    if mean(abs(ds.gz)) > 10, scale_gyro = pi/180; end

    for k = 1:N
        % 1. TAHMİN (Prediction) [cite: 232, 233]
        gz = ds.gz(k) * scale_gyro;
        psi = psi + gz * dt;
        ax = ds.ax(k) * scale_acc; ay = ds.ay(k) * scale_acc;
        aN = ax * cos(psi) - ay * sin(psi);
        aE = ax * sin(psi) + ay * cos(psi);
        u = [aN; aE];
        
        x_pred = A * x + B * u + q_hat; % q_hat eklenerek bias düzeltilir [cite: 232]
        P_pred = A * P * A' + Q;
        
        % 2. GÜNCELLEME VE ADAPTASYON
        if gps_idx <= num_gps && ds.time(k) >= ds.gps.time(gps_idx)
            [z_n, z_e] = lla2ned_sh(ds.gps.lat(gps_idx), ds.gps.lon(gps_idx), ref_lat, ref_lon);
            z = [z_n; z_e];
            
            % İnovasyon (Residual) [cite: 245]
            v = z - H * x_pred - r_hat;
            
            % Adaptasyon katsayısı (d_k) [cite: 254]
            d_k = (1 - b) / (1 - b^k);
            
            % --- R VE Q MATRİSLERİNİN GÜNCELLENMESİ  ---
            % Ölçüm Gürültü Matrisi Tahmini (R)
            R = (1 - d_k) * R + d_k * (v * v' - H * P_pred * H');
            R = diag(max(diag(R), 1e-4)); % Numerik stabilite için yarı-pozitif tanımlı tut
            
            % Kalman Kazancı ve Durum Güncelleme [cite: 241, 243]
            S = H * P_pred * H' + R;
            K = P_pred * H' / S;
            x = x_pred + K * v;
            P = (eye(4) - K * H) * P_pred;
            
            % Sistem Gürültü Matrisi Tahmini (Q)
            Q = (1 - d_k) * Q + d_k * (K * v * v' * K' + P - A * P * A');
            Q = diag(max(diag(Q), 1e-6)); 

            % Gürültü ortalamalarının güncellenmesi 
            r_hat = (1 - d_k) * r_hat + d_k * (z - H * x_pred);
            q_hat = (1 - d_k) * q_hat + d_k * (x - A * x);
            
            gps_idx = gps_idx + 1;
        else
            x = x_pred;
            P = A * P * A' + Q;
        end
        
        % Kayıt
        sh_out.n(k) = x(1); sh_out.e(k) = x(2);
        sh_out.vn(k) = x(3); sh_out.ve(k) = x(4);
        [sh_out.lat(k), sh_out.lon(k)] = ned2lla_sh(x(1), x(2), ref_lat, ref_lon);
    end
end

% --- YARDIMCI FONKSİYONLAR ---
function [n, e] = lla2ned_sh(lat, lon, lat0, lon0)
    Re = 6378137; deg2rad = pi/180;
    n = (lat - lat0) * deg2rad * Re;
    e = (lon - lon0) * deg2rad * Re * cos(lat0 * deg2rad);
end

function [lat, lon] = ned2lla_sh(n, e, lat0, lon0)
    Re = 6378137; rad2deg = 180/pi;
    lat = lat0 + (n / Re) * rad2deg;
    lon = lon0 + (e / (Re * cos(lat0 * pi/180))) * rad2deg;
end
