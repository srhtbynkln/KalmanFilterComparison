function kf_out = KF(ds)
    % Veri uzunluğu ve ön bellek tahsisi
    N = length(ds.time);
    kf_out.n = zeros(N, 1);
    kf_out.e = zeros(N, 1);
    kf_out.vn = zeros(N, 1);
    kf_out.ve = zeros(N, 1);
    kf_out.lat = zeros(N, 1);
    kf_out.lon = zeros(N, 1);
    kf_out.time = ds.time;
    
    % Referans koordinatlar (Metreye çevirmek için)
    ref_lat = ds.gps.lat(1);
    ref_lon = ds.gps.lon(1);
    
    % Durum Vektörü: x = [Kuzey (m); Doğu (m); Hız_Kuzey (m/s); Hız_Doğu (m/s)]
    x = [0; 0; 0; 0]; 
    
    % Yönelim (Yaw) açısı lineer filtrede durum vektöründe olamaz, dışarıda tutuyoruz
    psi = 0;
    if isfield(ds.gps, 'heading') && ~isnan(ds.gps.heading(1))
        psi = ds.gps.heading(1) * pi / 180;
    end
    
    dt = ds.dt;
    
    % --- SABİT LİNEER MATRİSLER ---
    % Durum Geçiş Matrisi (A)
    A = [1, 0, dt,  0; 
         0, 1,  0, dt; 
         0, 0,  1,  0; 
         0, 0,  0,  1];
     
    % Kontrol Giriş Matrisi (B)
    B = [0.5*dt^2,        0; 
                0, 0.5*dt^2; 
               dt,        0; 
                0,       dt];
            
    % Ölçüm Matrisi (H) - Sadece Kuzey ve Doğu konumunu ölçüyoruz
    H = [1, 0, 0, 0; 
         0, 1, 0, 0];
         
    % Kovaryans ve Gürültü Matrisleri
    P = diag([10, 10, 1, 1]); 
    %imu datasheetine göre Q matrisi
    acc_noise_std = 0.0006; 
    % Q = diag([0.01, 0.01, acc_noise_std^2, acc_noise_std^2]);
    %senin Q matrisin daha doğru sonuç veriyor ama nedenini anlamaıdm
    Q = diag([0.1, 0.1, 0.5, 0.5]); 
    %RTK'li olduğundan gnss verisini daha çok takip etmesini istedim r
    %matrisini küçülttüm
    R = diag([1, 1 ]); % Varsayılan 1m hata (1^2)
    
    gps_idx = 1;
    num_gps = length(ds.gps.time);
    
    % GPS verilerini baştan metreye çevir
    gps_n = zeros(num_gps, 1);
    gps_e = zeros(num_gps, 1);
    for i = 1:num_gps
        [gps_n(i), gps_e(i)] = lla2ned_local(ds.gps.lat(i), ds.gps.lon(i), ref_lat, ref_lon);
    end
    
    % Birim dönüşüm katsayıları
    % NOT: ivme zaten m/s^2 (telefon=lineer, modena=ham). g->m/s^2 olceklemesi YAPILMAZ.
    scale_acc = 1.0;
    scale_gyro = 1.0;
    if mean(abs(ds.gz)) > 10, scale_gyro = pi/180; end

    % --- LİNEER KF DÖNGÜSÜ ---
    for k = 1:N
        % 1. Yönelim Hesabı (Filtre dışında basit entegrasyon)
        gz = ds.gz(k) * scale_gyro;
        psi = psi + gz * dt;
        
        % 2. İvmeleri Kuzey-Doğu Eksenine Çevirme
        ax = ds.ax(k) * scale_acc; 
        ay = ds.ay(k) * scale_acc;
        
        aN = ax * cos(psi) - ay * sin(psi);
        aE = ax * sin(psi) + ay * cos(psi);
        
        u = [aN; aE]; % Kontrol Vektörü
        
        % 3. TAHMİN (PREDICTION) - Lineer Eşitlikler
        x = A * x + B * u;
        P = A * P * A' + Q;
        
        % 4. GÜNCELLEME (UPDATE)
        if gps_idx <= num_gps && ds.time(k) >= ds.gps.time(gps_idx)
            z = [gps_n(gps_idx); gps_e(gps_idx)];
            
            % Dinamik GNSS hatası (varsa)
            if isfield(ds.gps, 'hacc') && ds.gps.hacc(gps_idx) > 0
                pos_acc = ds.gps.hacc(gps_idx);
                R = diag([pos_acc^2, pos_acc^2]);
            end
            
            % Standart Kalman Kazancı ve Güncellemesi
            S = H * P * H' + R;
            K = P * H' / S;
            
            y = z - H * x; % İnovasyon
            x = x + K * y;
            P = (eye(4) - K * H) * P;
            
            gps_idx = gps_idx + 1;
        end
        
        % Sonuçları Kaydet
        kf_out.n(k) = x(1);
        kf_out.e(k) = x(2);
        kf_out.vn(k) = x(3);
        kf_out.ve(k) = x(4);
        
        [lat_k, lon_k] = ned2lla_local_func(x(1), x(2), ref_lat, ref_lon);
        kf_out.lat(k) = lat_k;
        kf_out.lon(k) = lon_k;
    end
end

% Yardımcı Fonksiyonlar
function [n, e] = lla2ned_local(lat, lon, lat0, lon0)
    Re = 6378137; deg2rad = pi/180;
    n = (lat - lat0) * deg2rad * Re;
    e = (lon - lon0) * deg2rad * Re * cos(lat0 * deg2rad);
end

function [lat, lon] = ned2lla_local_func(n, e, lat0, lon0)
    Re = 6378137; rad2deg = 180/pi;
    lat = lat0 + (n / Re) * rad2deg;
    lon = lon0 + (e / (Re * cos(lat0 * pi/180))) * rad2deg;
end
