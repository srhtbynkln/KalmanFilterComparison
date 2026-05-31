function plotResultsimu(ds, imu, outputDir)
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    imu.time = ds.time;
    % 1. Harita Çizimi (Enlem - Boylam)
    f1 = figure('Visible','on', 'Position', [100 100 1200 800]);
    
    
    geoplot(ds.gps.lat, ds.gps.lon, 'k.', 'DisplayName', 'GPS Raw', 'MarkerSize', 8); hold on;
    geoplot(imu.lat, imu.lon, 'g.-', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'imu');
    geobasemap('streets'); 
    title('Kalman Filter - Harita Üzerinde Konum');
    legend('Location', 'best');
    saveas(f1, fullfile(outputDir, 'imu_map.png'));
    
    % 2. X ve Y Eksenlerinde Metre Cinsinden Karşılaştırma Grafiği
    f2 = figure('Visible','on', 'Position', [150 150 1000 600]);
    
    % GPS verilerini karşılaştırma için metreye çevir
    Re = 6378137;
    deg2rad = pi/180;
    lat0 = ds.gps.lat(1);
    lon0 = ds.gps.lon(1);
    gps_n = (ds.gps.lat - lat0) * deg2rad * Re;
    gps_e = (ds.gps.lon - lon0) * deg2rad * Re .* cos(lat0 * deg2rad);
    
    % Kuzey (North)
    subplot(2,1,1);
    plot(imu.time, imu.n, 'g', 'LineWidth', 1.5); hold on;
    plot(ds.gps.time, gps_n, 'k:', 'LineWidth', 1.5);
    title('Kuzey (North) Ekseni Pozisyonu'); 
    xlabel('Zaman (s)'); ylabel('Mesafe (m)');
    legend('imu', 'GPS Raw', 'Location', 'best');
    grid on;
    
    % Doğu (East)
    subplot(2,1,2);
    plot(imu.time, imu.e, 'g', 'LineWidth', 1.5); hold on;
    plot(ds.gps.time, gps_e, 'k:', 'LineWidth', 1.5);
    title('Doğu (East) Ekseni Pozisyonu'); 
    xlabel('Zaman (s)'); ylabel('Mesafe (m)');
    legend('imu', 'GPS Raw', 'Location', 'best');
    grid on;
    
    saveas(f2, fullfile(outputDir, 'pos_comparison.png'));
end