function plotResultsAll(ds, kf, sh, stf, ekf, ekfSage, outputDir)
    % Create output directory if it doesn't exist
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    %% Figure 1: Geospatial Map Comparison
    f1 = figure('Color', 'w', 'Position', [100 100 1200 800]);
    
    % Plot Raw GPS
    geoplot(ds.gps.lat, ds.gps.lon, 'ko', 'MarkerSize', 3, 'DisplayName', 'GPS Raw'); 
    hold on;
    
    % Her filtre için farklı renk ve çizgi stilleri
    geoplot(ekf.lat, ekf.lon, 'g-', 'LineWidth', 2, 'DisplayName', 'Classic EKF');             % Yeşil
    geoplot(ekfSage.lat, ekfSage.lon, 'm--', 'LineWidth', 2, 'DisplayName', 'EKF Sage-Husa');   % Macenta (Eflatun)
    geoplot(kf.lat, kf.lon, 'c.-', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'KF');    % Siyan (Açık Mavi)
    geoplot(sh.lat, sh.lon, 'b-', 'LineWidth', 2, 'DisplayName', 'Sage-Husa AKF');           % Mavi
    geoplot(stf.lat, stf.lon, 'r-', 'LineWidth', 2, 'DisplayName', 'Strong Tracking Filter'); % Kırmızı
    
    % Görünümü netleştirmek için
    legend('show', 'Location', 'best');

    geobasemap('streets'); 
    title('Trajectory Comparison: Kalman Filter Variants');
    
    saveas(f1, fullfile(outputDir, 'Map_comparison.png'));
    
    %% Figure 2: North Position Time Series
    f2 = figure('Visible', 'on', 'Color', 'w'); % Set to 'on' to see it, or 'off' for headless
    
    % Plot all variants for a true comparison
    plot(ds.gps.time, ds.gps.lat, 'k:', 'LineWidth', 1.5, 'DisplayName', 'GPS Raw'); hold on;
    plot(ekf.time, ekf.n, 'g', 'LineWidth', 1.5, 'DisplayName', 'Classic EKF');
    plot(ekfSage.time, ekfSage.n, 'r', 'LineWidth', 1.5, 'DisplayName', 'EKF Sage-Husa');
    plot(EKFst.time, EKFst.n, 'b', 'LineWidth', 1.5, 'DisplayName', 'EKF Strong Tracking');
    
    title('North Position (m) Comparison');
    xlabel('Time (s)'); 
    ylabel('Distance (m)');
    legend('Location', 'northeast');
    grid on;
    
    saveas(f2, fullfile(outputDir, 'north_pos_comparison.png'));
    
    fprintf('Plots successfully saved to: %s\n', outputDir);
end