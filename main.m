%% ===================================================================
%  KalmanFilterComparison - ANA SCRIPT
%  GNSS/IMU fuzyonu ile konum & hiz kestirimi (tren + telefon veri setleri)
% ===================================================================
clear; clc; close all;

try
    scriptPath = fileparts(mfilename('fullpath'));
    if ~isempty(scriptPath), cd(scriptPath); end
catch
    disp('Uyari: dizin ayarlanamadi');
end
addpath(genpath(pwd));
if ~exist('Plot', 'dir'), mkdir('Plot'); end

%% --- VERI SETI SECIMI ---
% Telefon (Marmaray):  'data_marmaray/'
% Tren  (Modena RTK):  'data/'
dataDir = 'data_marmaray/';
fprintf('Veri dizini: %s\n', dataDir);

[dataStruct, T0] = loadAndSyncSensors(dataDir);
fprintf('IMU ornek: %d (%.0f Hz),  GPS ornek: %d (%.2f Hz),  Sure: %.0f s\n', ...
    numel(dataStruct.time), 1/dataStruct.dt, numel(dataStruct.gps.time), ...
    1/mean(diff(dataStruct.gps.time)), dataStruct.time(end)-dataStruct.time(1));

%% --- FILTRELERI CALISTIR ---
KF_out        = KF(dataStruct);
SH_out        = SageHusaKF(dataStruct);
STF_out       = STF(dataStruct);
EKF_out       = EKF(dataStruct);
EKF_sage_out  = EKF_sagehusa(dataStruct);
EKF_stf_out   = EKF_strong_tracking(dataStruct);   % outlier-gate'li (artik dahil)
imu_out       = imuLoc(dataStruct);                % saf DR (drift gosterir)
fusion_out    = fusionEKF(dataStruct);             % ANA COZUM: GPS pos+HIZ fuzyonu

%% --- DOGRULUK METRIKLERI ---
results = {KF_out, SH_out, STF_out, EKF_out, EKF_sage_out, EKF_stf_out, imu_out, fusion_out};
names   = {'KF','SageHusaKF','STF','EKF','EKF-SageHusa','EKF-StrongTrack','imuLoc(DR)','fusionEKF'};
metrics = computeMetrics(dataStruct, results, names);
writetable(metrics, 'metrics.csv');

%% --- TAHMIN CSV'leri ---
writetable(table(fusion_out.lat, fusion_out.lon, fusion_out.v, ...
    'VariableNames', {'Lat','Lon','Speed_ms'}), 'fusion_predictions.csv');
writetable(table(imu_out.lat, imu_out.lon, imu_out.v, ...
    'VariableNames', {'Lat','Lon','Speed_ms'}), 'imu_predictions.csv');

%% --- ODAK GRAFIK: fusion vs DR vs GPS (konum + hiz) ---
plotFusionVsDR(dataStruct, fusion_out, imu_out, 'Plot/');

%% --- (opsiyonel) mevcut tekil grafikler ---
try
    plotResultsKF(dataStruct, KF_out, 'Plot/');
    plotResultsSTF(dataStruct, STF_out, 'Plot/');
    plotResultsEKF(dataStruct, EKF_out, 'Plot/');
    plotResultsimu(dataStruct, imu_out, 'Plot/');
catch ME
    fprintf('Grafik uyarisi (Mapping Toolbox gerekebilir): %s\n', ME.message);
end

fprintf('\nBitti. Sonuc: metrics.csv ve Plot/ klasoru.\n');
fprintf('Ana cikti: fusionEKF -> konum + HIZ (bkz. REHBER.md)\n');

%% ===================================================================
function plotFusionVsDR(ds, fus, imu, outDir)
    if ~exist(outDir,'dir'), mkdir(outDir); end
    Re = 6378137; d2r = pi/180;
    lat0 = ds.gps.lat(1); lon0 = ds.gps.lon(1);
    gN = (ds.gps.lat-lat0)*d2r*Re;
    gE = (ds.gps.lon-lon0)*d2r*Re*cos(lat0*d2r);

    f = figure('Visible','on','Position',[80 80 1300 850]);

    % 1) Yorunge (NED)
    subplot(2,2,1); hold on; grid on; axis equal;
    plot(gE, gN, 'k.', 'DisplayName','GPS');
    plot(fus.e, fus.n, 'b-', 'LineWidth',1.3, 'DisplayName','fusionEKF');
    plot(imu.e, imu.n, 'r-', 'LineWidth',1.0, 'DisplayName','imuLoc (saf DR)');
    xlabel('Dogu (m)'); ylabel('Kuzey (m)'); title('Yorunge'); legend('Location','best');

    % 2) HIZ karsilastirmasi
    subplot(2,2,2); hold on; grid on;
    if isfield(ds.gps,'speed')
        plot(ds.gps.time, ds.gps.speed, 'k.', 'DisplayName','GPS speed');
    end
    plot(fus.time, fus.v, 'b-', 'LineWidth',1.2, 'DisplayName','fusionEKF');
    plot(imu.time, imu.v, 'r-', 'LineWidth',0.8, 'DisplayName','imuLoc (DR)');
    xlabel('Zaman (s)'); ylabel('Hiz (m/s)'); title('HIZ'); legend('Location','best');

    % 3) Konum hatasi (GPS'e gore)
    subplot(2,2,3); hold on; grid on;
    fn = interp1(fus.time, fus.n, ds.gps.time); fe = interp1(fus.time, fus.e, ds.gps.time);
    inx= interp1(imu.time, imu.n, ds.gps.time); iex= interp1(imu.time, imu.e, ds.gps.time);
    plot(ds.gps.time, sqrt((fn-gN).^2+(fe-gE).^2), 'b-', 'DisplayName','fusionEKF');
    plot(ds.gps.time, sqrt((inx-gN).^2+(iex-gE).^2), 'r-', 'DisplayName','imuLoc (DR)');
    xlabel('Zaman (s)'); ylabel('Konum hatasi (m)'); title('GPS''e gore konum hatasi');
    legend('Location','best');

    % 4) Kestirilen ivme bias'i (fusion)
    subplot(2,2,4); hold on; grid on;
    plot(fus.time, fus.bN, 'DisplayName','b_{N}');
    plot(fus.time, fus.bE, 'DisplayName','b_{E}');
    xlabel('Zaman (s)'); ylabel('Bias (m/s^2)'); title('fusionEKF - kestirilen ivme bias''i');
    legend('Location','best');

    saveas(f, fullfile(outDir,'fusion_vs_dr.png'));
end
