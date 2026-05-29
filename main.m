%%
clear; clc; close all;
%%
try
    scriptPath = fileparts(mfilename('fullpath'));
    if ~isempty(scriptPath)
        cd(scriptPath);
    end
catch
    disp('Warning');
end

addpath(genpath(pwd));
if ~exist('outputs', 'dir'), mkdir('outputs'); end
if ~exist('Plot', 'dir'), mkdir('Plot'); end

[dataStruct, T0] = loadAndSyncSensors('data/');

ref_lat = dataStruct.gps.lat(1);
ref_lon = dataStruct.gps.lon(1);
ref_alt = dataStruct.gps.h(1);
ref_yaw = dataStruct.gps.heading(1);
%% klasik kalman filtreleri
KF_out = KF(dataStruct);

SH_out = SageHusaKF(dataStruct);

STF_out = STF(dataStruct);

%% extended NonLinear kalman filtresi
EKF_out = EKF(dataStruct);

EKF_sage_out = EKF_sagehusa(dataStruct);

%%
plotResultsKF(dataStruct, KF_out, 'Plot/');

plotResultsSageHusa(dataStruct, SH_out, 'Plot/');

plotResultsSTF(dataStruct, STF_out, 'Plot/');


%%
kf_table = table(KF_out.lat, KF_out.lon, 'VariableNames', {'Lat', 'Lon'});
writetable(kf_table, 'kf_predictions.csv');

sh_table = table(SH_out.lat, SH_out.lon, 'VariableNames', {'Lat', 'Lon'});
writetable(sh_table, 'sh_predictions.csv');

stf_table = table(STF_out.lat, STF_out.lon, 'VariableNames', {'Lat', 'Lon'});
writetable(stf_table, 'stf_predictions.csv');

ekf_table = table(EKF_out.lat, EKF_out.lon, 'VariableNames', {'Lat', 'Lon'});
writetable(ekf_table, 'ekf_predictions.csv');

ekfsh_table = table(EKF_sage_out.lat, EKF_sage_out.lon, 'VariableNames', {'Lat', 'Lon'});
writetable(ekfsh_table, 'ekfsh_predictions.csv');

%%
plotResultsEKF(dataStruct, EKF_out, 'Plot/');

plotResultsEKF_sagehusa(dataStruct, EKF_sage_out, 'Plot/');

%%
plotResultsAll(dataStruct,KF_out,SH_out,STF_out ,EKF_out,EKF_sage_out, 'Plot/');
%%
functions = {@KF, @SageHusaKF, @STF, @EKF, @EKF_sagehusa};
names = {'Standard KF', 'Sage-Husa KF', 'Strong Tracking Filter', 'EKF', 'EKF-SageHusa'};
num_runs = 20; 
results = zeros(length(functions), 1);
N = length(dataStruct.time);

fprintf('=== Deney Başlatıldı ===\n');
total_tic = tic; % Toplam süre için

for i = 1:length(functions)
    fprintf('[%d/%d] %-25s çalışıyor... ', i, length(functions), names{i});
    
    iter_tic = tic;
    for r = 1:num_runs
        feval(functions{i}, dataStruct);
        disp(r);
    end
    results(i) = toc(iter_tic) / num_runs;
    
    fprintf('Bitti! (Ortalama: %.4f sn)\n', results(i));
end

total_duration = toc(total_tic);
fprintf('=== Tüm Deney Tamamlandı! Toplam Süre: %.2f sn ===\n\n', total_duration);

% Tablo oluşturma
Step_Time_ms = (results ./ N) * 1000; 
T = table(names', results, Step_Time_ms, 'VariableNames', {'Algorithm', 'Avg_Total_Time_s', 'Avg_Step_Time_ms'});
disp(T);