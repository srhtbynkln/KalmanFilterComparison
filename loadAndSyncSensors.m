function [ds, T0] = loadAndSyncSensors(dataDir)
% LOADANDSYNCSENSORS  Sensor verilerini yukler ve ortak bir 'ds' yapisina koyar.
%
% Iki veri formati otomatik algilanir:
%   - TELEFON  (Sensor Logger):  Accelerometer.csv, Gyroscope.csv, Location.csv
%   - MODENA   (RTK tren)     :  clean_imu.csv, clean_gnss.csv
%
% Cikis (ds) alanlari:
%   .time            IMU zaman vektoru (s, 0'dan)
%   .dt              ortalama ornekleme periyodu (s)
%   .dt_k            anlik ornekleme periyodu vektoru (s)  <-- jitter icin
%   .ax .ay .az      ivme (m/s^2)
%   .gx .gy .gz .gyr_z   acisal hiz (rad/s)
%   .gravity_present (logical)  ivme yercekimi iceriyor mu (true: modena, false: telefon)
%   .gps.time .lat .lon .h .hacc .heading        GNSS
%   .gps.speed       yatay hiz buyuklugu (m/s, Doppler)  [yoksa NaN]
%   .gps.vN .gps.vE  hizin Kuzey/Dogu bilesenleri (m/s)  [yoksa NaN]
%
% NOT (kritik): Sensor Logger 'Accelerometer' = LINEER ivme (yercekimsiz) ve
% birimi zaten m/s^2'dir. Eski kodda yapilan 9.81 ile carpma HATALIYDI
% (veriyi 9.81 kat sisiriyordu) ve kaldirilmistir.

    isNewFormat = exist(fullfile(dataDir, 'Accelerometer.csv'), 'file');

    if isNewFormat
        % ================= TELEFON (Sensor Logger) =================
        disp('Veri formati: TELEFON (Sensor Logger)');
        acc = readtable(fullfile(dataDir, 'Accelerometer.csv'));
        gyr = readtable(fullfile(dataDir, 'Gyroscope.csv'));
        loc = readtable(fullfile(dataDir, 'Location.csv'));

        ds.time = acc.seconds_elapsed;

        % Ivme: lineer (yercekimsiz), m/s^2 -> oldugu gibi kullan
        ds.ax = acc.x;  ds.ay = acc.y;  ds.az = acc.z;
        ds.gravity_present = false;

        % Jiroskop (rad/s): ivme ile FARKLI uzunlukta olabilir -> ivme zamanina
        % interpole ederek hizala (aksi halde indeks tasmasi olur).
        tg_imu = gyr.seconds_elapsed;
        ds.gx = interp1(tg_imu, gyr.x, ds.time, 'linear', 'extrap');
        ds.gy = interp1(tg_imu, gyr.y, ds.time, 'linear', 'extrap');
        ds.gz = interp1(tg_imu, gyr.z, ds.time, 'linear', 'extrap');
        ds.gyr_z = ds.gz;

        % --- GNSS ---
        ds.gps.time = loc.seconds_elapsed;
        ds.gps.lat  = loc.latitude;
        ds.gps.lon  = loc.longitude;
        ds.gps.h    = loc.altitude;
        ds.gps.hacc = loc.horizontalAccuracy;

        % Heading (bearing): -1 = gecersiz -> NaN
        hd = loc.bearing;  hd(hd < 0) = NaN;
        ds.gps.heading = hd;                 % derece

        % Doppler hiz: speed (m/s) + bearing(derece) -> vN, vE
        sp = loc.speed;  sp(sp < 0) = NaN;
        ds.gps.speed = sp;
        brg = loc.bearing * pi/180;          % bearing: Kuzey'den saat yonunde
        ds.gps.vN = sp .* cos(brg);
        ds.gps.vE = sp .* sin(brg);

        % Hiz dogrulugu (varsa)
        if ismember('speedAccuracy', loc.Properties.VariableNames)
            sa = loc.speedAccuracy;  sa(sa < 0) = NaN;
            ds.gps.speedAcc = sa;
        else
            ds.gps.speedAcc = nan(height(loc), 1);
        end

        T0 = 0;

    else
        % ================= MODENA (RTK tren) =================
        disp('Veri formati: MODENA (RTK tren)');
        imu_file  = fullfile(dataDir, 'clean_imu.csv');
        gnss_file = fullfile(dataDir, 'clean_gnss.csv');
        if ~exist(imu_file, 'file'),  error('Bulunamadi: %s', imu_file);  end
        if ~exist(gnss_file, 'file'), error('Bulunamadi: %s', gnss_file); end

        opts_imu = detectImportOptions(imu_file, 'VariableNamingRule', 'preserve');
        raw_imu  = readtable(imu_file, opts_imu);
        opts_gnss = detectImportOptions(gnss_file, 'VariableNamingRule', 'preserve');
        raw_gnss = readtable(gnss_file, opts_gnss);

        timeCol_imu  = findTimeCol(raw_imu);
        t_imu  = double(raw_imu.(timeCol_imu));
        timeCol_gnss = findTimeCol(raw_gnss);
        t_gnss = double(raw_gnss.(timeCol_gnss));

        T0 = min(t_imu(1), t_gnss(1));
        ds.time = t_imu - T0;

        ds.ax = findColumnData(raw_imu, {'acc_x', 'ax'});
        ds.ay = findColumnData(raw_imu, {'acc_y', 'ay'});
        ds.az = findColumnData(raw_imu, {'acc_z', 'az'});
        ds.gx = findColumnData(raw_imu, {'gyr_x', 'gx'});
        ds.gy = findColumnData(raw_imu, {'gyr_y', 'gy'});
        ds.gz = findColumnData(raw_imu, {'gyr_z', 'gz'});
        ds.gyr_z = ds.gz;
        ds.gravity_present = true;            % ham IMU: az ~ 9.8

        ds.gps.time = t_gnss - T0;
        ds.gps.lat  = findColumnData(raw_gnss, {'lat', 'latitude'});
        ds.gps.lon  = findColumnData(raw_gnss, {'lon', 'longitude'});
        ds.gps.h    = findColumnData(raw_gnss, {'alt', 'altitude', 'h'});
        ds.gps.hacc = findColumnData(raw_gnss, {'hAcc', 'hacc'});
        ds.gps.heading = findColumnData(raw_gnss, {'heading', 'bearing'});

        % Doppler hiz (modena clean_gnss: speed, velN, velE)
        ds.gps.speed = findColumnData(raw_gnss, {'speed'});
        ds.gps.vN = findColumnData(raw_gnss, {'velN'});
        ds.gps.vE = findColumnData(raw_gnss, {'velE'});
        ds.gps.speedAcc = findColumnData(raw_gnss, {'speedAcc'});
    end

    % hAcc guvenligi
    if all(ds.gps.hacc == 0) || all(isnan(ds.gps.hacc))
        ds.gps.hacc = ones(numel(ds.gps.time), 1) * 10;
    end

    % Ornekleme periyotlari
    if numel(ds.time) > 1
        d = diff(ds.time);
        ds.dt   = mean(d);
        ds.dt_k = [d; d(end)];            % son adim icin son dt tekrar
    else
        ds.dt = 0.01;  ds.dt_k = 0.01;
    end

    ds.acc = [ds.ax, ds.ay, ds.az];
    ds.datetime     = datetime('now') + seconds(ds.time);
    ds.gps.datetime = datetime('now') + seconds(ds.gps.time);
end

% ---------------- yardimci fonksiyonlar ----------------
function colName = findTimeCol(tbl)
    cols = tbl.Properties.VariableNames;
    idx = strcmpi(cols, 't') | strcmpi(cols, 'time') | strcmpi(cols, 'seconds_elapsed');
    if ~any(idx)
        idx = contains(cols, 'time', 'IgnoreCase', true);
    end
    if any(idx)
        colName = cols{find(idx, 1)};
    else
        error('Zaman sutunu bulunamadi (t/time/seconds_elapsed).');
    end
end

function colData = findColumnData(tbl, keywords)
    cols = tbl.Properties.VariableNames;
    colData = zeros(height(tbl), 1);         % bulunamazsa 0
    for i = 1:numel(keywords)
        idx = strcmpi(cols, keywords{i});    % once tam eslesme
        if any(idx)
            colData = tbl.(cols{find(idx,1)});  return;
        end
    end
    for i = 1:numel(keywords)
        idx = contains(cols, keywords{i}, 'IgnoreCase', true);
        if any(idx)
            colData = tbl.(cols{find(idx,1)});  return;
        end
    end
end
