function T = computeMetrics(ds, results, names)
% COMPUTEMETRICS  Filtre ciktilarinin DOGRULUGUNU GPS referansina gore olcer.
%
% Kullanim:
%   T = computeMetrics(ds, {KF_out, EKF_out, fusion_out}, {'KF','EKF','Fusion'});
%
% Metrikler (GPS ornek zamanlarinda, yalnizca iyi GPS hAcc<50 m):
%   PosRMSE_m   : konum RMSE (m)        -> dusuk = iyi
%   PosMax_m    : en buyuk konum hatasi (m)
%   CEP50_m     : hatalarin medyani (m)
%   CEP95_m     : %95 hata (m)
%   FinalDrift_m: son noktadaki konum hatasi (m)  -> DR drift'i icin onemli
%   VelRMSE_ms  : hiz RMSE (m/s) GPS speed'e gore (filtre hiz veriyorsa)
%
% NOT: Burada "referans" GPS'in kendisidir; dolayisiyla GPS'i de guncelleyen
% filtreler (KF/EKF/fusion) GPS'e yakin cikar. Saf-DR (imuLoc) ise gercek
% drift'i gosterir. Bu, mutlak dogruluk degil GORECELI karsilastirmadir;
% mutlak dogruluk icin bagimsiz bir referans (orn. RTK) gerekir.

    Re = 6378137; d2r = pi/180;
    ref_lat = ds.gps.lat(1);  ref_lon = ds.gps.lon(1);

    % GPS -> NED referans
    gN = (ds.gps.lat - ref_lat) * d2r * Re;
    gE = (ds.gps.lon - ref_lon) * d2r * Re * cos(ref_lat*d2r);
    tg = ds.gps.time;

    % iyi GPS maskesi
    good = ds.gps.hacc < 50 & ~isnan(gN) & ~isnan(gE);

    nF = numel(results);
    PosRMSE_m   = zeros(nF,1);
    PosMax_m    = zeros(nF,1);
    CEP50_m     = zeros(nF,1);
    CEP95_m     = zeros(nF,1);
    FinalDrift_m= zeros(nF,1);
    VelRMSE_ms  = nan(nF,1);

    for i = 1:nF
        o = results{i};
        % filtre konumunu GPS zamanlarina interpole et (ayni NED referansi)
        fn = interp1(o.time, o.n, tg, 'linear', NaN);
        fe = interp1(o.time, o.e, tg, 'linear', NaN);
        m  = good & ~isnan(fn) & ~isnan(fe);
        err = sqrt((fn(m)-gN(m)).^2 + (fe(m)-gE(m)).^2);

        if ~isempty(err)
            PosRMSE_m(i)  = sqrt(mean(err.^2));
            PosMax_m(i)   = max(err);
            CEP50_m(i)    = median(err);
            CEP95_m(i)    = prctile_local(err, 95);
            FinalDrift_m(i)= err(end);
        end

        % --- hiz ---
        fv = [];
        if isfield(o,'v') && numel(o.v)==numel(o.time)
            fv = o.v;
        elseif isfield(o,'vn') && isfield(o,'ve')
            fv = sqrt(o.vn.^2 + o.ve.^2);
        end
        if ~isempty(fv) && isfield(ds.gps,'speed')
            fvg = interp1(o.time, fv, tg, 'linear', NaN);
            sp  = ds.gps.speed;
            mv  = good & ~isnan(fvg) & ~isnan(sp) & sp>=0;
            if any(mv)
                VelRMSE_ms(i) = sqrt(mean((fvg(mv)-sp(mv)).^2));
            end
        end
    end

    T = table(names(:), PosRMSE_m, CEP50_m, CEP95_m, PosMax_m, FinalDrift_m, VelRMSE_ms, ...
        'VariableNames', {'Filtre','PosRMSE_m','CEP50_m','CEP95_m','PosMax_m','FinalDrift_m','VelRMSE_ms'});
    T = sortrows(T, 'PosRMSE_m');

    fprintf('\n=== DOGRULUK METRIKLERI (GPS referans, hAcc<50m) ===\n');
    disp(T);
end

function p = prctile_local(x, q)
    % Statistics Toolbox'siz yuzdelik
    x = sort(x(:));
    n = numel(x);
    if n == 1, p = x(1); return; end
    r = (q/100)*(n-1) + 1;
    lo = floor(r);  hi = ceil(r);
    p = x(lo) + (r-lo)*(x(hi)-x(lo));
end
