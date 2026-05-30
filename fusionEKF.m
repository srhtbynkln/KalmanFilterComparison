function out = fusionEKF(ds)
% FUSIONEKF  Gevsek-bagli (loosely-coupled) GNSS/INS fuzyonu.
%
% AMAC: Ivmeolcerden HIZ ve KONUM bilgisini DOGRU ve SINIRLI hatayla cikarmak.
% Bunu tek basina ivme ile yapmak imkansizdir (bias cift entegrasyonla t^2
% buyur). Cozum: GPS konumu + GPS Doppler HIZI ile ivmeyi birlikte fuzyonlamak.
% Ivme, GPS ornekleri ARASINI doldurur ve tunel/kesinti aninda kestirimi surdurur.
%
% Durum vektoru (6x1):
%   x = [ N ; E ; vN ; vE ; bax ; bay ]
%     N,E    : konum (m, yerel NED)
%     vN,vE  : hiz (m/s, Kuzey/Dogu)
%     bax,bay: ivmeolcer bias'i (m/s^2, GOVDE ekseninde)
%
% Yonelim (psi): GPS rotasi + jiroskop ile tamamlayici (complementary) filtre
% ile durum DISINDA tutulur (GPS hizi psi'yi cok iyi gozledigi icin yeterli).
%
% Olcumler: GPS konum (N,E) + GPS hiz (vN,vE).  Kotu GPS chi-square ile reddedilir.
%
% NOT (sinirlama): Veride telefon Orientation (rotation vector) olmadigindan
% govde->NED donusumu icin yalniz yaw (psi) kullanilir; telefonun yatay
% durdugu varsayilir. Orientation loglanirsa bu varsayim kalkar (bkz. REHBER.md).

    N = numel(ds.time);
    out.time = ds.time;
    out.n  = zeros(N,1);  out.e  = zeros(N,1);
    out.vn = zeros(N,1);  out.ve = zeros(N,1);
    out.v  = zeros(N,1);  out.psi = zeros(N,1);
    out.bax = zeros(N,1); out.bay = zeros(N,1);
    out.lat = zeros(N,1); out.lon = zeros(N,1);

    ref_lat = ds.gps.lat(1);
    ref_lon = ds.gps.lon(1);

    % GPS'i bir kez NED'e cevir
    num_gps = numel(ds.gps.time);
    gps_n = zeros(num_gps,1);  gps_e = zeros(num_gps,1);
    for i = 1:num_gps
        [gps_n(i), gps_e(i)] = lla2ned(ds.gps.lat(i), ds.gps.lon(i), ref_lat, ref_lon);
    end

    % ---- Baslangic durumu ----
    x = zeros(6,1);
    x(1) = gps_n(1);  x(2) = gps_e(1);
    % ilk gecerli GPS hizi varsa onu kullan
    if isfield(ds.gps,'vN') && ~isnan(ds.gps.vN(1))
        x(3) = ds.gps.vN(1);  x(4) = ds.gps.vE(1);
    end
    P = diag([5, 5, 2, 2, 0.5, 0.5].^2);

    % Baslangic yonelimi: ilk gecerli GPS heading'i, yoksa hiz vektoru yonu
    psi = init_heading(ds);

    % ---- Gurultu parametreleri ----
    sigma_a   = 0.5;     % ivme islem gurultusu (m/s^2) - telefon icin
    sigma_bias= 0.002;   % bias rastgele yuruyusu (m/s^2/sqrt(s))
    Kpsi      = 0.02;    % yonelim tamamlayici filtre kazanci (GPS rotasi)
    v_zupt    = 0.4;     % bu hizin altinda ZUPT denenir (m/s)
    a_zupt    = 0.3;     % yatay ivme bu altindaysa duruyor say (m/s^2)

    gps_idx = 1;
    I6 = eye(6);

    for k = 1:N
        dt = ds.dt_k(k);
        c = cos(psi);  s = sin(psi);

        % --- Yonelim tahmini: jiroskop entegrasyonu ---
        psi = psi + ds.gz(k) * dt;

        % --- Olculen govde ivmesi (yatay) ---
        ax = ds.ax(k);  ay = ds.ay(k);

        % --- Tahmin (Prediction) ---
        % aN = c*ax - s*ay - (c*bax - s*bay)
        % aE = s*ax + c*ay - (s*bax + c*bay)
        % Durum gecisi A (bias->hiz/konum eslemesi psi'ye bagli, zaman-degisken)
        A = [ 1 0 dt 0  -0.5*dt^2*c   0.5*dt^2*s ;
              0 1 0 dt  -0.5*dt^2*s  -0.5*dt^2*c ;
              0 0 1 0      -dt*c         dt*s     ;
              0 0 0 1      -dt*s        -dt*c     ;
              0 0 0 0       1            0        ;
              0 0 0 0       0            1        ];
        % Olculen ivmenin girisi (bias'siz kisim)
        aN_meas = c*ax - s*ay;
        aE_meas = s*ax + c*ay;
        Bu = [ 0.5*dt^2*aN_meas ;
               0.5*dt^2*aE_meas ;
                   dt*aN_meas    ;
                   dt*aE_meas    ;
                   0 ; 0 ];
        x = A*x + Bu;

        % Surec gurultusu Q
        q_pos = 0.25*dt^4*sigma_a^2;
        q_vel = dt^2*sigma_a^2;
        q_b   = dt*sigma_bias^2;
        Q = diag([q_pos, q_pos, q_vel, q_vel, q_b, q_b]);
        P = A*P*A' + Q;

        % --- GPS guncelleme(leri) ---
        while gps_idx <= num_gps && ds.gps.time(gps_idx) <= ds.time(k) + 1e-9
            hacc = ds.gps.hacc(gps_idx);
            % Kotu GPS'i tamamen atla (tunel / multipath)
            if ~(hacc > 0 && hacc < 50)
                gps_idx = gps_idx + 1;  continue;
            end

            % --- Konum olcumu (her zaman) ---
            zc = [gps_n(gps_idx); gps_e(gps_idx)];
            Hc = [1 0 0 0 0 0; 0 1 0 0 0 0];
            Rc = diag([hacc^2, hacc^2]);
            % NOT: outlier savunmasi hAcc<50 sert filtresiyle yapilir. Ayrica
            % chi-square kapisi EKLENMEZ; cunku erken bir sapmada kapi GPS'i
            % reddedip filtreyi kalici iraksamaya sokuyordu (test edildi).
            [x, P] = ekf_update(x, P, zc, Hc, Rc, I6, Inf);

            % --- Hiz olcumu (Doppler varsa) ---
            if isfield(ds.gps,'vN') && ~isnan(ds.gps.vN(gps_idx))
                zv = [ds.gps.vN(gps_idx); ds.gps.vE(gps_idx)];
                Hv = [0 0 1 0 0 0; 0 0 0 1 0 0];
                sv = 0.3;  % varsayilan hiz olcum std (m/s)
                if isfield(ds.gps,'speedAcc') && ~isnan(ds.gps.speedAcc(gps_idx)) ...
                        && ds.gps.speedAcc(gps_idx) > 0
                    sv = ds.gps.speedAcc(gps_idx);
                end
                Rv = diag([sv^2, sv^2]);
                [x, P] = ekf_update(x, P, zv, Hv, Rv, I6, Inf);

                % Yonelimi GPS rotasi ile duzelt (yeterince hizliysa)
                spd = ds.gps.speed(gps_idx);
                if ~isnan(spd) && spd > 2.0
                    psi_gps = atan2(ds.gps.vE(gps_idx), ds.gps.vN(gps_idx));
                    psi = psi + Kpsi * wrapToPiLocal(psi_gps - psi);
                end
            end
            gps_idx = gps_idx + 1;
        end

        % --- ZUPT: arac duruyorsa hizi sifira cek (bias'i da gozle) ---
        a_horz = sqrt(ax^2 + ay^2);
        spd_est = sqrt(x(3)^2 + x(4)^2);
        if a_horz < a_zupt && abs(ds.gz(k)) < 0.05 && spd_est < v_zupt
            Hz = [0 0 1 0 0 0; 0 0 0 1 0 0];
            [x, P] = ekf_update(x, P, [0;0], Hz, diag([0.05^2,0.05^2]), I6, Inf);
        end

        % --- Kayit ---
        out.n(k)=x(1); out.e(k)=x(2); out.vn(k)=x(3); out.ve(k)=x(4);
        out.bax(k)=x(5); out.bay(k)=x(6);
        out.v(k)=sqrt(x(3)^2+x(4)^2);
        out.psi(k)=psi;
        [out.lat(k), out.lon(k)] = ned2lla(x(1), x(2), ref_lat, ref_lon);
    end
end

% ================= yardimci fonksiyonlar =================
function [x, P] = ekf_update(x, P, z, H, R, I, gate)
    y = z - H*x;
    S = H*P*H' + R;
    if rcond(S) < 1e-15, return; end
    d2 = y' / S * y;                 % Mahalanobis^2
    if d2 > gate, return; end        % outlier -> guncelleme yok
    K = (P*H') / S;
    x = x + K*y;
    P = (I - K*H)*P*(I - K*H)' + K*R*K';   % Joseph formu
    P = 0.5*(P + P');
end

function psi0 = init_heading(ds)
    psi0 = 0;
    if isfield(ds.gps,'heading')
        h = ds.gps.heading(~isnan(ds.gps.heading));
        if ~isempty(h), psi0 = h(1) * pi/180;  return; end
    end
    if isfield(ds.gps,'vN')
        for i = 1:numel(ds.gps.vN)
            if ~isnan(ds.gps.vN(i)) && hypot(ds.gps.vN(i),ds.gps.vE(i)) > 2
                psi0 = atan2(ds.gps.vE(i), ds.gps.vN(i));  return;
            end
        end
    end
end

function a = wrapToPiLocal(a)
    a = atan2(sin(a), cos(a));
end

function [n, e] = lla2ned(lat, lon, lat0, lon0)
    Re = 6378137; d2r = pi/180;
    n = (lat - lat0) * d2r * Re;
    e = (lon - lon0) * d2r * Re * cos(lat0*d2r);
end

function [lat, lon] = ned2lla(n, e, lat0, lon0)
    Re = 6378137; r2d = 180/pi;
    lat = lat0 + (n / Re) * r2d;
    lon = lon0 + (e / (Re * cos(lat0*pi/180))) * r2d;
end
