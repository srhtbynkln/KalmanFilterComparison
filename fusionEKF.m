function out = fusionEKF(ds)
% FUSIONEKF  Gevsek-bagli (loosely-coupled) GNSS/INS fuzyonu.
%
% AMAC: Ivmeolcerden HIZ ve KONUM bilgisini DOGRU ve SINIRLI hatayla cikarmak.
% Tek basina ivme ile imkansiz (bias cift entegrasyonla t^2 buyur). Cozum:
% GPS konumu + GPS Doppler HIZI ile ivmeyi birlikte fuzyonlamak. Ivme, GPS
% ornekleri arasini doldurur ve tunel/kesintide kestirimi surdurur.
%
% Durum vektoru (6x1):  x = [ N ; E ; vN ; vE ; bN ; bE ]
%   N,E   : konum (m, yerel NED)
%   vN,vE : hiz (m/s)
%   bN,bE : ivme bias'i (m/s^2, NAV cercevesinde -> A matrisi SABIT olur)
%
% IVMENIN NAV CERCEVESINE DONUSU - iki mod (otomatik secilir):
%   1) ds.has_orientation = true  -> TAM 3B rotasyon (quaternion qw,qx,qy,qz).
%      Telefon egik/serbest olsa bile dogru. Ayrica yercekimi dunya
%      cercevesinde yalniz dikey eksende oldugundan yatay (N,E) ivmeyi
%      ETKILEMEZ; gravity sorunu otomatik cozulur.
%   2) Orientation yok -> YAW-ONLY fallback (eski davranis): yalniz jiroskop
%      yaw'i + GPS rotasi; telefon yatay varsayilir.
%
% KONVANSIYON NOTU (orientation modu): quaternion'in govde->DUNYA (ENU:
% X=Dogu, Y=Kuzey, Z=Yukari) rotasyonu oldugu varsayilir (Android rotation
% vector standardi). iOS/farkli export'ta eksen sirasi/isareti degisebilir;
% gerekirse quatRotate cikisindaki aE/aN eslemesini guncelleyin.

    N = numel(ds.time);
    out.time = ds.time;
    out.n  = zeros(N,1);  out.e  = zeros(N,1);
    out.vn = zeros(N,1);  out.ve = zeros(N,1);
    out.v  = zeros(N,1);  out.psi = zeros(N,1);
    out.bN = zeros(N,1);  out.bE = zeros(N,1);
    out.lat = zeros(N,1); out.lon = zeros(N,1);

    ref_lat = ds.gps.lat(1);
    ref_lon = ds.gps.lon(1);

    % GPS'i bir kez NED'e cevir
    num_gps = numel(ds.gps.time);
    gps_n = zeros(num_gps,1);  gps_e = zeros(num_gps,1);
    for i = 1:num_gps
        [gps_n(i), gps_e(i)] = lla2ned(ds.gps.lat(i), ds.gps.lon(i), ref_lat, ref_lon);
    end

    % ---- Baslangic durumu (ILK GUVENILIR FIX'E ERTELENIR) ----
    % GPS kaydi cold-start'ta (ve tunelden cikista) birkac saniye COP fix
    % verir: hAcc onlarca-binlerce metre. Filtreyi bu coplerden birinde
    % KUCUK P ile baslatmak, yanlis konuma "emin" baslamaya ve donen iyi
    % GPS'i zayif kazancla alip onlarca saniye yavas yakinsamaya yol acar
    % (marmaray'da tam bu oluyordu: ilk iyi fix'te 336 m hata).
    % COZUM: ilk hAcc<50 fix'ini bul, durumu ORADA baslat; o ana kadar
    % ciktiyi sabit tut, IMU'yu entegre etme.
    first_good = find(ds.gps.hacc > 0 & ds.gps.hacc < 50 & ~isnan(gps_n), 1);
    if isempty(first_good), first_good = 1; end
    init_time = ds.gps.time(first_good);

    x = zeros(6,1);
    x(1) = gps_n(first_good);  x(2) = gps_e(first_good);
    if isfield(ds.gps,'vN') && ~isnan(ds.gps.vN(first_good))
        x(3) = ds.gps.vN(first_good);  x(4) = ds.gps.vE(first_good);
    end
    P = diag([5, 5, 2, 2, 0.5, 0.5].^2);

    psi = init_heading(ds);

    initialized = false;     % ilk iyi fix'e ulasildi mi
    last_gps_t  = -inf;      % son KABUL edilen iyi GPS zamani (outage tespiti)

    % ---- parametreler ----
    sigma_a   = 0.5;     sigma_bias = 0.002;
    Kpsi      = 0.02;    v_zupt = 0.4;   a_zupt = 0.3;
    outage_T  = 1.5;     % s; bundan uzun GPS yoklugu -> 'outage' (tunel)

    useOri = isfield(ds,'has_orientation') && ds.has_orientation;
    if useOri, disp('fusionEKF: TAM 3B rotasyon (Orientation) kullaniliyor');
    else,      disp('fusionEKF: yaw-only mod (Orientation yok)'); end

    gps_idx = 1;
    I6 = eye(6);
    Hc = [1 0 0 0 0 0; 0 1 0 0 0 0];
    Hv = [0 0 1 0 0 0; 0 0 0 1 0 0];

    for k = 1:N
        dt = ds.dt_k(k);

        % --- Ilk guvenilir fix'e kadar BEKLE (deferred init) ---
        if ~initialized
            if ds.time(k) < init_time
                out.n(k)=x(1); out.e(k)=x(2);
                out.vn(k)=0; out.ve(k)=0; out.v(k)=0;
                out.bN(k)=x(5); out.bE(k)=x(6); out.psi(k)=psi;
                [out.lat(k), out.lon(k)] = ned2lla(x(1), x(2), ref_lat, ref_lon);
                while gps_idx <= num_gps && ds.gps.time(gps_idx) <= ds.time(k)+1e-9
                    gps_idx = gps_idx + 1;        % init oncesi cop fix'leri atla
                end
                continue;
            end
            initialized = true;  last_gps_t = init_time;
        end

        % outage (tunel) tespiti: bu adima girerken son iyi GPS'ten gecen sure
        outage = (ds.time(k) - last_gps_t) > outage_T;

        % --- nav-cerceve ivme (aN_meas, aE_meas) ---
        if useOri
            q  = [ds.qw(k) ds.qx(k) ds.qy(k) ds.qz(k)];
            aw = quatRotate(q, [ds.ax(k); ds.ay(k); ds.az(k)]);  % [aE; aN; aU]
            if ds.gravity_present, aw(3) = aw(3) - 9.81; end       % yatayi etkilemez
            aE_meas = aw(1);  aN_meas = aw(2);
        else
            c = cos(psi);  s = sin(psi);
            psi = psi + ds.gz(k) * dt;
            aN_meas = c*ds.ax(k) - s*ds.ay(k);
            aE_meas = s*ds.ax(k) + c*ds.ay(k);
        end

        % --- Tahmin (A SABIT: bias nav cercevesinde) ---
        A = [ 1 0 dt 0 -0.5*dt^2 0 ;
              0 1 0 dt 0 -0.5*dt^2 ;
              0 0 1 0  -dt      0 ;
              0 0 0 1   0      -dt ;
              0 0 0 0   1       0 ;
              0 0 0 0   0       1 ];
        Bu = [0.5*dt^2*aN_meas; 0.5*dt^2*aE_meas; dt*aN_meas; dt*aE_meas; 0; 0];
        x = A*x + Bu;

        Q = diag([0.25*dt^4*sigma_a^2, 0.25*dt^4*sigma_a^2, ...
                  dt^2*sigma_a^2,      dt^2*sigma_a^2, ...
                  dt*sigma_bias^2,     dt*sigma_bias^2]);
        P = A*P*A' + Q;

        % --- GPS guncelleme(leri) ---
        while gps_idx <= num_gps && ds.gps.time(gps_idx) <= ds.time(k) + 1e-9
            hacc = ds.gps.hacc(gps_idx);
            if ~(hacc > 0 && hacc < 50)         % kotu GPS (tunel/multipath) -> atla
                gps_idx = gps_idx + 1;  continue;
            end

            % Outage'dan (tunel) DONUSTE kovaryansi sis: kucuk P yuzunden donen
            % iyi GPS zayif kazancla alinir ve hata 10s'lerce yavas iner. Buyuk
            % P -> filtre ilk iyi fix'in konumunu hizla yakalar.
            if (ds.time(k) - last_gps_t) > 3.0
                P(1,1)=P(1,1)+2500; P(2,2)=P(2,2)+2500;   % konum (~50 m std)
                P(3,3)=P(3,3)+25;   P(4,4)=P(4,4)+25;      % hiz   (~5 m/s std)
            end

            % konum (outlier savunmasi hAcc filtresi; chi-square kapisi YOK,
            % erken sapmada filtreyi kalici iraksamaya sokuyordu - test edildi)
            zc = [gps_n(gps_idx); gps_e(gps_idx)];
            [x, P] = ekf_update(x, P, zc, Hc, diag([hacc^2, hacc^2]), I6);

            % Doppler hiz
            if isfield(ds.gps,'vN') && ~isnan(ds.gps.vN(gps_idx))
                zv = [ds.gps.vN(gps_idx); ds.gps.vE(gps_idx)];
                sv = 0.3;
                if isfield(ds.gps,'speedAcc') && ~isnan(ds.gps.speedAcc(gps_idx)) ...
                        && ds.gps.speedAcc(gps_idx) > 0
                    sv = ds.gps.speedAcc(gps_idx);
                end
                [x, P] = ekf_update(x, P, zv, Hv, diag([sv^2, sv^2]), I6);

                % yaw-only modda yonelimi GPS rotasi ile duzelt
                if ~useOri
                    spd = ds.gps.speed(gps_idx);
                    if ~isnan(spd) && spd > 2.0
                        psi_gps = atan2(ds.gps.vE(gps_idx), ds.gps.vN(gps_idx));
                        psi = psi + Kpsi * wrapToPiLocal(psi_gps - psi);
                    end
                end
            end
            last_gps_t = ds.time(k);            % iyi fix kabul edildi
            gps_idx = gps_idx + 1;
        end

        % --- ZUPT (arac gercekten duruyorsa) ---
        % Outage sirasinda UYGULANMAZ: GPS yokken durgun gorunen IMU, gercekte
        % hareket eden araci yanlislikla sifir hiza cakabilir (tunelde kritik).
        a_horz  = hypot(ds.ax(k), ds.ay(k));
        spd_est = hypot(x(3), x(4));
        if ~outage && a_horz < a_zupt && abs(ds.gz(k)) < 0.05 && spd_est < v_zupt
            [x, P] = ekf_update(x, P, [0;0], Hv, diag([0.05^2,0.05^2]), I6);
        end

        % --- Kayit ---
        out.n(k)=x(1); out.e(k)=x(2); out.vn(k)=x(3); out.ve(k)=x(4);
        out.bN(k)=x(5); out.bE(k)=x(6);
        out.v(k)=hypot(x(3),x(4));
        if useOri, out.psi(k)=atan2(x(4),x(3)); else, out.psi(k)=psi; end
        [out.lat(k), out.lon(k)] = ned2lla(x(1), x(2), ref_lat, ref_lon);
    end
end

% ================= yardimci fonksiyonlar =================
function [x, P] = ekf_update(x, P, z, H, R, I)
    y = z - H*x;
    S = H*P*H' + R;
    if rcond(S) < 1e-15, return; end
    K = (P*H') / S;
    x = x + K*y;
    P = (I - K*H)*P*(I - K*H)' + K*R*K';   % Joseph formu
    P = 0.5*(P + P');
end

function v = quatRotate(q, u)
    % q = [w x y z], govde->dunya rotasyonu; v = R(q)*u
    w=q(1); x=q(2); y=q(3); z=q(4);
    n = sqrt(w*w+x*x+y*y+z*z);  if n>0, w=w/n; x=x/n; y=y/n; z=z/n; end
    R = [ 1-2*(y^2+z^2),  2*(x*y-z*w),    2*(x*z+y*w) ;
          2*(x*y+z*w),    1-2*(x^2+z^2),  2*(y*z-x*w) ;
          2*(x*z-y*w),    2*(y*z+x*w),    1-2*(x^2+y^2) ];
    v = R*u;
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
