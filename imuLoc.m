function imu_out = imuLoc(ds)
% IMULOC  Saf ataletsel seyrusefer (dead-reckoning) - SADECE ivme + jiroskop.
%
% AMAC: "Ivmeden tek basina ne kadar konum/hiz cikar?" sorusunu gostermek.
% GPS yalniz BASLANGIC konumu/hizi/yonelimi icin kullanilir; sonra hic GPS yok.
%
% UYARI: Bu yontem KACINILMAZ olarak surulenir (drift). Ivmeolcerdeki kucuk
% bias cift entegrasyonla konumda ~ 0.5*bias*t^2 hata yaratir. Birkac on saniye
% iyi, sonra metrelerce/yuzlerce metre kayar. "Dogru" sonuc icin fusionEKF
% kullanin. Bu fonksiyon yalnizca KARSILASTIRMA/EGITIM amaclidir.
%
% Iyilestirmeler (ham entegrasyona gore):
%   - Bias cikarma: ilk duragan pencereden ivme ortalamasi bias kabul edilir
%   - Alcak geciren filtre (5 Hz) ile gurultu bastirma
%   - ZUPT: yatay ivme ve donus cok kucukse hiz sifirlanir (drift sinirlama)

    N = numel(ds.time);
    imu_out.time = ds.time;
    imu_out.n=zeros(N,1); imu_out.e=zeros(N,1);
    imu_out.vn=zeros(N,1); imu_out.ve=zeros(N,1);
    imu_out.v=zeros(N,1);  imu_out.psi=zeros(N,1);
    imu_out.lat=zeros(N,1); imu_out.lon=zeros(N,1);

    ref_lat = ds.gps.lat(1);
    ref_lon = ds.gps.lon(1);
    dt = ds.dt;

    % --- Alcak geciren (2. derece Butterworth, fc=5 Hz) filtre katsayilari ---
    fs = 1/dt;  fc = 5;
    g = tan(pi*fc/fs);  D = g^2 + sqrt(2)*g + 1;
    b = [g^2/D, 2*g^2/D, g^2/D];
    a = [1, 2*(g^2-1)/D, (g^2 - sqrt(2)*g + 1)/D];
    ax_f = filt2(ds.ax, b, a);
    ay_f = filt2(ds.ay, b, a);

    % --- Bias kestirimi: ilk 2 sn (veya duragan goründügü) pencereden ortalama ---
    win = min(N, max(50, round(2/dt)));
    bax = mean(ax_f(1:win));
    bay = mean(ay_f(1:win));

    % --- Baslangic yonelimi ve hizi (GPS'ten) ---
    psi = init_heading(ds);
    vN = 0; vE = 0;
    if isfield(ds.gps,'vN') && ~isnan(ds.gps.vN(1))
        vN = ds.gps.vN(1);  vE = ds.gps.vE(1);
    end
    pN = 0; pE = 0;

    a_zupt = 0.25;   % m/s^2
    for k = 1:N
        if k > 1, psi = psi + ds.gz(k)*dt; end

        ax = ax_f(k) - bax;     % bias cikar
        ay = ay_f(k) - bay;

        % ZUPT: duruyorsa hizi sifirla
        if sqrt(ax^2+ay^2) < a_zupt && abs(ds.gz(k)) < 0.05
            vN = 0; vE = 0;
        else
            aN = ax*cos(psi) - ay*sin(psi);
            aE = ax*sin(psi) + ay*cos(psi);
            pN = pN + vN*dt + 0.5*aN*dt^2;
            pE = pE + vE*dt + 0.5*aE*dt^2;
            vN = vN + aN*dt;
            vE = vE + aE*dt;
        end

        imu_out.n(k)=pN; imu_out.e(k)=pE;
        imu_out.vn(k)=vN; imu_out.ve(k)=vE;
        imu_out.v(k)=sqrt(vN^2+vE^2); imu_out.psi(k)=psi;
        [imu_out.lat(k), imu_out.lon(k)] = ned2lla(pN, pE, ref_lat, ref_lon);
    end
end

% ---------------- yardimci ----------------
function y = filt2(x, b, a)
    % basit IIR (filter benzeri), Signal Toolbox gerektirmez
    N = numel(x);  y = zeros(N,1);
    y(1) = x(1);  if N>1, y(2) = x(2); end
    for i = 3:N
        y(i) = b(1)*x(i) + b(2)*x(i-1) + b(3)*x(i-2) - a(2)*y(i-1) - a(3)*y(i-2);
    end
end

function psi0 = init_heading(ds)
    psi0 = 0;
    if isfield(ds.gps,'heading')
        h = ds.gps.heading(~isnan(ds.gps.heading));
        if ~isempty(h), psi0 = h(1)*pi/180; return; end
    end
    if isfield(ds.gps,'vN')
        for i = 1:numel(ds.gps.vN)
            if ~isnan(ds.gps.vN(i)) && hypot(ds.gps.vN(i),ds.gps.vE(i)) > 2
                psi0 = atan2(ds.gps.vE(i), ds.gps.vN(i)); return;
            end
        end
    end
end

function [lat, lon] = ned2lla(n, e, lat0, lon0)
    Re = 6378137; r2d = 180/pi;
    lat = lat0 + (n / Re) * r2d;
    lon = lon0 + (e / (Re * cos(lat0*pi/180))) * r2d;
end
