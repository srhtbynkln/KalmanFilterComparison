# -*- coding: utf-8 -*-
# Saf IMU (dead-reckoning) vs Saf GPS vs Fuzyon (EKF) karsilastirmasi.
# Hem telefon (Accelerometer/Gyroscope/Location/Orientation) hem marmaray formatini okur.
# Cikti: <dir>_compare.png  +  konsolda dogruluk metrikleri tablosu.
# Kullanim: python3 compare_methods.py data_m34/session_YYYYMMDD_HHMMSS
import sys, os, csv, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = (sys.argv[1] if len(sys.argv) > 1 else "data_marmaray").rstrip("/")
def load(n):
    p = os.path.join(d, n)
    return list(csv.DictReader(open(p, newline=""))) if os.path.exists(p) else []
def c(rows, k, default=None):
    out = []
    for r in rows:
        try: out.append(float(r[k]))
        except: out.append(default if default is not None else np.nan)
    return np.array(out)

acc = load("Accelerometer.csv"); gyr = load("Gyroscope.csv"); loc = load("Location.csv")
ori = load("Orientation.csv")
if not acc or not loc:
    print("Accelerometer/Location yok:", d); sys.exit(1)

t = c(acc, "seconds_elapsed"); ax = c(acc, "x"); ay = c(acc, "y"); az = c(acc, "z")
N = len(t); dt_k = np.append(np.diff(t), np.diff(t)[-1] if N > 1 else 0.01)
tgz = c(gyr, "seconds_elapsed"); gz = np.interp(t, tgz, c(gyr, "z")) if gyr else np.zeros(N)

HAS_ORI = bool(ori) and all(k in ori[0] for k in ("qw", "qx", "qy", "qz"))
if HAS_ORI:
    to = c(ori, "seconds_elapsed")
    qw = np.interp(t, to, c(ori, "qw")); qx = np.interp(t, to, c(ori, "qx"))
    qy = np.interp(t, to, c(ori, "qy")); qz = np.interp(t, to, c(ori, "qz"))

tg = c(loc, "seconds_elapsed"); lat = c(loc, "latitude"); lon = c(loc, "longitude")
hacc = c(loc, "horizontalAccuracy")
spd = c(loc, "speed"); spd[spd < 0] = np.nan
brg = c(loc, "bearing"); brg[brg < 0] = np.nan
vN_g = spd * np.cos(np.radians(brg)); vE_g = spd * np.sin(np.radians(brg))
Re = 6378137.0; d2r = math.pi / 180; rl = lat[0]
gN = (lat - lat[0]) * d2r * Re; gE = (lon - lon[0]) * d2r * Re * math.cos(rl * d2r)
ngps = len(tg)
HACC_OK = 50.0

def aNE(k, psi):
    if HAS_ORI:
        w, x, y, z = qw[k], qx[k], qy[k], qz[k]
        ux, uy, uz = ax[k], ay[k], az[k]
        aE = (1-2*(y*y+z*z))*ux + 2*(x*y-z*w)*uy + 2*(x*z+y*w)*uz
        aN = 2*(x*y+z*w)*ux + (1-2*(x*x+z*z))*uy + 2*(y*z-x*w)*uz
        return aN, aE
    c_, s_ = math.cos(psi), math.sin(psi)
    return c_*ax[k]-s_*ay[k], s_*ax[k]+c_*ay[k]

def init_heading():
    for i in range(ngps):
        if not np.isnan(brg[i]): return math.radians(brg[i])
    return 0.0

# ---------------- 1) SAF IMU (dead-reckoning) ----------------
def imu_only():
    psi = init_heading(); on = np.zeros(N); oe = np.zeros(N)
    vN = vE = pN = pE = 0.0
    if not np.isnan(vN_g[0]): vN, vE = vN_g[0], vE_g[0]
    for k in range(N):
        if k > 0 and not HAS_ORI: psi += gz[k]*dt_k[k]
        aN, aE = aNE(k, psi)
        if math.hypot(ax[k], ay[k]) < 0.3 and abs(gz[k]) < 0.05:
            vN = vE = 0.0
        else:
            pN += vN*dt_k[k] + 0.5*aN*dt_k[k]**2; pE += vE*dt_k[k] + 0.5*aE*dt_k[k]**2
            vN += aN*dt_k[k]; vE += aE*dt_k[k]
        on[k] = pN; oe[k] = pE
    return on, oe

# ---------------- 3) FUZYON (EKF: GPS pos+vel + IMU + bias) ----------------
def ekf_upd(x, P, z, H, R):
    y = z - H@x; S = H@P@H.T + R
    if np.linalg.cond(S) > 1e14: return x, P
    K = P@H.T@np.linalg.inv(S); x = x + K@y
    I = np.eye(6); P = (I-K@H)@P@(I-K@H).T + K@R@K.T
    return x, 0.5*(P+P.T)

def fusion():
    gi = [i for i in range(ngps) if 0 < hacc[i] < HACC_OK and not np.isnan(gN[i])]
    fg = gi[0] if gi else 0; it = tg[fg]
    x = np.zeros(6); x[0] = gN[fg]; x[1] = gE[fg]
    if not np.isnan(vN_g[fg]): x[2] = vN_g[fg]; x[3] = vE_g[fg]
    P = np.diag([25., 25, 4, 4, 0.25, 0.25]); psi = init_heading()
    sa, sb, Kpsi = 0.5, 0.002, 0.02
    Hc = np.array([[1,0,0,0,0,0],[0,1,0,0,0,0.]]); Hv = np.array([[0,0,1,0,0,0],[0,0,0,1,0,0.]])
    fn = np.zeros(N); fe = np.zeros(N); gp = 0; inited = False; lastg = -1e18
    for k in range(N):
        dt = dt_k[k]
        if not inited:
            if t[k] < it:
                fn[k] = x[0]; fe[k] = x[1]
                while gp < ngps and tg[gp] <= t[k]+1e-9: gp += 1
                continue
            inited = True; lastg = it
        outage = (t[k]-lastg) > 1.5
        if k > 0 and not HAS_ORI: psi += gz[k]*dt
        aN, aE = aNE(k, psi)
        A = np.array([[1,0,dt,0,-0.5*dt*dt,0],[0,1,0,dt,0,-0.5*dt*dt],
                      [0,0,1,0,-dt,0],[0,0,0,1,0,-dt],[0,0,0,0,1,0],[0,0,0,0,0,1.]])
        x = A@x + np.array([0.5*dt*dt*aN,0.5*dt*dt*aE,dt*aN,dt*aE,0,0.])
        qp = 0.25*dt**4*sa*sa; qv = dt*dt*sa*sa; qb = dt*sb*sb
        P = A@P@A.T + np.diag([qp,qp,qv,qv,qb,qb])
        while gp < ngps and tg[gp] <= t[k]+1e-9:
            h = hacc[gp]
            if not (0 < h < HACC_OK): gp += 1; continue
            if (t[k]-lastg) > 3.0:
                P[0,0]+=2500; P[1,1]+=2500; P[2,2]+=25; P[3,3]+=25
            x, P = ekf_upd(x, P, np.array([gN[gp], gE[gp]]), Hc, np.diag([h*h, h*h]))
            if not np.isnan(vN_g[gp]):
                x, P = ekf_upd(x, P, np.array([vN_g[gp], vE_g[gp]]), Hv, np.diag([.09, .09]))
                if not HAS_ORI and not np.isnan(spd[gp]) and spd[gp] > 2:
                    pg = math.atan2(vE_g[gp], vN_g[gp]); psi += Kpsi*math.atan2(math.sin(pg-psi), math.cos(pg-psi))
            lastg = t[k]; gp += 1
        if not outage and math.hypot(ax[k], ay[k]) < 0.3 and abs(gz[k]) < 0.05 and math.hypot(x[2], x[3]) < 0.4:
            x, P = ekf_upd(x, P, np.array([0, 0.]), Hv, np.diag([.0025, .0025]))
        fn[k] = x[0]; fe[k] = x[1]
    return fn, fe

iN, iE = imu_only()
fN, fE = fusion()

# ---------------- METRIKLER (referans = iyi GPS noktalari) ----------------
good = (hacc < HACC_OK) & ~np.isnan(gN)
if not good.any():
    print(f"Gecerli GPS noktasi yok (hAcc<{HACC_OK:.0f} m kosulunu saglayan fix yok): {d}")
    print("Metrik uretilemez; daha temiz GPS'li bir oturum gerekir."); sys.exit(0)
def metrics(name, on, oe):
    pn = np.interp(tg, t, on); pe = np.interp(tg, t, oe)
    err = np.hypot(pn[good]-gN[good], pe[good]-gE[good])
    return name, np.sqrt(np.mean(err**2)), np.median(err), err.max(), err[-1]

rows = [metrics("Saf IMU", iN, iE), metrics("Fuzyon (EKF)", fN, fE)]
print("="*72)
print(f"VERI: {d}   (GPS nokta: {int(good.sum())}, hAcc med: {np.median(hacc[good]):.0f} m, "
      f"orientation: {'VAR' if HAS_ORI else 'YOK (yaw-only)'})")
print("-"*72)
print(f"{'Yontem':16s} {'PosRMSE':>9s} {'CEP50':>8s} {'Maks':>9s} {'SonHata':>9s}")
for n, r, c50, mx, fin in rows:
    print(f"{n:16s} {r:8.1f}m {c50:7.1f}m {mx:8.1f}m {fin:8.1f}m")
print("(referans: GPS'in kendisi; saf GPS hatasi tanim geregi 0)")
# Basari orani = hata esigin altinda kalan GPS noktalarinin yuzdesi (availability)
pnF = np.interp(tg, t, fN); peF = np.interp(tg, t, fE)
ef = np.hypot(pnF[good]-gN[good], peF[good]-gE[good])
pnI = np.interp(tg, t, iN); peI = np.interp(tg, t, iE)
ei = np.hypot(pnI[good]-gN[good], peI[good]-gE[good])
print("-"*72)
print(f"{'Basari orani (hata<esik %)':28s} {'<5m':>7s} {'<10m':>7s} {'<20m':>7s}")
print(f"{'Saf IMU':28s} {100*np.mean(ei<5):6.0f}% {100*np.mean(ei<10):6.0f}% {100*np.mean(ei<20):6.0f}%")
print(f"{'Fuzyon (EKF)':28s} {100*np.mean(ef<5):6.0f}% {100*np.mean(ef<10):6.0f}% {100*np.mean(ef<20):6.0f}%")

# ---------------- GRAFIK ----------------
fig, axp = plt.subplots(1, 2, figsize=(15, 7))
a0 = axp[0]
a0.plot(gE[good], gN[good], "o-", color="#2ca02c", ms=4, lw=1, label="Saf GPS (referans)")
a0.plot(iE, iN, "-", color="#ff7f0e", lw=1.5, label="Saf IMU (dead-reckoning)")
a0.plot(fE, fN, "-", color="#1f77b4", lw=2, label="Füzyon (EKF)")
a0.scatter(gE[good][0], gN[good][0], c="k", s=80, marker="o", zorder=5, label="Başlangıç")
a0.set_aspect("equal", "datalim"); a0.grid(alpha=.3); a0.legend()
a0.set_xlabel("Doğu (m)"); a0.set_ylabel("Kuzey (m)")
a0.set_title(f"Yörünge karşılaştırması — {os.path.basename(d)}")

a1 = axp[1]
pnI = np.interp(tg, t, iN); peI = np.interp(tg, t, iE)
pnF = np.interp(tg, t, fN); peF = np.interp(tg, t, fE)
a1.plot(tg[good], np.hypot(pnI[good]-gN[good], peI[good]-gE[good]), color="#ff7f0e", label="Saf IMU hata")
a1.plot(tg[good], np.hypot(pnF[good]-gN[good], peF[good]-gE[good]), color="#1f77b4", label="Füzyon hata")
a1.set_xlabel("t (s)"); a1.set_ylabel("GPS'e uzaklık (m)"); a1.grid(alpha=.3); a1.legend()
a1.set_title("Konum hatası (GPS referansına göre)")

out = d + "_compare.png"
plt.tight_layout(); plt.savefig(out, dpi=130); print("Grafik:", out)
