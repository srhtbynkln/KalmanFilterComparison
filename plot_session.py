# -*- coding: utf-8 -*-
# Telefon kaydini (KalmanLogger oturum klasoru) gorsel ozetler:
#   1) GPS track (yorunge haritasi, GPSLogger gibi)
#   2) Hiz - zaman
#   3) Kat edilen mesafe - zaman
#   4) Yonelim: baslangica gore donus (jiro, deadband) + GPS rotasi
# Kullanim: python3 plot_session.py data_m34/session_YYYYMMDD_HHMMSS
import sys, os, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sess = sys.argv[1] if len(sys.argv) > 1 else "data_m34/session_20260530_122331"
sess = sess.rstrip("/")

def load(name):
    p = os.path.join(sess, name)
    if not os.path.exists(p): return []
    with open(p, newline="") as f:
        return list(csv.DictReader(f))
def col(rows, k):
    return np.array([float(r[k]) for r in rows])

loc = load("Location.csv")
gyr = load("Gyroscope.csv")

fig = plt.figure(figsize=(13, 9))
fig.suptitle(f"Kayit ozeti — {os.path.basename(sess)}", fontsize=14, fontweight="bold")

# ---------- GPS turetilmis ----------
if loc:
    t = col(loc, "seconds_elapsed"); t = t - t[0]
    lat = col(loc, "latitude"); lon = col(loc, "longitude")
    hacc = col(loc, "horizontalAccuracy")
    spd = col(loc, "speed"); spd[spd < 0] = np.nan
    brg = col(loc, "bearing"); brg[brg < 0] = np.nan
    Re = 6378137.0; d2r = math.pi / 180
    E = (lon - lon[0]) * d2r * Re * math.cos(lat[0] * d2r)
    N = (lat - lat[0]) * d2r * Re
    seg = np.hypot(np.diff(E), np.diff(N))
    dist = np.concatenate([[0], np.cumsum(seg)])
else:
    t = E = N = spd = brg = dist = hacc = np.array([])

# ---------- 1) TRACK ----------
ax1 = fig.add_subplot(2, 2, 1)
if loc and len(E) > 1:
    sc = ax1.scatter(E, N, c=t, cmap="viridis", s=18)
    ax1.plot(E, N, "-", color="#999", lw=0.6, zorder=0)
    ax1.scatter(E[0], N[0], c="green", s=90, marker="o", label="Baslangic", zorder=3)
    ax1.scatter(E[-1], N[-1], c="red", s=90, marker="s", label="Bitis", zorder=3)
    ax1.set_aspect("equal", "datalim")
    ax1.legend(loc="best", fontsize=8)
    plt.colorbar(sc, ax=ax1, label="t (s)")
    ax1.set_title(f"GPS Track  (toplam {dist[-1]:.0f} m, hAcc med {np.median(hacc):.0f} m)")
else:
    ax1.text(0.5, 0.5, "GPS verisi yok / yetersiz", ha="center", va="center")
    ax1.set_title("GPS Track")
ax1.set_xlabel("Dogu (m)"); ax1.set_ylabel("Kuzey (m)"); ax1.grid(alpha=0.3)

# ---------- 2) HIZ ----------
ax2 = fig.add_subplot(2, 2, 2)
if loc:
    ax2.plot(t, spd, "-o", ms=3, color="#1f77b4")
    ax2.set_title(f"Hiz (GPS)  — maks {np.nanmax(spd):.1f} m/s")
ax2.set_xlabel("t (s)"); ax2.set_ylabel("m/s"); ax2.grid(alpha=0.3)

# ---------- 3) MESAFE ----------
ax3 = fig.add_subplot(2, 2, 3)
if loc:
    ax3.plot(t, dist, "-", color="#2ca02c", lw=2)
    ax3.set_title(f"Kat edilen mesafe — {dist[-1]:.0f} m")
ax3.set_xlabel("t (s)"); ax3.set_ylabel("m"); ax3.grid(alpha=0.3)

# ---------- 4) YONELIM (baslangica gore donus) ----------
ax4 = fig.add_subplot(2, 2, 4)
if gyr:
    tg = col(gyr, "seconds_elapsed"); tg = tg - tg[0]
    gz = col(gyr, "z")
    rel = np.zeros(len(tg))  # jiro entegrali + deadband (uygulamadaki mantik)
    for i in range(1, len(tg)):
        dt = tg[i] - tg[i - 1]
        rel[i] = rel[i - 1] + (-gz[i] * dt if abs(gz[i]) > 0.03 and dt < 0.5 else 0.0)
    ax4.plot(tg, np.degrees(rel), "-", color="#9467bd", lw=2, label="Donus (jiro, baslangic=0)")
    ax4.axhline(0, color="#bbb", lw=0.8)
if loc and np.isfinite(brg).any():
    ax4.plot(t, brg, ".", ms=4, color="#ff7f0e", alpha=0.6, label="GPS rotasi (mutlak)")
ax4.set_title("Yonelim  (+ saga, - sola)")
ax4.set_xlabel("t (s)"); ax4.set_ylabel("derece"); ax4.grid(alpha=0.3); ax4.legend(fontsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.97])
out = sess + "_track.png"
plt.savefig(out, dpi=130)
print("Kaydedildi:", out)
