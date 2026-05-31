# -*- coding: utf-8 -*-
# GPS track'i GERCEK HARITA uzerinde gosterir (GPSLogger gibi).
#  - <session>_map.html : OpenStreetMap uzerinde interaktif track (tarayicida ac)
#  - <session>_map.png  : harita zeminli statik goruntu (contextily; internet gerekir)
# Kullanim: python3 plot_track_map.py data_m34/session_YYYYMMDD_HHMMSS
import sys, os, csv, math
import numpy as np

sess = (sys.argv[1] if len(sys.argv) > 1 else "data_m34/session_20260530_151415").rstrip("/")

def load(name):
    p = os.path.join(sess, name)
    return list(csv.DictReader(open(p, newline=""))) if os.path.exists(p) else []
def col(rows, k): return np.array([float(r[k]) for r in rows])

loc = load("Location.csv")
if not loc:
    print("Location.csv yok/bos:", sess); sys.exit(0)

lat = col(loc, "latitude"); lon = col(loc, "longitude")
hacc = col(loc, "horizontalAccuracy")
spd = col(loc, "speed"); spd[spd < 0] = np.nan
pts = list(zip(lat, lon))
base = os.path.basename(sess)

# ---------- 1) Interaktif OSM haritasi (folium) ----------
import folium
m = folium.Map(location=[float(lat.mean()), float(lon.mean())], zoom_start=17,
               tiles="OpenStreetMap")
folium.PolyLine(pts, color="#1f77b4", weight=4, opacity=0.8,
                tooltip=f"{base}").add_to(m)
# her nokta: dogrulukla orantili daire
for la, lo, h in zip(lat, lon, hacc):
    folium.Circle([la, lo], radius=max(2, float(h)), color="#ff7f0e",
                  fill=True, fill_opacity=0.06, weight=1).add_to(m)
folium.Marker(pts[0], tooltip="Baslangic",
              icon=folium.Icon(color="green", icon="play")).add_to(m)
folium.Marker(pts[-1], tooltip="Bitis",
              icon=folium.Icon(color="red", icon="stop")).add_to(m)
m.fit_bounds([[lat.min(), lon.min()], [lat.max(), lon.max()]])
html = sess + "_map.html"
m.save(html)
print("Interaktif harita:", html)

# ---------- 2) Statik harita zeminli PNG (contextily) ----------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import contextily as ctx
    R = 6378137.0
    x = np.radians(lon) * R
    y = np.log(np.tan(np.pi / 4 + np.radians(lat) / 2)) * R
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.plot(x, y, "-", color="#1f77b4", lw=3, zorder=3)
    ax.scatter(x, y, c=np.arange(len(x)), cmap="viridis", s=20, zorder=4)
    ax.scatter(x[0], y[0], c="green", s=120, marker="o", zorder=5, label="Baslangic")
    ax.scatter(x[-1], y[-1], c="red", s=120, marker="s", zorder=5, label="Bitis")
    pad = max(60, (x.max() - x.min()), (y.max() - y.min())) * 0.3 + 30
    ax.set_xlim(x.min() - pad, x.max() + pad); ax.set_ylim(y.min() - pad, y.max() + pad)
    ax.set_aspect("equal")
    ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:3857")
    ax.set_axis_off(); ax.legend(loc="best")
    ax.set_title(f"GPS Track — {base}  (hAcc med {np.median(hacc):.0f} m)")
    png = sess + "_map.png"
    plt.savefig(png, dpi=140, bbox_inches="tight")
    print("Statik harita:", png)
except Exception as e:
    print("Statik harita atlandi (internet/tile gerekebilir):", e)
