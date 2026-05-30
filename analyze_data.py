# -*- coding: utf-8 -*-
import csv, math, statistics as st

def load(path):
    with open(path, newline='') as f:
        r = csv.DictReader(f)
        rows = list(r)
    return rows

def col(rows, name):
    out = []
    for x in rows:
        try: out.append(float(x[name]))
        except: out.append(float('nan'))
    return out

def stats(v):
    v = [x for x in v if not math.isnan(x)]
    if not v: return "bos"
    return f"n={len(v)} mean={st.mean(v):+.4f} std={st.pstdev(v):.4f} min={min(v):+.3f} max={max(v):+.3f}"

print("="*70)
print("MARMARAY (telefon / Sensor Logger)")
print("="*70)
acc = load('data_marmaray/Accelerometer.csv')
gyr = load('data_marmaray/Gyroscope.csv')
loc = load('data_marmaray/Location.csv')

t = col(acc, 'seconds_elapsed')
dts = [t[i+1]-t[i] for i in range(len(t)-1)]
print(f"\n-- Zaman/orneklem --")
print(f"sure: {t[-1]-t[0]:.1f} s ({(t[-1]-t[0])/60:.1f} dk),  dt mean={st.mean(dts):.5f}s -> {1/st.mean(dts):.1f} Hz, dt std={st.pstdev(dts):.5f}")

ax,ay,az = col(acc,'x'),col(acc,'y'),col(acc,'z')
mag = [math.sqrt(ax[i]**2+ay[i]**2+az[i]**2) for i in range(len(ax))]
print(f"\n-- Ivme (ham, olceksiz) --")
print(f"ax: {stats(ax)}")
print(f"ay: {stats(ay)}")
print(f"az: {stats(az)}")
print(f"|a|: {stats(mag)}   <-- ~9.8 ise yercekimi VAR, ~0 ise lineer(yercekimsiz)")
print(f"ilk 200 ornek |a| ort: {st.mean(mag[:200]):.4f}  (arac duruyorsa bu ~yercekimi veya ~0)")

gx,gy,gz = col(gyr,'x'),col(gyr,'y'),col(gyr,'z')
gmag=[math.sqrt(gx[i]**2+gy[i]**2+gz[i]**2) for i in range(len(gx))]
print(f"\n-- Jiroskop (ham) --")
print(f"gz: {stats(gz)}")
print(f"|g|: {stats(gmag)}   <-- max>10 ise deg/s, kucukse rad/s")

lt = col(loc,'seconds_elapsed')
ldt=[lt[i+1]-lt[i] for i in range(len(lt)-1)]
hacc=col(loc,'horizontalAccuracy'); spd=col(loc,'speed'); brg=col(loc,'bearing')
lat=col(loc,'latitude'); lon=col(loc,'longitude')
print(f"\n-- GPS/Location --")
print(f"orneklem: {st.mean(ldt):.3f}s -> {1/st.mean(ldt):.2f} Hz")
print(f"hAcc: {stats(hacc)}")
hb=[h for h in hacc if not math.isnan(h)]
print(f"  hAcc>15 olan oran: {100*sum(1 for h in hb if h>15)/len(hb):.1f}%   hAcc>50: {100*sum(1 for h in hb if h>50)/len(hb):.1f}%")
print(f"speed: {stats(spd)}   <-- -1 = gecersiz")
sv=[s for s in spd if s>=0]
print(f"  gecerli speed orani: {100*len(sv)/len(spd):.1f}%  gecerliyse max={max(sv) if sv else 0:.2f} m/s")
print(f"bearing: {stats(brg)}  (-1=gecersiz)")
# kat edilen mesafe
R=6378137; d2r=math.pi/180
dist=0
for i in range(1,len(lat)):
    if any(math.isnan(x) for x in (lat[i],lat[i-1],lon[i],lon[i-1])): continue
    dn=(lat[i]-lat[i-1])*d2r*R
    de=(lon[i]-lon[i-1])*d2r*R*math.cos(lat[i]*d2r)
    dist+=math.sqrt(dn*dn+de*de)
print(f"\nGPS ile kat edilen toplam mesafe: {dist:.0f} m ({dist/1000:.2f} km)")
print(f"ort hiz (mesafe/sure): {dist/(lt[-1]-lt[0]):.2f} m/s = {3.6*dist/(lt[-1]-lt[0]):.1f} km/h")

# IMU-GPS zaman ortusmesi
print(f"\n-- Zaman hizalama --")
print(f"IMU  t: {t[0]:.2f} .. {t[-1]:.2f}")
print(f"GPS  t: {lt[0]:.2f} .. {lt[-1]:.2f}")

print("\n"+"="*70)
print("MODENA (tren / RTK) - karsilastirma")
print("="*70)
try:
    imu=load('data/clean_imu.csv'); gn=load('data/clean_gnss.csv')
    ti=col(imu,'t'); dti=[ti[i+1]-ti[i] for i in range(min(50000,len(ti)-1))]
    print(f"IMU: {len(imu)} satir, ~{1/st.mean(dti):.0f} Hz")
    axm,aym,azm=col(imu,'acc_x'),col(imu,'acc_y'),col(imu,'acc_z')
    print(f"acc_z: {stats(azm[:50000])}  <-- ~9.8 ise yercekimi var")
    h2=col(gn,'hAcc'); s2=col(gn,'speed')
    print(f"GNSS hAcc: {stats(h2)}")
    print(f"GNSS speed: {stats(s2)}")
except Exception as e:
    print("modena okunamadi:", e)
