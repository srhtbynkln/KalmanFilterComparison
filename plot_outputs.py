# -*- coding: utf-8 -*-
# GPS + ivme verisi ve fuzyon ciktilarinin grafigi (marmaray)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import validate as v   # dogrulanmis fusionEKF/imuLoc ve veriler

# --- filtreleri calistir ---
fn, fe, fv, bx, by = v.fusionEKF()
inx, iex, iv = v.imuLoc()

t, tg = v.t, v.tg
ax, ay, az = v.ax, v.ay, v.az
gN, gE, spd, hacc = v.gN, v.gE, v.spd, v.hacc

# konum hatasi (GPS zamanlarinda)
ffn = np.interp(tg, t, fn); ffe = np.interp(tg, t, fe)
iin = np.interp(tg, t, inx); iie = np.interp(tg, t, iex)
good = (hacc < 50)
err_f = np.sqrt((ffn-gN)**2 + (ffe-gE)**2)
err_i = np.sqrt((iin-gN)**2 + (iie-gE)**2)

fig, axs = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Marmaray (telefon) — GPS + Ivme verisi ve Fuzyon Ciktilari', fontsize=15, fontweight='bold')

# 1) Yorunge
a = axs[0,0]
a.plot(gE, gN, 'k.', ms=3, label='GPS (ham)')
a.plot(fe, fn, 'b-', lw=1.3, label='fusionEKF')
a.plot(iex, inx, 'r-', lw=0.8, alpha=0.8, label='imuLoc (saf ivme)')
a.set_xlabel('Dogu (m)'); a.set_ylabel('Kuzey (m)'); a.set_title('Yorunge')
a.axis('equal'); a.grid(True, alpha=0.3); a.legend()

# 2) Yorunge (sadece GPS vs fusion, yakin plan)
a = axs[0,1]
a.plot(gE, gN, 'k.', ms=4, label='GPS (ham)')
a.plot(fe, fn, 'b-', lw=1.3, label='fusionEKF')
a.set_xlabel('Dogu (m)'); a.set_ylabel('Kuzey (m)')
a.set_title('Yorunge — GPS vs fusionEKF (yakin)')
a.axis('equal'); a.grid(True, alpha=0.3); a.legend()

# 3) HIZ
a = axs[0,2]
a.plot(tg, spd, 'k.', ms=3, label='GPS speed (Doppler)')
a.plot(t, fv, 'b-', lw=1.0, label='fusionEKF')
a.plot(t, iv, 'r-', lw=0.6, alpha=0.7, label='imuLoc')
a.set_xlabel('Zaman (s)'); a.set_ylabel('Hiz (m/s)'); a.set_title('HIZ')
a.set_ylim(0, 30); a.grid(True, alpha=0.3); a.legend()

# 4) Ham ivme
a = axs[1,0]
a.plot(t, ax, lw=0.4, label='ax', alpha=0.8)
a.plot(t, ay, lw=0.4, label='ay', alpha=0.8)
a.plot(t, az, lw=0.4, label='az', alpha=0.6)
a.set_xlabel('Zaman (s)'); a.set_ylabel('Ivme (m/s^2)')
a.set_title('Ham ivmeolcer (lineer, yercekimsiz)'); a.set_ylim(-5, 5)
a.grid(True, alpha=0.3); a.legend()

# 5) Konum hatasi
a = axs[1,1]
a.plot(tg, err_f, 'b-', lw=1.0, label='fusionEKF')
a.plot(tg, err_i, 'r-', lw=1.0, label='imuLoc (saf ivme)')
a.set_xlabel('Zaman (s)'); a.set_ylabel('Konum hatasi (m)')
a.set_title("GPS'e gore konum hatasi (log)"); a.set_yscale('log')
a.grid(True, alpha=0.3, which='both'); a.legend()

# 6) Kestirilen ivme bias'i
a = axs[1,2]
a.plot(t, bx, 'g-', lw=1.0, label='b_ax')
a.plot(t, by, 'm-', lw=1.0, label='b_ay')
a.set_xlabel('Zaman (s)'); a.set_ylabel('Bias (m/s^2)')
a.set_title('fusionEKF — kestirilen ivme bias''i')
a.grid(True, alpha=0.3); a.legend()

plt.tight_layout(rect=[0,0,1,0.97])
out = 'Plot/cikti_grafikleri.png'
plt.savefig(out, dpi=130)
print('Kaydedildi:', out)

# ozet
rmse_f = np.sqrt(np.mean(err_f[good]**2))
rmse_i = np.sqrt(np.mean(err_i[good]**2))
print(f'fusionEKF PosRMSE = {rmse_f:.1f} m   |   imuLoc PosRMSE = {rmse_i:.0f} m')
