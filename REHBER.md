# Telefon İvmesinden Hız ve Konum Çıkarma — Rehber

Bu doküman, **basit bir telefonun GPS + ivmeölçer verisinden** doğru hız ve konum
çıkarma hedefine yönelik yapılan geliştirmeleri, fizikteki temel sınırı ve nasıl
çalıştırılacağını anlatır.

---

## 1. Önce kritik gerçek: ivmeölçer TEK BAŞINA konum veremez

İvmeölçer ivmeyi ölçer. Hız = ivmenin 1 kez, konum = 2 kez integralidir. Ölçümdeki
en küçük **sabit hata (bias)** bile:

- Hızda **bias × t** kadar (zamanla doğrusal),
- Konumda **0.5 × bias × t²** kadar (zamanla karesel)

hata biriktirir. Telefon MEMS ivmeölçerinde bias tipik olarak ~0.01–0.2 m/s²'dir.
Örnek: 0.1 m/s² bias → 60 saniyede konumda 0.5·0.1·60² = **180 m** hata.

**Bu yüzden "saf ivme" yöntemi (dead-reckoning) kaçınılmaz olarak sürüklenir.**
Bunu projedeki gerçek Marmaray verisinde ölçtük:

| Yöntem | Konum RMSE | Medyan (CEP50) | Maks. hata | Son nokta hatası | Hız RMSE |
|---|---|---|---|---|---|
| **fusionEKF** (GPS+IMU füzyonu) | **10.6 m** | **7.4 m** | **37 m** | **8.6 m** | **1.81 m/s** |
| imuLoc (saf ivme / DR) | 13.722 m | 11.678 m | 23.901 m | 23.901 m | 11.76 m/s |

> 36 dakikalık, 23.8 km'lik bir Marmaray yolculuğu. Saf-DR **13.7 km** kayarken,
> füzyon **metre** seviyesinde kalıyor. Sonuç tek cümlede: **doğru konum/hız için
> GPS ile füzyon şarttır; ivme yalnızca GPS örnekleri arasını doldurur ve GPS
> kesildiğinde (tünel) kısa süre köprü kurar.**
>
> **Güncelleme (tünel/başlatma iyileştirmesi):** Önceki maks. hata 336 m, GPS
> kaydının cold-start'ından geliyordu — GNSS logu çöp fix'lerle (hAcc 1437 m'ye
> kadar) başlıyor, filtre bu çöpte küçük P ile "emin" başlatılıp ilk iyi fix'e
> 30+ saniye yavaş yakınsıyordu. Filtre artık **ilk güvenilir (hAcc<50) fix'te**
> başlatılıyor; outage'dan dönüşte kovaryans şişiriliyor ve outage sırasında
> ZUPT kapatılıyor. Sonuç: **RMSE 25→10.6 m, maks. 336→37 m** (validate.py).

---

## 2. Doğru çözüm: `fusionEKF.m` (gevşek-bağlı GNSS/INS)

Durum vektörü: `x = [N, E, vN, vE, bax, bay]`
- `N, E` : konum (m, yerel Kuzey/Doğu)
- `vN, vE`: hız (m/s)
- `bax, bay`: **ivmeölçer bias'i** (filtre bunu canlı kestirir ve çıkarır)

Çalışma mantığı:
1. **Tahmin:** İvmeyi (bias çıkarılmış) yönelimle Kuzey/Doğu'ya çevirip hız ve
   konumu ilerletir.
2. **Düzeltme:** Her GPS örneğinde hem **konumu** hem de **Doppler hızını (speed)**
   ölçüm olarak kullanır. GPS hızı, hız durumunu doğrudan gözlemleyerek ivme
   sürüklenmesini sınırlar — projedeki en kritik eklenti budur.
3. **Yönelim (psi):** Jiroskopla ilerletilir, GPS rotası (course-over-ground) ile
   düzeltilir.
4. **Kötü GPS reddi:** `hAcc > 50 m` olan örnekler atlanır (Marmaray tünelinde GPS
   1400 m'ye fırlıyor). Bu sırada konum yalnızca IMU ile yürür.
5. **ZUPT:** Araç gerçekten durduğunda (ivme+jiro+hız ~0) hız sıfıra çekilir.

> Tasarım notu: Başta istatistiksel (chi-square) bir outlier kapısı denendi ama
> erken bir sapmada GPS'i reddedip filtreyi kalıcı ıraksamaya soktu; bu yüzden
> kaldırıldı, outlier savunması `hAcc` filtresine bırakıldı.

---

## 3. Temel sınır: telefon ORYANTASYONU

İvmeölçer telefonun **gövde eksenindedir**. İvmeyi Kuzey/Doğu'ya çevirmek için
telefonun 3B yönelimi (roll/pitch/yaw) gerekir. Marmaray verisinde **Orientation
(rotation vector) kaydı yok**; bu yüzden yalnızca yaw kullanıldı ve telefonun
yaklaşık yatay durduğu varsayıldı. Telefon cepte/elde eğikse bu varsayım bozulur.

### Daha iyi veri toplamak için (Sensor Logger, ücretsiz)
Bir sonraki kayıtta şu sensörleri de aç:
- **Orientation** veya **Rotation Vector** (en kritik — eğim telafisi için)
- Accelerometer (lineer), Gyroscope, Location (zaten var)
- Mümkünse Gravity ve Magnetometer
- Başta/sonda ~20 sn telefonu **hareketsiz** tut (bias kalibrasyonu için)

### Orientation desteği ARTIK KODDA HAZIR (otomatik)
Kod oryantasyon verisini kendiliğinden algılar:
- `data*/Orientation.csv` (veya `RotationVector.csv`) varsa ve `qw,qx,qy,qz`
  sütunlarını içeriyorsa → `loadAndSyncSensors.m` okur, IMU zamanına interpole
  eder, `ds.has_orientation = true` yapar.
- O zaman `fusionEKF.m` ve `imuLoc.m` ivmeyi **tam 3B quaternion rotasyonuyla**
  Kuzey/Doğu'ya çevirir (eğim/serbest tutuş telafi edilir; yerçekimi dünya
  çerçevesinde yalnız dikey eksende olduğundan yatayı etkilemez → gravity otomatik).
- Orientation **yoksa** kod eskisi gibi **yaw-only** çalışır.

> Konvansiyon: quaternion gövde→DÜNYA (ENU: X=Doğu, Y=Kuzey, Z=Yukarı) kabul edilir
> (Android rotation vector standardı). iOS export'unda eksen farklıysa `quatRotate`
> çıkışındaki `aE/aN` eşlemesini güncellemek yeterli.

---

## 4. Nasıl çalıştırılır

```matlab
% main.m içinde veri setini seç:
dataDir = 'data_marmaray/';   % telefon (Marmaray)
% dataDir = 'data/';          % tren (Modena RTK)

>> main
```

Çıktılar:
- Konsolda **doğruluk metrikleri tablosu** (`metrics.csv`)
- `Plot/fusion_vs_dr.png` : yörünge, hız, konum hatası, kestirilen bias
- `fusion_predictions.csv`, `imu_predictions.csv` : Lat/Lon/Hız

Python ile hızlı doğrulama (MATLAB'sız, sayıları görmek için):
```bash
python validate.py
```

---

## 5. Yapılan düzeltmeler (özet)

| Dosya | Düzeltme |
|---|---|
| `loadAndSyncSensors.m` | Telefon ivmesindeki hatalı `×9.81` kaldırıldı (veri zaten m/s²); GPS **speed/vN/vE** okunuyor; anlık `dt`; `gravity_present` bayrağı |
| `KF/SageHusaKF/STF` | Yerçekimsiz veriyi 9.81 kat şişiren `scale_acc` heuristiği kaldırıldı |
| `STF.m` | İçine yanlışlıkla yazılmış `imuLoc` kopyası yerine **doğru Strong Tracking Filter** geri getirildi |
| `EKF.m`, `EKF_strong_tracking.m` | Yerçekimi çıkarma `gravity_present`'a göre koşullu (telefonda çıkarılmaz) |
| `fusionEKF.m` | **YENİ** — GPS konum+hız füzyonu, bias kestirimi, ZUPT (asıl çözüm) |
| `fusionEKF.m` | **TÜNEL/BAŞLATMA** — ilk güvenilir fix'te başlat (cold-start çöpünü atla); outage dönüşünde kovaryans şişir; outage sırasında ZUPT kapat → maks. hata 336→37 m |
| `imuLoc.m` | Dürüst saf-DR: bias çıkarma + alçak geçiren filtre + ZUPT (karşılaştırma) |
| `computeMetrics.m` | **YENİ** — tüm filtreler için konum/hız RMSE tablosu |
| `main.m` | Tüm filtreler + metrikler + odak grafiği bağlandı |

---

## 6. Beklentiler (ne gerçekçi?)

- **Hız:** GPS Doppler + IMU ile **~1–2 m/s** doğruluk beklenir (elde edildi: 1.79 m/s).
- **Konum (GPS varken):** GPS doğruluğu mertebesinde, **~5–10 m** (CEP50 7.5 m).
- **Konum (GPS kesik, tünel):** IMU köprüsü; saniyeler içinde onlarca metre büyür.
  Başlatma/outage iyileştirmesinden sonra Marmaray'da maks. **37 m** (eski 336 m,
  GPS cold-start'tan geliyordu). Gerçek uzun tünellerde drift yine süreye bağlı
  büyür; **Orientation (eğim telafisi) + daha iyi bias modeli** bunu daha da azaltır.
- **Saf ivme (GPS'siz):** dakikalar içinde kilometrelerce kayar — **konum için
  kullanılamaz**, yalnızca çok kısa süreli (saniyeler) köprüleme için anlamlıdır.
