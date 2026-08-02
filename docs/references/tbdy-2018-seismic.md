# TBDY 2018 sismik tasarım spektrumu + response spectrum — formülasyon (Faz D4a)

KATAI 2D sismik izinin **TBDY 2018 (Türkiye Bina Deprem Yönetmeliği)** entegrasyonu: yatay elastik
tasarım spektrumu + yerel zemin etki katsayıları + ivme-kaydından response spectrum. Bu, v1.0 bayrak
farklılaştırıcısının (TBDY-native sismik) çekirdeğidir. Matematik **resmî AFAD TBDY 2018 metninden**
birebir kilitlendi. İlgili: [[project-plaxis-parity-roadmap]], [[dynamic-seismic-formulation]].

**Kaynak (kilitli, YETKİLİ):** *Türkiye Bina Deprem Yönetmeliği 2018* (AFAD, Bakanlar Kurulu 2018/11275),
Bölüm 2 "Deprem Etkisi Altında Binaların Tasarımı İçin Deprem Yer Hareketi" — §2.3.2 (yerel zemin etki
katsayıları, Tablo 2.1/2.2), §2.3.4.1 (yatay elastik tasarım spektrumu, Denk. 2.2–2.3). Resmî PDF:
https://www.afad.gov.tr/kurumlar/afad.gov.tr/2309/files/TBDY_2018.pdf. Değerler PDF'ten `pdftotext` ile
çıkarılıp teyit edildi (ikincil web kaynaklarında hatalar vardı — ör. tablolar 6 sütun, 5 değil).

## 1. Tasarım spektral ivme katsayıları (§2.3.2.2)
Harita spektral ivme katsayıları S_S (kısa periyot) ve S_1 (1.0 s), AFAD Türkiye Deprem Tehlike
Haritası'ndan (https://tdth.afad.gov.tr) DD-1…DD-4 düzeyi + koordinat için alınır. Yerel zemin
etki katsayılarıyla **tasarım** spektral ivme katsayıları:
```
S_DS = S_S · F_S          S_D1 = S_1 · F_1
```

## 2. Yerel zemin etki katsayıları (§2.3.2.1, Tablo 2.1/2.2) — RESMÎ DEĞERLER
Zemin sınıfları ZA–ZF (Vs30 / (N60)30 / (cu)30, Tablo 16.1). ZF = sahaya-özel analiz (tablosuz).
Ara S_S/S_1 için **doğrusal enterpolasyon**; aralık dışında uçtaki değere sabitlenir (clamp).

**Tablo 2.1 — F_S (kısa periyot):** sütunlar S_S = 0.25, 0.50, 0.75, 1.00, 1.25, 1.50
```
ZA:  0.8  0.8  0.8  0.8  0.8  0.8
ZB:  0.9  0.9  0.9  0.9  0.9  0.9
ZC:  1.3  1.3  1.2  1.2  1.2  1.2
ZD:  1.6  1.4  1.2  1.1  1.0  1.0
ZE:  2.4  1.7  1.3  1.1  0.9  0.8
```
**Tablo 2.2 — F_1 (1.0 s periyot):** sütunlar S_1 = 0.10, 0.20, 0.30, 0.40, 0.50, 0.60
```
ZA:  0.8  0.8  0.8  0.8  0.8  0.8
ZB:  0.8  0.8  0.8  0.8  0.8  0.8
ZC:  1.5  1.5  1.5  1.5  1.5  1.4
ZD:  2.4  2.2  2.0  1.9  1.8  1.7
ZE:  4.2  3.3  2.8  2.4  2.2  2.0
```

## 3. Yatay elastik tasarım spektrumu (§2.3.4.1, Denk. 2.2–2.3)
Yatay elastik tasarım spektral ivmesi S_ae(T) [g cinsinden], doğal titreşim periyodu T'ye bağlı:
```
S_ae(T) = (0.4 + 0.6·T/T_A)·S_DS        (0 ≤ T ≤ T_A)
S_ae(T) = S_DS                           (T_A ≤ T ≤ T_B)
S_ae(T) = S_D1 / T                       (T_B ≤ T ≤ T_L)
S_ae(T) = S_D1·T_L / T²                  (T_L < T)
```
Köşe periyotları (Denk. 2.3):  **T_A = 0.2·S_D1/S_DS**,  **T_B = S_D1/S_DS**.  Sabit-yerdeğiştirme
geçiş periyodu **T_L = 6 s** (§2.3.4.1). Yerdeğiştirme spektrumu (Denk. 2.4): S_de(T) = (T²/4π²)·g·S_ae(T).
(Düşey spektrum §2.3.5: T_AD=T_A/3, T_BD=T_B/3, TLD=T_L, 0.32+0.48·T/T_AD kolu — sonraki faz.)

## 4. Response spectrum (ivme-kaydından) — dinamik motora bağ
Bir ivme kaydı a_g(t)'nin **%5 sönümlü elastik response spectrum**'u: her periyot T için tek-serbestlik
(SDOF) salınıcı (m=1, ω=2π/T, c=2ξω) −m·a_g(t) ile sürülür (Newmark; [[dynamic-seismic-formulation]]
solve_newmark yeniden kullanılır); **pseudo-spektral ivme S_a(T) = ω²·max|u_rel(t)|**. PLAXIS/sismik
standart çıktı; TBDY tasarım spektrumuyla kıyas + spektral-eşleştirme (D4) için.
Kapalı-form (kararlı-hal, sürekli sinüs a_g=A sin ω_0 t): S_a(T) = A·Rd(r,ξ), r=ω_0/ω=ω_0·T/2π,
Rd=1/√((1−r²)²+(2ξr)²). Limitler: T→0 (rijit) → S_a=PGA=A; T=T_0=2π/ω_0 (rezonans) → S_a=A/(2ξ).

## 5. KATAI çekirdek + doğrulama (D4a — test_tbdy_seismic)
- `analysis/response_spectrum.hpp`: `tbdy_elastic_spectrum(S_DS,S_D1,T[,T_L])` (§3); `tbdy_site_coefficients
  (SiteClass,S_S,S_1)→{F_S,F_1}` (§2, 6-sütun tablo + doğrusal enterpolasyon/clamp); `response_spectrum
  (accel,dt,periods,ξ)→S_a(T)` (§4, SDOF solve_newmark).
- **Doğrulama `test_tbdy_seismic`:** (a) **spektrum** köşeleri T=0→0.4 S_DS, plato S_DS, T_A/T_B geçiş
  sürekliliği, S_D1/T ve T_L sonrası S_D1 T_L/T² — worked example ile birebir; (b) **zemin katsayıları**
  tablo çapa değerleri + enterpolasyon + clamp; (c) **response spectrum** sürekli-sinüs girdisiyle S_a(T)=
  A·Rd(r,ξ): PGA limiti (T→0→A), rezonans A/(2ξ), ara-r. Tümü <%2.
- **Sonraki (D4b):** GUI — Dynamic faz tipi (PhaseType), TBDY parametreleri (S_S/S_1/zemin sınıfı veya
  doğrudan S_DS/S_D1), ivme-kaydı girişi, response-spectrum + tasarım-spektrumu çıktı grafiği; 2D dinamik
  çözüm sürücüsü (base input + absorbing/free-field sınırlar D3/D3b). NOT: ivme zaman-serisi + ara-adım
  saklama → ertelenen P3 zaman-playback ortak altyapı.


## EC8 EKİ — EN 1998-1:2004 yatay elastik tepki spektrumu (KATAI ec8_elastic_spectrum)

**Formüller (§3.2.2.2, Denk. 3.2-3.5):** dört dal — [0,T_B]: a_g·S·[1+(T/T_B)(2.5η−1)];
[T_B,T_C]: 2.5η·a_g·S (plato); [T_C,T_D]: 2.5η·a_g·S·T_C/T; [T_D,∞): 2.5η·a_g·S·T_C·T_D/T².
Sönüm düzeltmesi (Denk. 3.6): η = √(10/(5+ξ%)) ≥ 0.55 (ξ=%5'te η=1). a_g = γ_I·a_gR.

**Parametreler (EN ÖNERİLEN değerler; ulusal ekler değiştirebilir — GUI ve rapor bunu beyan eder):**
Tip 1 (M_s>5.5, Tablo 3.2) {S,T_B,T_C,T_D}: A{1.0,0.15,0.4,2.0} B{1.2,0.15,0.5,2.0}
C{1.15,0.2,0.6,2.0} D{1.35,0.2,0.8,2.0} E{1.4,0.15,0.5,2.0}. Tip 2 (Tablo 3.3):
A{1.0,0.05,0.25,1.2} B{1.35,0.05,0.25,1.2} C{1.5,0.1,0.25,1.2} D{1.8,0.1,0.3,1.2}
E{1.6,0.05,0.25,1.2}.

**Kaynak doğrulaması (2026-07-19):** tablo değerleri bağımsız mühendislik dokümantasyonuyla
(SCIA Engineer EN1998-1 yardım sayfası, help.scia.net) satır satır çapraz-doğrulandı; eğitim
bilgisiyle %100 örtüştü. V&V: test_tbdy_seismic (f) — tablolar testte BAĞIMSIZ yeniden yazılır
(çift transkripsiyon), köşe süreklilikleri, dört dalın kapalı-formları, η değerleri, γ_I lineerliği;
test_dynamic_gui EC8 bindirme kablolaması 41 periyotta bit-birebir.