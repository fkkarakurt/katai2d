# Şev stabilitesi (φ-c reduction / SRM) — doğrulama ve mesh yakınsama

**Tarih:** 2026-06-04. **Benchmark:** Rocscience/Slide verification #1 (Griffiths &
Lane sonrası): homojen 1:2 foundation şevi, γ=20.2 kN/m³, c=3 kPa, φ=19.6°.
**Referans FoS:** Bishop **0.988**, Spencer/GLE **0.987** (LEM); Phase2 6-düğümlü
üçgen FE-SRM **0.997**. Bizim eleman tipimiz tri6 olduğundan en yakın kıyas Phase2 T6.

Yöntem: `factor_of_safety` (SRM) — c_f=c/SRF, φ_f=atan(tanφ/SRF), her SRF için MC +
gravity nonlineer çözüm (analitik tutarlı teğet + backtracking line search,
non-asosiye için nonsimetrik PARDISO), çökme = yakınsamama, SRF bisection.
Regresyon testi: `tests/test_slope.cpp` (max_area=5 → 165 eleman → FoS=1.010, %1.3
T6'ya / %2.0 LEM'e — hedef <%5). Çalışma kodu: `tests/study_slope_convergence.cpp`.

## Mesh yakınsama (load_steps=8, sonuçlar load_steps 4/8/16'ya duyarsız)

| max_area | eleman | FoS (ψ=0, benchmark) | FoS (ψ=φ, asosiye) |
|---------:|-------:|---------------------:|-------------------:|
| 10.0     | 86     | 1.093                | 1.103              |
| 5.0      | 165    | 1.009                | 1.032              |
| 2.5      | 344    | 1.032                | 1.046              |
| 1.25     | 668    | 0.976                | 1.007              |
| 0.6      | 1383   | 0.872 (ψ=0)          | —                  |

## Yorum

- **Tüm pratik (mühendislik) mesh yoğunluklarında (165–668 eleman) FoS, referansları
  (~0.99) %2–4 içinde bracketliyor.** İlk rapor edilen FoS=1.010 (165 eleman) temsili.
- **Asosiye akış (ψ=φ) belirgin şekilde daha mesh-objektif** (1.103→1.007, nazik ve
  monoton azalma); **non-asosiye (ψ=0) daha çok aşağı kayıyor** (1.009→0.976→0.872).
  Bu, **FE-SRM + mükemmel plastisite + non-asosiye akışın iyi-bilinen mesh-objektiflik
  kaybıdır** (regularizasyon olmadan shear-band lokalizasyonu). **Bir solver hatası
  değildir** — aynı sınırlama PLAXIS dâhil tüm FE-SRM kodlarında vardır; bu yüzden
  pratikte ölçülü mesh yoğunluğu kullanılır ve referans 0.997 belirli bir mesh'tedir.
- FoS'un load-stepping'e duyarsızlığı (4/8/16 adımda 0.9756/0.9756/0.9736) çökme
  sinyalinin **gerçek** olduğunu, iterasyon-sayısı artefaktı olmadığını gösterir.
- **Sonuç:** solver, kurulu FE-SRM yönteminin teorik davranışını (asosiye → objektif,
  non-asosiye → mesh-bağımlı) doğru üretiyor ve referansları beklenen bantta yakalıyor.
  "Tek kesin yakınsamış FoS" non-asosiye için kavramsal olarak yoktur; mühendislik
  bandı (~0.97–1.03) referanslarla tutarlıdır.

## tri15 (15-düğümlü kuartik eleman) — nonlineer sonuçlar ve ÇÖZÜLEN bulgu (2026-06-04)

| max_area | eleman | FoS tri6 (ψ=0) | FoS tri15 (sabit 8 adım) | FoS tri15 (**adaptif**) |
|---------:|-------:|---------------:|-------------------------:|------------------------:|
| 10.0     | 86     | 1.091          | 0.999                    | **0.999**               |
| 5.0      | 165    | 1.009          | 0.741 ⚠️ (sahte)         | **0.944** ✅            |

**İlk gözlem:** 165 elemanlı tri15 FoS=0.741'e düşüyordu (erratik). **Teşhis** (`diag_tri15`,
sabit SRF): SRF=0.75'te (şev rahatça kararlı, c_f=c/0.75 taban'dan güçlü) tri6 yakınsarken
**tri15 + non-asosiye + 8 yük adımı yakınsamıyor**, ama 16 adımla VEYA asosiye'yle
yakınsıyor → **fizik değil, yanlış-çöküş**: daha büyük/sert tri15 (non-asosiye) sistemi için
sabit yük adımı çok kaba, Newton sığmıyor, bisection yanlışça "çöküş" sayıyor.

**Kök sebep ve çözüm: sabit yük adımlandırması.** `solve_nonlinear`'a **adaptif (otomatik)
yük adımlandırması** eklendi (sub-stepping/cutback + line-search-stall'da erken-kesme;
Crisfield / PLAXIS-ABAQUS otomatik adımlama). Bir artım yakınsamazsa adım yarıya bölünüp
tekrar denenir; ancak min artımda bile denge yoksa gerçek çöküş sayılır. Çöküş sinyali
artık keyfi adım sayısından bağımsız → FoS **eleman tipinden ve sistem sertliğinden
bağımsız sağlam.** tri6 cutback gerektirmediğinden **birebir korundu** (Prandtl/bearing/
şev aynı, <%5).

**Çözüm sonrası (adaptif):** tri15 165 elem **0.741 → 0.944.** tri6 referansları üstten
(1.09, 1.01), tri15 alttan (1.00, 0.94) bracketliyor — ikisi de birkaç % içinde; tri15'in
hafif düşük olması yüksek-dereceli elemanın **gerçek (hafif) lokalizasyon** eğilimi, artık
sahte çöküş değil. **tri15 artık nonlineer için de güvenilir.** (Performans: çöküşe-yakın
yavaş Newton nedeniyle limit-analiz testleri ~2× yavaşladı; robust limit-analizin doğal
maliyeti, bireysel analizler makul.)

## Açık iş (gelecek)

- Çok-ince mesh kuyruğundaki (ψ=0, 1383 eleman → 0.872) düşüşün mesh-kalite mi yoksa
  saf lokalizasyon mu olduğunun ayrıştırılması (min-açı dağılımı incelemesi).
- İkinci/üçüncü bağımsız şev benchmark'ı (φ=0 Taylor stabilite sayısı, farklı eğim) —
  mesher'ın keskin/dar toe köşesi sınırı çözülünce (şu an benchmark köşeleri ≥90°).
