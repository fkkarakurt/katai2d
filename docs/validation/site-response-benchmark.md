# Benchmark Dalga-1 — Çok-Tabaka 1D Saha Tepkisi (SHAKE-sınıfı katmanlı çözüm kıyası)

**Track 7 / ROADMAP §4 v0.6 kapısı, ilk bağımsız DİNAMİK benchmark.** Test: `test_site_response_benchmark`
(ctest, her build'de). Tarih: 2026-07-19. Risk sınıfı: **Vital** (sahaya özel deprem büyütmesi —
tasarım spektrumunun ve tüm sismik talebin girdisi).

## Ne kanıtlanıyor

KATAI'nin 2D dinamik GUI yolu (mesh → `PhaseType::Dynamic` → Newmark+Rayleigh), **katmanlı** bir
zemin profilinin frekans-bağımlı büyütmesini (transfer fonksiyonu |T(ω)| = |a_yüzey/a_taban|),
SHAKE ailesinin lineer çekirdeği olan **katmanlı viskoelastik transfer-matris çözümüne** karşı
üretir. Bu, tek-tabaka kapalı-form rezonans testlerinin (test_dynamics/test_dynamic_gui) ötesine
geçen ilk gerçek çok-tabaka kıyasıdır: empedans kontrastı rezonansları hiçbir tek-tabaka formülün
öngöremeyeceği yerlere taşır (aşağıdaki f₁ tartışması).

## Bağımsız oracle (çözücüyle SIFIR kod paylaşımı)

Rijit taban + within harmonik giriş a_g = A·e^{iωt}; tabaka başına, göreli çerçevede, FE'nin
Rayleigh sönümünün (C = αM + βK) **birebir sürekli-ortam karşılığıyla**:

    −ω²ρ* û − G* û″ = −ρA,      G* = G(1 + iωβ),      ρ* = ρ(1 − iα/ω)

Genel çözüm (tabaka j, yerel z tabaka tabanından): û_j = a_j cos(k_j z) + b_j sin(k_j z) + P_j,
k_j = ω√(ρ*_j/G*_j), P_j = Aρ_j/(ω²ρ*_j). Koşullar: tabanda û=0; ara yüzeylerde û ve τ=G*û′
sürekli; yüzeyde τ=0. Küçük kompleks lineer sistem (2 bilinmeyen/tabaka, Eigen LU). Transfer:
T(ω) = 1 − ω²û(yüzey)/A.

**Kritik metodolojik nokta — sönüm birebir eşlendi:** SHAKE frekans-BAĞIMSIZ histeretik sönüm
(G(1+2iξ)) kullanır; KATAI Newmark+Rayleigh kullanır (frekans-bağımlı ξ(ω)=α/2ω+βω/2). İkisini
karşılaştırırken "yaklaşık eşdeğer ξ" bandına saklanmak yerine, oracle'a Rayleigh'in TAM kompleks
karşılığı (yukarıdaki G*, ρ*) kondu → kıyas like-for-like ve tek meşru fark uzay/zaman
ayrıklaştırmasıdır. (Histeretik-vs-Rayleigh farkı bir MODELLEME seçimidir; SHAKE'in kendi sönüm
modeliyle kıyas, eşdeğer-lineer iterasyonla birlikte gelecekteki genişlemedir — aşağıda.)

**Oracle'ın kendi pinleri (testin içinde):**
1. Tek tabaka, sönümsüz → klasik T = 1/cos(ωH/Vs)'e **6.0e-16** ile iner (işaret/BC kanıtı).
2. İki-tabaka temel frekansı, elle türetilebilir karakteristik denklemi sağlar:
   **tan(k₁h₁)·tan(k₂h₂)·α_imp = 1.0005** (α_imp = ρ₁Vs₁/ρ₂Vs₂; grid çözünürlüğü içinde 1).
3. Rayleigh α,β testte bağımsız kapalı-formdan yeniden hesaplanır (α=2ξω₁ω₂/(ω₁+ω₂),
   β=2ξ/(ω₁+ω₂)); çekirdek farklı eşleseydi tüm sönümlü kıyaslar patlardı.

## Profiller ve sonuçlar

**(b) Tek tabaka, sönümlü** — H=20 m, Vs=126.5 m/s, ρ=2.0, ξ=%5 @ (f₁, 3f₁):

| f | |T|_FE | |T|_oracle | hata |
|---|---|---|---|
| 0.5·f₁ = 0.791 Hz | 1.4134 | 1.4126 | **+%0.05** |
| f₁ = 1.581 Hz (rezonans) | 12.743 | 12.766 | **−%0.18** |
| 2·f₁ = 3.162 Hz | 1.0076 | 0.9907 | **+%1.7** |

**(c) İki tabaka (yumuşak-üstte, gerçek empedans kontrastı α_imp=0.343)** —
üst: 8 m, Vs=120, ρ=1.8; alt: 12 m, Vs=300, ρ=2.1; ξ=%5:

| Nokta | f [Hz] | |T|_FE | |T|_oracle | hata | bant |
|---|---|---|---|---|---|
| Temelin altı (flank) | 1.502 | 1.4764 | 1.4755 | **+%0.06** | %4 |
| 1. rezonans (tepe) | 3.005 | 16.967 | 16.982 | **−%0.09** | %5 |
| Modlar arası çukur | 4.850 | 2.2440 | 2.2296 | **+%0.65** | %6 |
| 2. rezonans (tepe) | 6.695 | 7.7676 | 7.7381 | **+%0.38** | %7 |

**Fizik dersi (testin öncülünü düzeltti):** çeyrek-dalga seyahat-süresi kuralı
f₁ ≈ 1/(4Σh/Vs) = 2.34 Hz verir; GERÇEK temel 3.005 Hz'tir — kuvvetli kontrastta yumuşak tabaka
kendi rijit-taban çeyrek-dalgasına (3.75 Hz) yaslanır. Bu yüzden test kaba kuralı değil, kesin
karakteristik denklemi pinler. (Pratik çıkarım: kuvvetli kontrastlı sahada seyahat-süresi kuralına
göre seçilmiş Rayleigh hedef frekansları temelden %20+ şaşabilir — mühendise not.)

## Ürünleşen ders — otomatik f₁ + Rayleigh bandı korkuluğu

Yukarıdaki %22'lik seyahat-süresi sapması bir ürün korumasına çevrildi: **her Dynamic faz artık
modelin kendi temel frekansını hesaplar** ((K,M) üzerinde ters-kuvvet iterasyonu; 2D geometri +
tabakalar + yapılar dahil; `SolveResult.dyn_model_f1`) ve faz mesajında f₁ ile oradaki etkin
Rayleigh sönümünü (ξ_eff = α/4πf₁ + βπf₁) raporlar. Bant f₁'i kapsamıyorsa veya ξ_eff hedeften
>%20 sapıyorsa açık uyarı verir; düzgün köşeli bant (f_R2 ≈ 3f_R1) bant-içi sapmayı ≤~%13'te
tuttuğu için doğru kullanımda uyarı HİÇ tetiklenmez (ölçü, panik değil).

**V&V (iki bağımsız yol aynı sayıya):** üniform kolonda otomatik f₁ = Vs/4H **+%0.00**;
iki-tabakada otomatik f₁ = 3.0047 Hz vs katmanlı-oracle kutbu 3.0050 Hz → **−%0.01** (montajlı 2D
K,M öz-değeri ile transfer-matris rezonansı hiçbir kod paylaşmaz). Ayrımcılık: köşeli bant →
uyarı YOK; [0.5, 1.0] Hz'e sıkışmış bant (gerçek f₁=3.0 dışarıda) → uyarı VAR.

**Pratik sonucun büyüklüğü (bu profilde, dürüst):** seyahat-süresi hedefli bant bile f₂=3f_R1
sayesinde f₁'i KAPSAR; etkin sönüm %5 yerine %4.53 olur (≈%9 az sönüm ≈ rezonans genliğinde ~%10)
— fark benchmark'ta oracle'a birebir yansıtıldığı için kıyası etkilemez, ama mühendisin niyetinden
sapmadır; korkuluk tam bu sapmayı görünür kılar.

## Ölçüm yöntemi ve dürüst sınırlar

- FE ölçümü: 30 çevrim, 60 adım/çevrim; genlik kaydın SON çeyreğinden (ξ=%5'te geçici rejim
  ~e^{−0.05·2π·22} ≈ 1e-3'e düşer). |T|_FE = max|a_yüzey,toplam|/A.
- 2D kolon SH koşullarıyla (yanlar u_y=0) 1D kesme kolonuna indirgenir; mesh 0.5 m, λ_min≈17 m
  → dalga çözünürlüğü bol.
- **Kapsam dışı (bu benchmark'ın iddia ETMEDİĞİ):** eşdeğer-lineer iterasyon (G/G₀−γ, ξ−γ
  eğrileri; SHAKE'in nonlineer tarafı), gerçek deprem kaydı girişi, outcrop/deconvolution
  (absorbing taban v0.6 kaleminde), frekans-bağımsız histeretik sönümle kıyas. Bunlar Dalga-1'in
  sonraki kalemleri; bu kayıt yalnız **lineer katmanlı büyütme** iddiasını kapatır.

## Reprodüksiyon
`cmake --build build/msvc-rwdi --target test_site_response_benchmark` → çalıştır. Formülasyon
temeli: `docs/references/dynamic-seismic-formulation.md`; ilgili: Kramer (1996) §7 (katmanlı
transfer fonksiyonu), SHAKE/DEEPSOIL lineer çekirdek yöntemi.
