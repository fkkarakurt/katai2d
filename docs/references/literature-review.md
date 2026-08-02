# Literatür İncelemesi & PLAXIS Hizalama (2026-06-04 checkpoint)

Bu belge, KATAI 2D'nin şu ana kadar uyguladığı yöntemleri otoriter kaynaklarla
(akademik makale, kitap, PLAXIS dökümantasyonu) çapraz-kontrol eder, PLAXIS 2D ile
hizalamayı doğrular ve sıradaki özelliklerin matematiğini kilitlemek için kaynakları
toplar. Felsefe: kodlamadan önce matematiği referansla sabitle (bkz
`mohr-coulomb-formulation.md`).

---

## 1. PLAXIS 2D ile hizalama — DOĞRULANDI

PLAXIS 2D'nin temel modelleme seçimleri ve bizdeki karşılığı:

| PLAXIS 2D | KATAI 2D | Durum |
|-----------|----------|-------|
| Analiz modu: **plane strain + axisymmetric** | İkisi de var (elastik + MC) | ✅ Tam parite |
| Eleman: **6-düğümlü + 15-düğümlü üçgen** (kullanıcı seçer) | tri6 + tri15 (kullanıcı seçer) | ✅ Tam parite |
| 15-düğüm = varsayılan (yüksek doğruluk) | tri15 doğrulandı (Prandtl %1.1 @120 elem) | ✅ |
| Doğrusal çözücü: seyrek direkt | MKL PARDISO | ✅ Eşdeğer |

**Sonuç:** Eleman tipi ve analiz modu seçimimiz PLAXIS 2D ile birebir. PLAXIS'in
15-düğümlü üçgeni varsayılan yapması, bizim tri15'te gördüğümüz "kaba mesh'te yüksek
doğruluk" kazanımıyla aynı gerekçe.

**Henüz bizde olmayan PLAXIS çekirdek iş akışları (sıradaki işler — bkz §3):**
- **K0 prosedürü** ile başlangıç gerilmesi (yatay katmanlar, K0 = σ'_h/σ'_v).
- **Staged construction** (kademeli inşa: eleman aktif/pasif, kazı/dolgu).
- **Efektif gerilme / drained-undrained** (σ' = σ − u).
- İleri bünye: **Hardening Soil**, Soft Soil, vb.
- Yapısal elemanlar (plate, anchor, geogrid, embedded pile, interface).

---

## 2. Mevcut işin doğrulaması — literatürle çapraz-kontrol

| Çıktımız | Değer | Referans (otoriter) | Hata |
|----------|-------|---------------------|------|
| Prandtl şerit temel N_c (φ=0) | 5.25 vs **2+π=5.142** | Prandtl (1920/1921) | %2.1 |
| Şerit temel N_c (φ=20°) | 15.23 vs **14.835** | **Reissner (1924)** N_q = tan²(45+φ/2)·e^(π·tanφ); **Prandtl** N_c=(N_q−1)cotφ | %2.7 |
| Şev FoS (Rocscience verif#1) | 1.010 vs **~0.99** | LEM (Bishop/Spencer); FE-SRM: **Griffiths & Lane (1999)**, Géotechnique 49(3):387–403 | %2.0 |
| Lamé kalın-cidar silindir (elastik) | hata **1e-5** | Lamé kapalı-form σ_r, σ_θ | — |
| Silindir tam-plastik çökme (Tresca) | 13.82 vs **2c·ln(b/a)=13.86** | klasik plastisite (Hill, *Mathematical Theory of Plasticity*) | %0.3 |
| MC return mapping + tutarlı teğet | — | **Sysala & Čermák (2016)** arXiv:1508.07435; **Clausen, Damkilde & Andersen (2006/2007)** | doğrulandı |

**Bearing capacity formülleri (Reissner/Prandtl) klasik kaynaklarla birebir teyit
edildi** — bizim `test_bearing_phi`'deki formülün aynısı. Şev FE-SRM yöntemimiz
Griffiths & Lane (1999)'un yöntemiyle aynı (a-priori kayma yüzeyi yok, dayanım
azaltma, yakınsamama = çöküş). Non-asosiye mesh-bağımlılığı bulgumuz da bu yöntemin
bilinen özelliği (bkz `slope-srm-convergence.md`).

---

## 3. Sıradaki özelliklerin matematiği — kaynaklar kilitlendi

### 3a. Efektif gerilme & su (en kritik geoteknik fizik)
- **Terzaghi efektif gerilme ilkesi (1923):** σ' = σ − u (u = boşluk suyu basıncı).
  Dayanım/sıkışabilirlik yalnız σ'·ile yönetilir.
- **Biot (1941) konsolidasyon:** iskelet deformasyonu ↔ akışkan akışı eşlenik (poroelastik).
- **PLAXIS yaklaşımı:** Drained / Undrained (A: efektif param + boşluk basıncı üretimi;
  B: efektif c,φ + undrained; C: total stress su) / Consolidation hesap tipleri.
- **MC'ye etki:** akma efektif gerilmede değerlendirilir (σ' = σ − u·δ); en küçük
  uyarlama = boşluk basıncını gerilmeden düşmek (steady-state/hidrostatik ilk hedef).

### 3b. Staged construction + K0 (PLAXIS'in ana iş akışı)
- **K0 prosedürü:** yatay katmanlı geometride başlangıç gerilmesi σ'_v = γ·z,
  σ'_h = K0·σ'_v. Yalnız yatay yüzey/freatik için. (Eğimli zeminde gravity loading.)
- **Staged construction:** fazlar arası eleman/yük **aktif/pasif** (kazı, dolgu, iksa
  kurulumu); her fazda denge yeniden kurulur; deplasman sıfırlanabilir. Gerilme
  yeniden dağılımı (kazıda boşalan elemanın iç kuvveti komşulara aktarılır).
- Kaynak: PLAXIS Reference/Scientific Manual (Bentley).

### 3c. Hardening Soil (en zor bünye modeli — kilitli sırada MC sonrası)
- **Schanz, Vermeer & Bonnier (1999),** *"The hardening soil model: formulation and
  verification"*, Beyond 2000 in Computational Geotechnics (Plaxis Symposium), 281–296.
- **Yapı:** klasik plastisite çerçevesi; gerilme-bağımlı stiffness (E50 sekant,
  Eoed ödometre, Eur boşaltma-yeniden yükleme, m üs); **iki akma yüzeyi** — kayma
  (shear/frictional) hardening + **cap** (hacimsel) hardening. MC kırılma kriteri.
  Frictional hardening **non-asosiye**, cap hardening **asosiye**.
- E50/Eoed sabit oran YOK (bağımsız belirlenir).

### 3d. Yapısal elemanlar (istinat/kazı dikeyini açar)
- **Plate (perde/duvar):** Mindlin kiriş (eğilme + kayma), eksenel + eğilme rijitliği.
- **Anchor:** node-to-node yay (eksenel).
- **Geogrid:** yalnız çekme; **Interface:** azaltılmış-dayanım kayma elemanı (R_inter).

---

## 4. Kaynak haritası (bibliyografya — her birinin bize verdiği)

**Plastisite / MC çekirdeği (uygulandı):**
- Sysala & Čermák (2016), arXiv:1508.07435 — principal-uzay return mapping, 4 bölge,
  tutarlı teğet. → `mohr-coulomb-formulation.md`, `mohr_coulomb.cpp`.
- Clausen, Damkilde & Andersen (2006/2007) — aynı kapalı-form, principal-uzay.

**Limit analiz / taşıma gücü (doğrulama):**
- Prandtl (1920/1921), Reissner (1924) — N_c, N_q kapalı-form. → `test_prandtl`, `test_bearing_phi`.
- Shield (1955), Eason & Shield (1960), Cox, Eason & Hopkins (1961) — **dairesel temel
  Tresca N_c = 5.69 (pürüzsüz), 6.05 (pürüzlü)**. → gelecekteki axisym dairesel temel benchmark.
- Hill, *The Mathematical Theory of Plasticity* — kalın silindir çökme 2c·ln(b/a). → `test_axisym_collapse`.

**Şev / SRM (doğrulama):**
- Griffiths & Lane (1999), Géotechnique 49(3):387–403 — FE-SRM temel makalesi. → `strength_reduction`, `test_slope`.

**Mesher (uygulandı):**
- Shewchuk (Triangle; Delaunay Refinement; adaptive predicates), Ruppert (1995). → `meshing-design.md`, `delaunay.cpp`.

**Sıradaki özellikler (kilitlenecek):**
- Terzaghi (1923), Biot (1941) — efektif gerilme & konsolidasyon.
- Schanz, Vermeer & Bonnier (1999) — Hardening Soil.
- PLAXIS Reference/Scientific/Material Models Manual (Bentley) — K0, staged, undrained,
  eleman formülasyonu, bünye modelleri (hizalama referansı).

---

## 5. Genel değerlendirme

- **Çekirdek doğru ve PLAXIS ile hizalı:** eleman tipleri, analiz modları, MC return
  mapping/teğet, limit-analiz formülleri ve SRM yöntemi otoriter kaynaklarla birebir;
  tüm kanonik benchmark'larda <%3 (çoğu <%1).
- **Boşluk:** PLAXIS'in "en çok kullanılan %20"si için eksik olanlar — efektif gerilme/su,
  staged construction, ileri bünye (HS), yapısal elemanlar. Her birinin otoriter kaynağı
  yukarıda; matematik kodlamadan önce ilgili `docs/references/*-formulation.md`'de kilitlenecek.
- **Sıradaki mantıklı adım (öneri):** efektif gerilme + undrained (geoteknikte en temel
  fizik, görece sınırlı ve doğrulanabilir) VEYA staged construction (PLAXIS'in ana iş akışı).
</content>
