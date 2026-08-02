# Yapısal Elemanlar — Formülasyon (P2.4): Plate (Mindlin kiriş), Anchor, Geogrid

Desteklenmiş kazı dikeyinin (perde duvar + ankraj/strut + interface) omurgası. PLAXIS 2D'nin
en yaygın gerçek-mühendislik kullanımı. Felsefe: matematiği otoriter kaynakla kilitle, en küçük
parçayı analitik doğrula (bkz `mohr-coulomb-formulation.md`, `hardening-soil-formulation.md`).

**Kaynaklar (kilitli):**
- **PLAXIS 2D Material Models Manual, Bölüm 18 "Structural Behaviour"** (Eq 18-1…18-9) — BİRİNCİL.
  Plate=Mindlin kiriş (kayma-deformasyonlu): N=EA·ε, M=EI·κ, Q=(kEA/(2(1+ν)))·γ\*, k=5/6.
- **Zienkiewicz & Taylor (FEM), Timoshenko kiriş elemanı** — kayma kilitlenmesi (shear locking)
  ve azaltılmış/seçmeli integrasyon (reduced/selective integration).
- Goodman/Day&Potts (interface, sonraki parça), Sluis (embedded pile, sonraki parça) —
  `analysis-and-structural.md` §5-6.

## 0. Konvansiyon
2D düzlem-şekildeğiştirme. Plate = mesh kenarına oturan hat elemanı (tri6 kenarı = 3 düğüm:
2 köşe + 1 orta). Çözücü tension-pozitif (zemin ile aynı). Plate kuvvetleri **birim genişlik
başına** (per unit width, düzlem-dışı): N [kN/m], Q [kN/m], M [kNm/m]. Eksen: s = plate boyunca
(yerel x'), n = plate normali (yerel y').

## 1. Plate = Mindlin (Timoshenko) kiriş — kinematik
Her düğümde **3 DOF: (u_x, u_y, φ)** — iki öteleme (zeminle PAYLAŞILAN) + bir dönme φ (plate'e
ÖZGÜ ek DOF). Mindlin kinematiği (dönme, enine öteleme eğiminden BAĞIMSIZ → kayma deformasyonu):
```
ε  = du_s/ds                 (eksenel şekildeğiştirme; u_s = yerel eksenel öteleme)
κ  = dφ/ds                   (eğrilik)
γ  = du_n/ds − φ             (kayma şekildeğiştirmesi; Mindlin: φ ≠ du_n/ds)
```
u_s, u_n = global (u_x,u_y)'nin plate eksenine izdüşümü: u_s=u_x·c+u_y·s, u_n=−u_x·s+u_y·c
(c=cosθ, s=sinθ, θ=plate açısı).

## 2. Bünye — yapısal kuvvet-şekildeğiştirme (PLAXIS MMM Eq 18-6…18-9)
```
N = EA·ε                              (eksenel kuvvet, Eq 18-6)
M = EI·κ                              (eğilme momenti, Eq 18-9)
Q = (k·EA/(2(1+ν)))·γ = kGA·γ         (kesme kuvveti, Eq 18-8; k=5/6 kayma düzeltme faktörü)
```
- **k = 5/6** (dikdörtgen kesit kayma düzeltme faktörü).
- G = E/(2(1+ν)); plate'te EA verilir, kesme rijitliği kGA' = k·EA/(2(1+ν)) (Eq 18-8 biçimi —
  EA üzerinden, ν ile). İzotropik plate'te EA=E·d, EI=E·d³/12, **eşdeğer kalınlık d_eq=√(12·EI/EA)**.
- H (hoop, düzlem-dışı) = EA₂·ε_H; düzlem-şekildeğiştirmede ε_H=0 ⇒ H=0 (Eq 18-7). İlk uygulamada
  izotropik (EA₂=EA) ve düzlem-şekildeğiştirme ⇒ H göz ardı.
- ν: plate Poisson; izotropik. (Anizotropik plate'te ν=0 varsayılır, Eq 18-5.)

## 3. Yerel rijitlik matrisi (3-düğümlü kuadratik Timoshenko)
Yerel DOF vektörü (düğüm i için [u_s, u_n, φ]): 3 düğüm × 3 = 9 DOF. Şekil fonksiyonları N_i(ξ)
1D kuadratik (ξ∈[−1,1]; köşe ξ=∓1, orta ξ=0). Yerel rijitlik:
```
K_local = ∫_L ( EA·Bε^T Bε + EI·Bκ^T Bκ + kGA·Bγ^T Bγ ) ds
   Bε  = [dN/ds  uygulanan u_s'e]        (eksenel)
   Bκ  = [dN/ds  uygulanan φ'ye]         (eğilme)
   Bγ  = [dN/ds  u_n'e ,  −N  φ'ye]      (kesme: du_n/ds − φ)
```
- **KAYMA KİLİTLENMESİ (shear locking):** ince plate'te (d≪L) tam integrasyon kayma terimini
  aşırı-rijit yapar (γ→0 sahte kısıtı). **Çözüm: kesme terimine AZALTILMIŞ integrasyon** (bending/
  axial tam; shear bir mertebe düşük Gauss). 3-düğümlü kuadratik için: bending+axial 3-nokta,
  shear 2-nokta (seçmeli azaltılmış integrasyon — Zienkiewicz&Taylor). Bu, ince ve kalın plate'i
  birden doğru yapar.

## 4. Yerel → global dönüşüm
Yerel [u_s,u_n,φ] ile global [u_x,u_y,φ] arasında düğüm-bloğu dönüşümü T_i:
```
[u_s, u_n, φ]_i = R · [u_x, u_y, φ]_i ,   R = [[c, s, 0],[−s, c, 0],[0,0,1]]
```
K_global = T^T K_local T (T = blok-diyagonal R). φ DOF dönmez (skalar). Düz plate için θ sabit;
eğri plate'te eleman-bazlı θ.

## 5. Zemin ile birleştirme (DOF mimarisi)
- Öteleme DOF'ları (u_x,u_y) plate düğümünde **zeminle PAYLAŞILIR** (aynı global denklem indeksi)
  → plate rijitliği zemin rijitliğine eklenir (montaj toplama). Plate, zeminle uyumlu (kenar düğümleri).
- Dönme DOF'u φ plate'e ÖZGÜ → DofMap'e **EK DOF** olarak eklenir (zemin 2-DOF bloğundan sonra,
  global denklem sistemine; bkz DofMap genişletmesi). Sadece plate düğümlerinde.
- Plate ucu serbest/ankastre: φ DOF sabitlenebilir (moment-bağlantı) ya da serbest (mafsal).

## 6. Yapısal kuvvet geri-kazanımı
Gauss noktasında ε,κ,γ → N=EA·ε, M=EI·κ, Q=kGA·γ. Düğüme ekstrapolasyon (sonuç/tasarım çıktısı).
Plastisite (sonraki): M_p (max moment), N_p (max eksenel), elmas etkileşimi (Eq 18-1: |N/N_p|+|M/M_p|<1).

### 6a. Output: iç-kuvvet DİYAGRAMLARI + tablo (PLAXIS Output muadili; UYGULANDI)
`analysis/structural_forces.hpp` — yakınsamış global çözümden (`NewtonResult.displacement`, tam global-DOF
uzayı) plate/anchor/geogrid iç kuvvetlerinin DAĞILIMINI + konumunu üretir. Saf SON-İŞLEM, MKL'e bağlı değil.
PLAXIS Reference Manual "Output → Structures → Bending moments M / Shear forces Q / Axial forces N"
(birim genişlik başına kN·m/m, kN/m). `wall_force_envelope`'un yalnız max |M|,|Q|,|N| döndüren hâlinin
GENEL (DofMap-tabanlı, hem gömülü perde hem bonded plate-in-soil) + DAĞILIMLI karşılığı.

- **`plate_force_diagram`** (3- ve 5-düğüm): plate zincirini her elemanda eşit-aralıklı ξ∈[−1,+1]'de
  örnekler; her istasyonda `ForceStation{s, x, y, N, Q, M}` — yay-uzunluğu s fiziksel konumlardan birikir
  (eğik/kavisli perde de doğru). Öteleme DOF'u trans_dof≥0 ise o (bağımsız perde), aksi hâlde
  global_dof(node,c) (mesh paylaşımlı). N=EA·ε ve M=EI·κ TAM integrasyondan (lineer alanlar → her ξ'de
  süperyakınsar). **Q = BARLOW geri-kazanımı:** seçmeli azaltılmış integrasyonda kesme şekildeğiştirmesi
  yalnız AZALTILMIŞ Gauss (Barlow/optimal stress) noktalarında süperyakınsar; düğüm/uçlarda sahte salınır
  (Zienkiewicz & Taylor; Prathap). Q'yu o noktalarda (3-düğüm: ξ=±1/√3; 5-düğüm: 4-nokta Gauss) örnekle,
  Lagrange ile istasyona interpole et → sahte salınım yok. *Doğrulandı `test_struct_forces`:* konsol
  uçtan yük → M(s)=P(L−s) LİNEER (max_M_err 2.7e−9), Q(s)=P SABİT (Barlow ile err 5.9e−10), N≈0.
- **`anchor_force`**: solver'ı (nonlinear_solver.cpp anchor döngüsü) BİREBİR yansıtır:
  N=clamp(kk(U−U_p), −F_c, +F_t), kk=EA/L, U=Σg_i·u (sabit DOF dışlanır, solver gibi). U_p =
  NewtonResult.anchor_plastic[i] (elastik ⇒ 0). Çıktı `{N, yielded}`.
- **`geogrid_force_diagram`**: `geogrid::axial_return`'ü birebir yansıtır (tension-only + N_p); plastik
  durum 2 Gauss noktasında tanımlı → diyagram o 2 noktada (PLAXIS gerilme-noktaları gibi). ε_p =
  NewtonResult.geogrid_plastic[2·gi..].
- **`force_envelope(diag)`**: bir diyagramın max |N|,|Q|,|M| zarfı (wall_force_envelope ile uyumlu).

## 7. Anchor — tek-yönlü eksenel yay (PLAXIS MMM Eq 18-1; UYGULANDI)
```
N = (EA/L)·U                          (U = uzama; düğümler arası eksenel)
```
İki düğüm arası eksenel yay; yalnız öteleme DOF'ları (φ YOK). Rijitlik K=kk·g·gᵀ (kk=EA/L,
g=∂U/∂u=[−dir,+dir]); 4×4 (node-to-node) ya da 2×2 (fixed-end). Çözücüde ayrı döngü (zemin yolu
değişmez), elastik f=K·u_total. **PLAXIS Reference Manual:** node-to-node ve fixed-end AYNI eleman
teknolojisi (tek-yönlü yay):
- **node-to-node:** iki mesh düğümü arası (strut/iç destek). L = eşdeğer uzunluk (≤0 ⇒ geometrik
  mesafe; eşdeğer ≠ geometrik ise PLAXIS ölçek faktörü = doğrudan EA/L_eq).
- **fixed-end:** bir mesh düğümü + sabit uzak-uç (yön + eşdeğer uzunluk anchor özelliği). Zemin
  ankrajı (bağ bölgesi mesh dışında, sabit). Çözücüde node_b<0 ⇒ fixed_point.
Doğrulandı `test_anchor`: fixed-end yay N=EA/L·U≈F + u≈−FL/EA (yumuşak zeminde); node-to-node strut
duvar yanal sehimini ~600× azaltır, strut kuvveti ≈ uygulanan yük (denge).

### 7a. Anchor — elastoplastik (F_max,tens / F_max,comp; PLAXIS MMM §18.1, UYGULANDI)
PLAXIS ankrajı elastik-mükemmel-plastik tek-yönlü yay: eksenel kuvvet N elastik EA/L rijitlikle
gelişir ama **maksimum çekme F_max,tens ve maksimum basınç F_max,comp** ile sınırlıdır (her ikisi de
pozitif büyüklük; admissible N ∈ [−F_max,comp, +F_max,tens]). Akmadan sonra kuvvet platoda kalır,
**kalıcı uzama** (plastik elongasyon U_p) birikir — destek kapasitesi aşıldığında perde/zemin
serbestçe ilerler. Bu, desteklenmiş kazının kritik davranışıdır (struts/ankrajlar kapasitelidir).

**1D return mapping** (state = committed plastik elongasyon U_p^c; U = güncel toplam elongasyon):
```
N_tr = kk·(U − U_p^c)        kk = EA/L              (elastik öngörücü)
N_tr > F_t      :  N = F_t,  U_p = U − F_t/kk,  D_t = 0     (çekme akması)
N_tr < −F_c     :  N = −F_c, U_p = U + F_c/kk,  D_t = 0     (basınç akması)
aksi (elastik)  :  N = N_tr, U_p = U_p^c,       D_t = kk
```
F_t=F_max,tens, F_c=F_max,comp. F_max≤0 ⇒ sınırsız (saf elastik; eski davranış birebir korunur).
Teğet akmada 0 (mükemmel-plastik); zemin+perde rijitliği sistemi tekil-olmaktan kurtarır. Yoldan-
bağımlı: U_p committed/trial olarak çözücüde taşınır (Gauss state gibi; adım yakınsayınca commit).
Doğrulama `test_anchor`: strut'a F_max'ı aşan yük → N=F_max'ta plato, ekstra sehim, boşaltmada
kalıcı U_p (klasik elastoplastik yay).

## 8. Geogrid (tension-only eksenel membran; PLAXIS MMM §18.2, Eq 18-2; UYGULANDI)
Plate'in eğilmesiz (EI=0, kesme yok, φ-DOF yok) hâli — yalnız eksenel membran (donatı/geosentetik).
Mesh kenarına oturur (3-düğümlü eksenel eleman, 6 DOF = 3 düğüm × öteleme). N = EA·ε.
- **TENSION-ONLY:** geogrid yalnız çekme taşır; basınçta gevşer (slack) → N=0. Bu **tersinir**
  (nonlineer-elastik): tekrar gerilince hemen çekme taşır, kalıcı gevşeklik yok. N = max(0, EA·ε).
- **N_p (max çekme kuvveti, opsiyonel elastoplastik):** N N_p'ye ulaşınca akar → **kalıcı** plastik
  uzama (ε_p birikir). N_p≤0 ⇒ sınırsız (saf tension-only nonlineer-elastik, state'siz).

**Return mapping** (Gauss-noktası başına; state = committed plastik eksenel şekildeğiştirme ε_p^c):
```
N_tr = EA·(ε − ε_p^c)
N_tr ≥ N_p (N_p>0)  :  N = N_p,  ε_p = ε − N_p/EA,  D_t = 0    (çekme akması — KALICI)
N_tr ≤ 0            :  N = 0,    ε_p = ε_p^c,        D_t = 0    (basınç kesimi — TERSİNİR slack)
aksi (elastik)      :  N = N_tr, ε_p = ε_p^c,        D_t = EA
```
Basınç kesimi ε_p'yi DEĞİŞTİRMEZ (gevşeklik kalıcı değil, tersinir); yalnız N_p akması plastik
birikir. Eksenel B satırı = plate'in eksenel kısmı (Be, c,s gömülü); 2-nokta Gauss (kuadratik kenarda
ε lineer → tam). Montaj: f += ∫ Be^T N ds, K += ∫ Be^T D_t Be ds. ε_p committed/trial çözücüde taşınır.
Doğrulama `test_geogrid`: çekmede geogrid yükü taşır (sehim azalır), basınçta devre-dışı (N=0);
N_p akma platosu + boşaltmada kalıcı uzama.

## 9. Doğrulama planı (basitten karmaşığa)
1. **Konsol kiriş (plate, zeminsiz):** uçtan yük P → uç sehimi δ = PL³/(3EI) + PL/(kGA') (eğilme+kesme);
   saf eksenel NL/EA; rastgele yönelimde dönme-değişmezliği. *(Bu — `test_plate`.)*
2. **Basit kiriş / moment patch:** sabit moment → sabit eğrilik κ=M/EI.
3. **Plate-in-soil:** zemine gömülü perde, yanal yük → zemin-yapı etkileşimi (DOF paylaşımı doğru).
4. **Anchor:** elastoplastik yay; F_max sınırı.
5. **Desteklenmiş kazı (entegre):** perde + ankraj + kazı fazları → PLAXIS dersi.

## 10. Plate — elastoplastik M-N mafsalı (Mp/Np elması; PLAXIS MMM §18.3, UYGULANDI)

**Birincil kaynak (PLAXIS 2D Material Models Manual 2025.1, §18.3 "2D plates", s.222-223 — birebir):**
> "When the combination of a bending moment and an axial force occurs in a plate, then the actual
> bending moment or axial force at which plasticity occurs is lower than respectively Mp or Np.
> The relationship between Mp and Np is visualised in Figure 18-1. **The diamond shape represents
> the ultimate combination of forces for which plasticity will occur.** Force combinations inside
> the diamond will result in elastic deformations only."
> "Bending moments and axial forces are **calculated at the stress points** of the beam elements.
> When yield function is violated, **stresses are redistributed according to the theory of
> plasticity**, so that the maxima are complied with. This will result in **irreversible
> deformations**. Output of bending moments and axial forces is given in the nodes, which requires
> **extrapolation of the values at the stress points. Nodal forces are not checked** against the
> maximum forces."
> "By default the maximum moment is set to 1·10¹⁵ units if the material type is set to elastic."

Çıkarımlar: (a) akma yüzeyi **ELMAS** — köşeler (±N_p, 0) ve (0, ±M_p), kenarlar doğru parçası;
(b) kontrol **gerilme (Gauss) noktasında**, düğüm çıktısı gerilme-noktalarından EKSTRAPOLE (aşabilir
— PLAXIS de kontrol etmez); (c) Q (kesme) akma yüzeyinde DEĞİL (elastik kalır); (d) elastik = kapak
sınırsız. M-κ diyagramlı mod (Fig 18-2) bu kapsamda DEĞİL (ayrı iş kalemi).

### 10.1 Genelleştirilmiş-gerilme plastisitesi
Genelleştirilmiş gerilme s=(N,M), şekildeğiştirme e=(ε,κ), elastik D=diag(EA,EI), plastik durum
e_p=(ε_p,κ_p) (eğilme/eksenel Gauss noktası başına 2 skaler; kesme γ elastik, durumsuz):
```
N = EA·(ε − ε_p) ,   M = EI·(κ − κ_p) ,   Q = kGA'·γ            (kesme kapaksız)
f(s) = |N|/N_p + |M|/M_p − 1 ≤ 0        (1/N_p := 0 sınırsızsa; 1/M_p := 0 sınırsızsa)
ė_p = λ̇·∂f/∂s = λ̇·(sign(N)/N_p , sign(M)/M_p)                  (asosiye akış; Koiter köşede)
```
Mükemmel plastisite (pekleşme yok). Geri-Euler dönüşü = enerji-normunda (½ΔNᵀD⁻¹ΔN) dışbükey
çokgene EN-YAKIN-NOKTA izdüşümü → tek ve sürekli. Kapalı form (s_N=sgn(N_tr), s_M=sgn(M_tr)):

**Yüzey dönüşü** (elmasın (s_N,s_M) çeyreğindeki kenarı; a=(s_N/N_p, s_M/M_p)):
```
Δλ = f_tr / h ,   h = EA/N_p² + EI/M_p²
N = N_tr − Δλ·EA·s_N/N_p ,   M = M_tr − Δλ·EI·s_M/M_p
GEÇERLİ ⇔ s_N·N ≥ 0 VE s_M·M ≥ 0     (dönüş kenar parçasının üstünde kaldı)
```
**Köşe dönüşleri** (yalnız İKİ kapak da sonluyken; Koiter normal-konisi, SS kenar dönüşü deseni):
```
V_M = (0, s_M·M_p)  GEÇERLİ ⇔ (|M_tr| − M_p)·M_p/EI ≥ |N_tr|·N_p/EA
V_N = (s_N·N_p, 0)  GEÇERLİ ⇔ (|N_tr| − N_p)·N_p/EA ≥ |M_tr|·M_p/EI
```
(Türetim: Δe_p = D⁻¹(s_tr − V), koni {α·a₊ + β·a₋ : α,β ≥ 0}, a_± = (±1/N_p, s_M/M_p) →
s_M·Δκ_p·M_p ≥ |Δε_p|·N_p; V_N simetrik.) Geçerlilik-kademeli seçim: yüzey → V_M → V_N.
Plastik durum güncellemesi her dalda `e_p = e − D⁻¹s` (dönüş noktasından geri okunur).

**Tutarlı teğet** (akma yüzeyi DOĞRUSAL + mükemmel plastisite ⇒ tutarlı ≡ continuum):
```
elastik:  D_ep = diag(EA, EI)
yüzey:    D_ep = D − (D·a)(aᵀD)/(aᵀD·a)      (rank-1 indirgeme; N_p=∞ ⇒ diag(EA, 0))
köşe:     D_ep = 0                            (iki bileşen de sabitlendi — mafsal)
```
Tek-kapak hâlleri kendiliğinden: 1/N_p=0 → yüzey dönüşü tam M-kelepçesi (M=s_M·M_p, N elastik),
köşe yok; simetrik olarak 1/M_p=0. İki kapak da sınırsız → f=−1, saf elastik (eski yol).

### 10.2 Eleman entegrasyonu ve tesisat
- Eksenel+eğilme Gauss noktalarında (3-düğüm: 3-nokta; 5-düğüm: 5-nokta) return mapping; kesme
  AZALTILMIŞ noktalarında (2/4-nokta) her zaman elastik — §3'ün seçmeli azaltılmış şeması korunur.
- f_e = Σ w·J·(Beᵀ·N + Bkᵀ·M) + Σ w·J·kGA'·Bgᵀ·(Bg·u);  K_e aynı noktalarda D_ep ile.
- Durum: eleman başına [ε_p,κ_p]×Gauss — çözücüde committed/trial (anchor U_p deseni), fazlar arası
  `StructuralInit.plate_plastic/plate5_plastic` ile taşınır (Track 1a sözleşmesi). Return map
  committed durumun SAF fonksiyonu → line search güvenli.
- `props.plastic()` kapalıyken montaj ESKİ f=K·u yolunda (bit-birebir); elastoplastik eleman her
  zaman kuadratür yolunda (karma Gauss durumları için gerekli).
- Diyagram (PLAXIS kuralı birebir): N,M **Gauss değerlerinden** istasyonlara Lagrange açılımı
  (Q'nun Barlow deseniyle aynı); ekstrapole istasyon M_p'yi hafifçe AŞABİLİR — PLAXIS de düğüm
  değerini kontrol etmez (yukarıdaki alıntı; dürüst beyan).
- LİNEER Dinamik faz plate'i ELASTİK çözer (D6b: rapor da elastik zarf) — Mp/Np orada yalnız
  talep/kapasite uyarısı üretir (arayüz utilisation deseni); plastisite statik ailede + nonlineer
  dinamikte etkindir.

### 10.3 V&V (test_plate_plastic — hepsi ÖLÇÜLDÜ ve geçti)
1. Dönüş haritası taraması (81×81, 3 kapasite-katına dek): f ≤ 0 KESİN (worst 0.00e+00),
   idempotens, iç-nokta bit-kimliği, FD-teğet (bölge içleri parçalı-lineer → makine kesinliği),
   **bağımsız oracle** = elmas sınırının 4000-nokta örneklenmesiyle enerji-normu en-yakın-nokta:
   kapalı form ≡ CPP (fazla 0.00e+00). Tek-kapak halleri: kapaksız bileşen bit-birebir elastik,
   öbürü tam kelepçe, durum sürüklenmez.
2. Elastik-limit kimliği: eleman düzeyinde kuadratür f == K·u (1e-12); BVP'de kapaklı-ama-akmasız
   koşu kapaksız ikizle 1.3e-16'da özdeş; akma altında plastik durum BİT-SIFIR.
3. BVP akma BAŞLANGICI (yumuşak-zemin basit kirişi, test_plate_soil harness'i): kapalı form
   P_onset = 2·M_p/s_g (s_g = kritik GAUSS noktasının mesnet kolu — kontrol GERİLME NOKTASINDA
   olduğundan kapalı form Gauss konumunu içerir; PLAXIS ile aynı özellik). 0.9× → durum bit-sıfır,
   1.1× → mafsal (braket kapandı). 1.5× → kritik Gauss |M| TAM M_p'de (100.0000000000) + sehim
   orantılı elastiğin ≫10 katı (yumuşama). **ÖLÇÜLEN DERS: elastik zemin yatağı üzerindeki
   kirişte tek mafsal MEKANİZMA DEĞİLDİR — fazla yük zemine akar, çözüm yakınsar; Prandtl
   çökme-load_factor deseni burada uygulanmaz** (ilk tasarım öyleydi, ölçüm düzeltti).
4. M-N etkileşimi (BVP): eksenel ön-çekme T = N_p/2 + büyük eğilme → kritik Gauss ELMAS
   YÜZEYİNDE oturur (f = 0.000e+00) ve moment platosu ~M_p·(1 − N/N_p); N_p sınırsız kontrol
   aynı yüklerde tam M_p'ye çıkar → düşüş gerçekten etkileşimden.
5. Diyagram: kapaklı Gauss→istasyon Lagrange tepesi M_p·1.031 (ekstrapolasyon aşımı — PLAXIS
   düğüm kuralıyla tutarlı, beyanlı); BOŞ plastik durum = kapaksız ELASTİK rapor (D6b: lineer
   dinamik zarf kırpılmaz — mafsallı deplasmanda M ≫ M_p görünür, ölçüldü 464×).
6. tri15 (5-düğümlü plate): aynı onset braketi + TAM M_p doyması (100.0000000000).

**BİLİNEN GİZLİ SINIR (beyanlı, GUI'den erişilemez):** verilen-deplasman (prescribed ū) rampası
yalnız ZEMİN eleman döngüsüne ulaşır; yapısal eleman döngüleri sabit DOF'u okumaz (katkı 0).
Deformasyon ū bugün yalnız çekirdek/test yoludur (GUI'de deformasyon BC olarak yok); GUI'ye
eklenirse yapısal elemanlara toplam λ·ū taşınmalı (internal_forces.hpp Ramp notu).

İlgili: [[analysis-and-structural]], [[literature-review]], MC [[mohr-coulomb-formulation]],
HS [[hardening-soil-formulation]].
