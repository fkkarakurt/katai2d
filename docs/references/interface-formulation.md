# Interface (zemin-yapı ara yüzü) — Formülasyon (P2.4)

Desteklenmiş kazı dikeyinin son kritik parçası: perde/duvar ile zemin arasındaki **göreli kayma
ve ayrılma** (sliding/separation). Onsuz yapı zemine "yapışık" → fazla rijit, gerçek-dışı moment/
sehim. PLAXIS'in `R_inter` mukavemet-azaltma faktörü tam burada devreye girer.

**Kaynaklar (kilitli):**
- **Goodman, Taylor & Brekke (1968)** — sıfır-kalınlık (zero-thickness) joint elemanı; göreli
  deplasman ↔ traksiyon, normal & kayma rijitliği (k_n, k_s).
- **Day & Potts (1994)**, *IJNAMG* 18:689–708 — sıfır-kalınlık interface'in sayısal kararlılığı:
  yüksek k_n'in yarattığı gerilme salınımını önlemek için **Newton-Cotes (nodal) integrasyon**
  (Gauss değil) önerilir → düğüm-çiftleri ayrık yay gibi davranır, salınım yok.
- **PLAXIS 2D Reference & Material Models Manual** — bilineer **Coulomb** interface modeli; mukavemet
  azaltma `R_inter`; sanal kalınlık (virtual thickness) ile rijitlik; interface ν_i=0.45.

## 0. Konvansiyon
Çözücü **tension-pozitif** (tüm kod tabanıyla aynı). Interface = mesh kenarına oturan sıfır-kalınlık
hat elemanı (tri6 kenarı = 3 düğüm: 2 köşe + 1 orta). İki taraf: **zemin-tarafı** (mesh düğümleri,
taban DOF) ve **yapı-tarafı** (çakışık düğümler, `DofMap::add_extra_dof` ile 2'şer ek DOF ux,uy).
Yapı (plate/yük) yapı-tarafı DOF'larına bağlanır; interface iki tarafı k_n,k_s + Coulomb ile bağlar.
Yerel eksen: s = kenar teğeti (ŝ=(c,s)), n = kenar normali (n̂=(−s,c)).

## 1. Kinematik — göreli deplasman (jump)
Düğüm-çifti i için sıçrama [[u]]_i = u_(yapı,i) − u_(zemin,i). Gauss/Newton-Cotes noktasında:
```
Δu_s = Σ_i N_i ( c·[[u_x]]_i + s·[[u_y]]_i )      (teğetsel göreli kayma)
Δu_n = Σ_i N_i (−s·[[u_x]]_i + c·[[u_y]]_i )      (normal göreli açılma)
```
N_i(ξ) = 1B kuadratik şekil fonksiyonları (ξ=∓1 köşe, ξ=0 orta). B (2×12) bu ilişkiyi 12 DOF'a
(6 zemin + 6 yapı) bağlar; rotasyon (c,s) B'ye gömülü.

## 2. Bünye — bilineer Coulomb (PLAXIS interface, tension-pozitif)
Yerel traksiyon t=(τ, σ_n). Elastik:
```
τ   = k_s · Δu_s^e          σ_n = k_n · Δu_n          (k_s kayma, k_n normal rijitlik [kN/m³])
```
**Normal tension cut-off:** σ_n = k_n·Δu_n ama σ_n ≤ σ_t (σ_t = çekme dayanımı, varsayılan 0).
Açılma (Δu_n>0 → σ_n>0) σ_t'yi aşarsa interface ayrılır (gap): σ_n = σ_t, normal rijitlik düşer.
**Coulomb kayma kriteri** (σ_n tension-pozitif; basınç σ_n<0 ⇒ −σ_n tanφ_i>0 kapasite ekler):
```
τ_max = c_i − σ_n · tan φ_i        (≥0'a clamp)
|τ| < τ_max  ⇒ elastik (yapışık);   |τ| = τ_max ⇒ kayma (plastik slip)
```
**Return mapping** (Newton-Cotes düğüm-çifti başına; state = committed plastik kayma Δu_s^p,c):
```
σ_n = clamp(k_n·Δu_n, −∞, σ_t)                         (normal; tension cut-off)
τ_tr = k_s·(Δu_s − Δu_s^p,c) ;   τ_max = max(0, c_i − σ_n tanφ_i)
|τ_tr| ≤ τ_max :  τ = τ_tr,            Δu_s^p = Δu_s^p,c,        D_s = k_s   (elastik)
aksi          :  τ = sign(τ_tr)·τ_max, Δu_s^p = Δu_s − τ/k_s,    D_s = 0     (kayma — KALICI)
```
Mükemmel-plastik (ψ_i=0, PLAXIS varsayılanı R_inter<1): dilatasyon yok ⇒ kayma normal'i etkilemez,
teğet ayrık (D=diag(D_s,k_n)) — simetrik, SPD-dostu. Δu_s^p committed/trial taşınır (gauss_states gibi).

## 3. Mukavemet azaltma — R_inter (PLAXIS)
Bitişik zeminin MC parametrelerinden:
```
tan φ_i = R_inter · tan φ_soil ,   c_i = R_inter · c_soil ,   ψ_i = 0 (R_inter<1; aksi ψ_soil)
G_i     = R_inter² · G_soil  (≤ G_soil)
```
**Rijitlik (sanal kalınlık t_i):** k_s = G_i / t_i, k_n = E_oed,i / t_i; ν_i=0.45 sabit (sayısal),
E_oed,i = 2G_i(1−ν_i)/(1−2ν_i); t_i = (virtual thickness factor δ≈0.1)·ortalama-eleman-boyu. Bizim
elemanda k_n,k_s GİRDİ (caller PLAXIS bağıntısıyla hesaplar); helper `interface_stiffness` sağlar.

## 4. Eleman rijitliği + integrasyon (Newton-Cotes)
```
K = ∫_L Bᵀ D B ds ,  f = ∫_L Bᵀ t ds ,  D = diag(D_s, k_n)
```
**Newton-Cotes (nodal) integrasyon** (Day & Potts 1994): noktalar 3 düğümde (ξ=−1,0,+1), ağırlıklar
w=(1/3, 4/3, 1/3), ds=J dξ, J=L/2 (düz kenar). Her düğümde yalnız N_i=1 → düğüm-çiftleri **ayrık**
(zemin-yapı yayı k·w_i·J) → yüksek-k_n salınımı yok (Gauss'un sorunu). Coulomb her çiftte bağımsız.

## 5. Doğrulama planı (basitten karmaşığa)
1. **Elastik kayma/normal:** sabit zemin tarafı, yapı tarafına normal basınç + kayma yükü (kapasite
   altı) → Δu_s=τ/k_s, Δu_n=σ_n/k_n round-off. *(test_interface)*
2. **Coulomb zarfı:** kayma yükünü rampala → limit = τ_max·L = (c_i − σ_n tanφ_i)·L; σ_n değiştir →
   doğrusal sürtünme doğrusu (c_i kesişim, tanφ_i eğim). *(test_interface)*
3. **Perde-zemin ara yüzü:** R_inter ile perde sehimi/moment "yapışık"a göre artar → tie-back kazı
   benchmark'ında PLAXIS ile karşılaştır (sonraki).

## 6. Başlangıç gerilmesi — K0 install (wished-in-place perde)
Gömülü perde (bariyer) split mesh'te modellenir: perde hattı boyunca düğümler ÇİFTLENİR (sol/sağ ayrık),
perde plate'i bağımsız ek-DOF'larda, her tarafa bir interface. **SORUN:** K0 jeostatik durumda yatay
gerilme σ_h=K0·σ'_v perde hattında SÜREKSİZdir; interface deplasman-türevli (σ_n=k_n·Δu_n) ve sıfır
normal gerilmeyle başlar → perde "kurulum"da büyük SAHTE hareketle yerleşir (kazısız analiz = kazılı
analiz aynı sehim; FE artefaktı). PLAXIS K0 prosedürü interface gerilmesini bitişik zemin gerilmesinden
ALIR; biz bunu **başlangıç normal gerilmesi σ_n0** ile yaparız (tension-poz, basınç<0):
```
σ_n = σ_n0 + k_n·Δu_n      (tension cut-off + Coulomb zarfı TOPLAM σ_n üzerinde)
σ_n0 = K0·σ'_v             (her Newton-Cotes noktasında, derinliğe göre; σ'_v=−γ'(z_surf−y))
```
Δu_n=0'da (u=0) σ_n=σ_n0 → interface K0 yatay gerilmesini taşır → wished-in-place perde u=0'da ÖZ-DENGEDE
(sahte kurulum hareketi YOK). **NORMAL YÖNÜ KRİTİK:** sol ve sağ interface aynı geometrik kenardadır
(çiftlenmiş düğümler çakışık) → edge_frame ikisine de aynı (c,s)→aynı normal verir. Tek formül σ_n0=K0σ'_v
ile her ikisinin basınç olması için **sol interface'in düğüm sırası (A↔B) ters çevrilir** (build_embedded_wall),
böylece normali sol zeminden DIŞA bakar → her iki tarafta basınç σ_n<0 (tension cut-off doğru, seed kırpılmaz).
Aksi halde (aynı sıra) sol interface σ_n0>0 (tension) çıkar → cut-off seed'i sıfırlar. **Staged release**
(kazı): baseline B = f_int(u=0) = ∫_active Bᵀσ0 (zemin) + Σ interface σ_n0 kuvvetleri; target(λ)=B+λ(grav_active−B)
→ λ=0'da residual 0 (seeded denge), λ=1'de tam kazı boşalması. **Doğrulama (test_wall_k0_excavation):**
analitik residual u=0'da 4e-12; solver max|perde ux|=0 (install); kazıyla esnek bariyer kazıya sehim
(kazı-kaynaklı, ≫ install). Kaynak: PLAXIS K0 procedure → interface initial stress. **(P2.4, commit 9e786b3.)**

## 7. GENEL interface — keyfi yönelim + zemin-içi (standalone) (commit'lenecek)
§6'daki perde interface'i YALNIZ dikey perde (tri6) içindi. PLAXIS interface'i HER doğrultuda, her yapı
hattında VEYA zemin-içi serbest slip yüzeyi olabilir. Bunu **genel** kılan üç parça (validated §1-6
çekirdeğini değiştirmeden):
1. **`split_mesh_at_segment`** (staged_construction.hpp) — `split_mesh_at_wall`'un genellemesi: keyfi
   segment (ax,ay)→(bx,by) boyunca, segmentin POZİTİF tarafındaki (işaretli mesafe d=(p−a)·n̂≥0) düğümler
   orijinal, NEGATİF taraf ikiz (yeni) düğümlere bağlanır. Segment mesh kenarlarına oturmalı (yapısal/
   interface hattı PSLG kısıtı → düğümler hizalı). **KRİTİK HİZA DERSİ:** seam köşe düğümüyle BAŞLAMALI
   (köşe,orta,köşe… tri6 / köşe,q,orta,q,köşe… tri15); alt köşe atlanıp orta-düğümle başlanırsa interface
   üçlüleri kayar → kalıcı slip, kaynaklı limit sürekli ortamı kurtarmaz (test'te yakalandı).
2. **`build_soil_interface`** (general_interface.hpp) — seam'den soil-soil `InterfaceElement`(3)/
   `InterfaceElement5`(5): `soil_nodes`=pozitif-taraf orijinal, `struct_dof`=NEGATİF-taraf İKİZ düğümlerin
   **mesh DOF'ları** (ek-DOF değil; ikizler gerçek mesh düğümü). §1-4 montaj/Coulomb çekirdeği değişmeden
   yeniden kullanılır. Yapı-tarafı bir plate yerine karşı zemin → soil-soil joint.
3. **Yönelim-duyarlı K0 seed** — `σ_n0 = (K0·n_x² + n_y²)·σ'_v` (geostatik σ_xx=K0σ'_v, σ_yy=σ'_v,
   σ_xy=0 ⇒ σ_n=n̂ᵀσn̂). Dikey (n=±x) → K0σ'_v (§6 ile birebir); yatay (n=±y) → σ'_v.
**Wiring:** build_mesh Interface hatlarını PSLG kısıtı yapar; build_problem standalone `StructKind::Interface`'i
böler+kurar+seed'ler (tri6+tri15), fixed-boundary'ye ulaşan dup düğümlerin BC'sini yayar; GUI'de "Interface"
çizme aracı. **DOĞRULAMA: (test_general_interface)** kaynaklı interface (yüksek k + σ_t=∞) = bonded sürekli
ortam **%0.00** (mekanizma doğru, slip→0), yumuşak interface (sonlu k_s, σ_t=0) belirgin slip+açılma (var≠yok);
tri6+tri15. **(test_interface_gui)** GUI-yol: çıplak kısıt-hattı etkisi **%0.00** (mesh aynı) vs interface
etkisi **%28** → interface gerçekten monte ediliyor (eski davranışta düşürülüp %0 olurdu). **PLATE interface
(embedded wall) GENELLEŞTİRİLDİ (aynı commit):** dikey-tri6 dar durumu → her yönelim + tri6/tri15. build_problem
duvar ön-geçişi: dikey → doğrulanmış x-sütun split'i (`split_mesh_at_wall`, BİT-BİREBİR korundu); non-vertical →
`split_mesh_at_segment` (seam `s` arc-length ile anahtarlanır → builder'ın y-sıralaması doğru kalır); order'a göre
`build_embedded_wall`(6)/`build_embedded_wall5`(15); seed yönelim-duyarlı `(K0·nx²+ny²)·σ'_v` (dikey→K0·σ'_v birebir).
toe = derin uç (paylaşımlı). **DOĞRULAMA (test_interface_gui):** tri15 dikey duvar iface-bayrağı sonucu %29 değiştirir;
tri6 eğik duvar %29 değiştirir (eskiden bonded kalıp %0 olurdu); test_wall_k0_excavation (tri6 dikey) bit-birebir korundu.
**KALAN:** tri15 duvar M/Q/N kuvvet-diyagramı çıktısı (plate5::forces; solve doğru, yalnız diyagram); interface'in
fazlar-arası taşınması (split mesh sabit-geometri).

İlgili: [[structural-plate-formulation]] (plate/anchor/geogrid), [[mohr-coulomb-formulation]],
[[analysis-and-structural]] §5 (Goodman/Day&Potts), [[literature-review]].
