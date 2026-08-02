# GUI Tasarımı — PLAXIS 2D temelli, modern Dear ImGui (Faz GUI)

KATAI 2D arayüzü Dear ImGui (docking) + GLFW + OpenGL 4.1 üzerine kurulur. Amaç: ImGui'nin demode
varsayılan görünümünden kurtulup **modern, profesyonel** bir CAD/geoteknik arayüzü; **UTF-8/Türkçe**
tam destek; çalışma akışı **PLAXIS 2D Input program** mantığına dayalı. Kaynak: PLAXIS 2D 2025.1
Reference Manual §3 (Input program), §3.6 (Explorers), §7 (Meshing/calculation), §8 (Output).

## 1. Modernleştirme — PLAXIS-tarzı AÇIK tema, TAM İNGİLİZCE (uluslararası)
- **Dil:** Arayüz **tamamen İngilizce** (uluslararası program). Türkçe yalnız iç dokümanlarda.
- **Tema:** PLAXIS gibi **AÇIK (light) profesyonel tema** — açık-gri paneller, beyaz girdiler, koyu
  metin, PLAXIS-mavisi vurgu/seçim, ince kenarlıklar, Windows-tarzı düşük yuvarlatma (~3px).
  `apply_plaxis_style()`. (Koyu tema PLAXIS'ten uzaktı; açık tema esas.)
- **Font:** sistem Segoe UI ~18px (ImGui 1.92 dinamik atlas). UTF-8 destekli (σ/τ + Latin Extended).
- **Yerleşim:** docking + DockSpaceOverViewport(PassthruCentralNode) → merkez şeffaf "drawing area",
  GL sahnesi merkeze scissor'lı çizilir; paneller kenarlara doklanır (DockBuilder varsayılan PLAXIS
  yerleşimi). **Drawing area: beyaz zemin + soluk grid + üst/sol CETVELLER (rulers, kamera-eşlemeli).**

## 2. PLAXIS 2D çalışma akışı (Input program modları — Ref §3.4)
PLAXIS 5 modda çalışır; KATAI mod çubuğu (üst sekme/toolbar) bunları yansıtır:
1. **Soil (Zemin)** — borehole'lar, zemin katmanları, su koşulları, başlangıç koşulları (Ref §4).
2. **Structures (Yapılar)** — geometri, yükler, yapısal elemanlar (plate/anchor/geogrid/interface/
   embedded beam), hidrolik koşullar (Ref §5).
3. **Mesh (Ağ)** — mesh üretimi, yerel iyileştirme (Ref §7.1).
4. **Flow conditions (Akış)** — su seviyeleri, akış sınır koşulları (Ref §7.10).
5. **Staged construction (Kademeli inşaat)** — fazlar, eleman aktif/pasif, hesap (Ref §7.3-7.5).

## 3. Paneller (PLAXIS Explorers — Ref §3.6)
- **Model explorer** (sol): geometri + yapısal + zemin nesneleri ağacı.
- **Selection explorer / Properties** (sağ): seçili nesnenin özellikleri (malzeme, geometri).
- **Phases explorer** (Staged modunda): hesap fazları.
- **Viewport** (merkez): model + mesh + sonuç (deplasman/gerilme kontur, deforme mesh).
- **Log / Command** (alt): hesap ilerlemesi, mesajlar.
- **Menü çubuğu**: File, Edit, Soil, Structures, Mesh, Flow, Phases, Options, Help (Ref §3.3).

## 4. Output (sonuç) — Ref §8/§9
Deforme mesh, deplasman/gerilme kontur, yapısal kuvvetler (N/Q/M), load-displacement eğrisi.
Mevcut FieldView (GL kontur) çekirdeği bunun temeli.

## 5. Aşamalı GUI yol haritası
1. **Modern kabuk** [bu adım]: tema + UTF-8 font + docking layout + mod çubuğu + paneller (iskelet).
2. Soil mode: katman/malzeme tanımı + viewport çizim.
3. Structures mode: geometri/yük/yapısal eleman editörü.
4. Mesh mode: mesher entegrasyonu + görselleştirme.
5. Staged + hesap: faz tanımı + solve + sonuç görselleştirme (mevcut FieldView).
6. Output: kontur/eğri/yapısal-kuvvet görünümleri.

İlgili: [[plaxis-gap-analysis]], `docs/ARCHITECTURE.md`; render `katai_render`.
