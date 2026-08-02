# third_party — vendored bağımlılıklar

Hepsi **permissive** lisanslı (ticari-güvenli, bkz `docs/DECISIONS.md` D7).
Drop-in: kaynaklar repoya kopyalandı, iç `.git` dizinleri silindi (T-ENV-4).
Ağsız, tekrarlanabilir build; tek build sistemi (kök CMake) hepsini derler.

| Klasör | Kütüphane | Sürüm (vendor tarihi 2026-06-02) | Lisans | Kullanım |
|--------|-----------|----------------------------------|--------|----------|
| `eigen/` | Eigen | dev/master (≈3.4.90, WORLD=3) | MPL-2 | yoğun/küçük matris cebir (header-only) |
| `glfw/`  | GLFW | 3.5.0 (dev) | zlib | pencere + GL context |
| `glad/`  | glad2 üretimi | GL **4.1 core** (D17) | MIT / public-domain | GL fonksiyon loader |
| `imgui/` | Dear ImGui | 1.92.9 WIP (**docking** dalı) | MIT | GUI chrome (D16) |

## Güncelleme / yeniden üretme
- **glad** (üretilmiş, ağ gerektirir):
  `python -m glad --api "gl:core=4.1" --out-path third_party/glad c`
- **eigen / glfw / imgui**: ilgili upstream'den shallow clone → `.git` sil.
  GLFW örnek/test/docs CMake'te kapalı; ImGui yalnızca core + glfw/opengl3 backend.

> Not: oneMKL **vendored değil** — sistemde kurulu (T-ENV-1) ve `LinearSolver`
> soyutlaması arkasından link'lenecek (Faz 0'da Eigen ile başlanıyor).
