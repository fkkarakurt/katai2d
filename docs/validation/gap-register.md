# Gap register — "no example we cannot run"

Roadmap section 6.6: every case in the PLAXIS 2D official Validation Manual and Tutorial
Manual, plus the selected external verification sets, is either reproduced or explicitly
declared. Four statuses, and nothing without one:

- **reproduced** — run in this tree, deviation measured, pinned by a ctest regression or a
  recorded study;
- **reproducible — not yet run** — every needed capability exists; the case is queued work,
  not a limitation;
- **blocked** — a named missing capability; each such row is a queue item in roadmap
  section 8.4;
- **out of declared scope** — with the reason, so the boundary is a statement rather than
  an omission.

**Completion state, honestly:** this register currently holds the rows the existing record
supports. Enumerating the *complete* table of contents of the two official manuals against
it is the next §6.6 action; until that pass is done, absence of a case from this file means
*not yet triaged*, never *covered*. (Official manuals: PLAXIS 2D 2025.1 set; the Validation
Manual cases below cite the Version 8 case numbering used when they were reproduced.)

## 1. PLAXIS Validation Manual — cases triaged so far

| Case | Status | Evidence |
|---|---|---|
| §2.1 Smooth rigid strip footing, elastic (Giroud) | **reproduced** | KV-FND-001, `test_plaxis_validation`; +1.4% vs analytic, +0.8% vs PLAXIS's published number |
| §2.2 Strip load on Gibson soil (E = 299 z) | **reproduced** | KV-FND-002, `test_plaxis_validation`; −4.0% vs PLAXIS's finite-layer value, shared −6% bias vs half-space pinned |
| §2.3 Beam bending | reproducible — not yet run as the manual's case | two-sided analytic parity recorded: cantilever closed form 0.01% (Blevins), quartic plate 2.5e-13 (`test_plate`, `test_plate5`); the manual's exact configuration not yet mirrored |
| §2.4 Axisymmetric plate bending | **blocked — axisymmetric structural elements** | structural elements are plane-strain only today; capability queued (v1.x), see the comparison table's honest-gap row |
| §2.6 Updated mesh / large deformation | **out of declared scope** | small-strain formulation is a declared boundary (ROADMAP scope matrix); revisit only if the scope decision changes |
| §3.1 Circular footing bearing capacity (Cox) | **reproduced** | KV-FND-003, `test_plaxis_validation`; +3.9% vs the slip-line value at the 0.25 m mesh, psi = phi stated |
| §3.2 Strip footing on c(z) clay (Davis–Booker) | **reproduced** | KV-FND-004, `test_plaxis_validation`; +2.8% vs analytic (smooth case); the rough variant is queued |
| §3.3 Interface Coulomb slip (block) | reproducible — not yet run as the manual's case | interface Coulomb closed forms complete (`test_interface`, `test_general_interface`); the manual's block configuration not yet mirrored |
| §3.4 Cylindrical cavity expansion | reproducible — not yet run as the manual's case | fully plastic cylinder p = 2c ln(b/a) at 0.3% and Lamé at 1e-5 (`test_axisym_collapse`, `test_axisym_cylinder`); the manual's exact case not yet mirrored |
| §4.1 One-dimensional consolidation (Terzaghi) | **reproduced** (same analytic, both sides) | KV-CON-001, `test_consolidation` + GUI path `test_consolidation_gui`; U(Tv) within 0.03 |
| §4.2 Unconfined flow with a free surface (Dupuit) | reproducible — not yet run as the manual's case | Dupuit + Charny seepage-face + Darcy closed forms verified (`test_seepage`, `test_seepage_gui`) |
| §4.3 Confined flow under a wall (Harr) | reproducible — not yet run as the manual's case | under-dam flow-net and uplift closed forms verified (`test_seepage`) |

## 2. PLAXIS Tutorial Manual — cases triaged so far

| Case | Status | Evidence |
|---|---|---|
| Lesson 1 — circular footing on sand | **reproduced** (study) | `study_footing_plaxis`, `docs/validation/footing-benchmark-plaxis.md`; +3.2% on the published footing force |
| §17.8.5 dynamic base conventions | **reproduced** (conventions locked) | compliant-base implementation locks the Tutorial/Scientific-manual conventions verbatim; radiation oracle −0.02% (`test_compliant_base`, seismic-verification.md) |
| Excavation lessons (anchored/strutted walls) | reproducible — not yet run | all needed capabilities exist (staged construction, embedded walls, anchors, interfaces); queued as EXC-class matrix cases |

## 3. External verification sets — triaged so far

| Set | Status | Evidence |
|---|---|---|
| El Centro 1940 NS record | **reproduced** | KV-DYN-001, `test_real_record`; shipped digitisation PGA 0.31882 g vs published 0.319 g, spectrum inside the published band |
| SHAKE-class layered site response | **reproduced** (independent oracle) | `test_site_response_benchmark`; transfer-matrix solution sharing no code, deviations 0.06–0.18% (site-response-benchmark.md) |
| Sheet-pile wall in cohesive soil (published PLAXIS numbers, IJCRT 2024) | **reproduced** (study) | `study_wall_benchmark`, wall-benchmark-plaxis.md; tri15 M_max −3% |
| ACADS slope-stability set | reproducible — not yet run | phi-c reduction verified against Bishop-class values (`test_slope`); the ACADS cases themselves not yet mirrored |
| Schweiger Berlin excavation (multi-code) | reproducible — not yet run | HS-small, staged excavation, walls and anchors all exist; a large modelling effort, queued as the flagship EXC case |
| OpenSees cross-runs (dynamics, site response) | not yet run — planned | roadmap §6.5: runnable locally, upgrades the record from citation to repeatable cross-run |
| KRATOS GeoMechanics cross-runs (consolidation, staged excavation) | not yet run — planned | roadmap §6.5, same rationale |

## 4. Declared refusals that belong to the interaction register, not here

SoftSoilCreep × Consolidation and NonPorous × Flow are **refused by design with the
refusal itself tested** (R2). They are recorded here only to say where they live: the
schema-interaction register of roadmap §6 (`test_capability_coverage`, future work), not
the external-case register.

> Update discipline: a row changes status only with evidence (a test id, a study document,
> or a named capability), and the next full pass over the manuals' tables of contents may
> only add rows — it never removes one.
