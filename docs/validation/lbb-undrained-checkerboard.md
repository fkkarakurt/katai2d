# Equal-order u–p Biot elements: LBB / checkerboard characterization (undrained)

**Status:** characterized (was "unproven" in the SSI/water audit). Reproduce with
`cmake --build build/msvc-rwdi --target study_lbb_undrained` then run
`build/msvc-rwdi/bin/study_lbb_undrained.exe`. Source: `tests/study_lbb_undrained.cpp`.

## The concern

`analysis/consolidation.hpp` interpolates displacement and pore pressure with the **same** shape
functions (tri6 → P2–P2, tri15 → P4–P4). Equal-order u–p pairs violate the LBB / inf-sup condition,
so at the near-incompressible undrained limit they can admit spurious checkerboard pore-pressure
modes. The storage block `S = ∫(n/K_w) NᵀN` stabilizes the system only by `O(1/K_w)`, which vanishes
as the pore fluid becomes incompressible. All prior consolidation tests are 1‑D columns (which barely
excite the 2‑D mode) and skip the early-time boundary layer, so the 2‑D undrained pore field was
never checked.

## The experiment

A laterally-confined, base-fixed, **sealed (undrained)** 12×8 m saturated block (LE, E=10 MPa, ν=0.3),
with a central 4 m **strip** pressure (q=100 kPa) applied at t=0⁺. The partial-width load gives a 2‑D
non-uniform undrained pore field with an edge stress-concentration — where checkerboard appears.
Metric: `osc[n] = |p[n] − mean(p over element-neighbours)| / range` (range-normalized; smooth ≈ 0.02–0.05,
checkerboard ≈ O(0.1–1)). Discriminators: realistic vs soft `K_w/n`; a `dt` sweep; and the settlement field.

## Results

| element | probe | pore rms-osc | checker nodes | settlement rms-osc |
|---|---|---|---|---|
| **tri6** | t=0⁺ (dt→0), realistic K_w/n=6.7e6 | **0.154** | 554 | 0.006 |
| tri6 | t=0⁺, soft K_w/n=100 (control) | 0.028 | 30 | 0.006 |
| tri6 | dt=1e‑2, realistic K_w | 0.015 | 7 | — |
| tri6 | dt=1.0, realistic K_w | 0.006 | 0 | — |
| **tri15** | t=0⁺ (dt→0), realistic K_w | **0.024** | 56 | 0.010 |
| tri15 | t=0⁺, soft K_w (control) | 0.019 | 38 | — |

## Conclusions

1. **The LBB checkerboard is REAL but confined to tri6 at the near-pure-undrained instant.** At
   realistic (incompressible) K_w the tri6 pore field oscillates node-to-node (rms-osc 0.154, 554
   nodes, pore swinging −175…+113 kPa under a compressive strip); making the fluid compressible drops
   it 5.5×. That K_w-dependence is the LBB signature.
2. **Settlements are NOT affected.** The displacement field stays smooth (rms-osc 0.006, identical to
   the compressible control). The spurious mode lives in the pressure null-space; the engineering
   output (settlement/consolidation curve) is unpolluted. The validated 1‑D Terzaghi U–T_v and final
   settlement results are unaffected.
3. **Any normal time step damps it.** The backward-Euler diffusion `Δt·H` smooths high-frequency
   pore modes: at dt≥~1e‑2 (here) the checkerboard is essentially gone (0 checker nodes at dt≥1).
   It only appears when the first step is pathologically small (short duration ÷ many steps on tri6).
4. **tri15 is practically immune** (rms-osc 0.024, K_w ratio 1.2×) — higher-order pore interpolation
   does not excite the mode here.

## Practical guidance / mitigations

- Trust **settlements and consolidation curves** on any element — they are LBB-insensitive.
- For clean **early-time pore-pressure fields** on tri6, use an adequate first time step, or switch to
  **tri15**. A full inf-sup-stable fix (reduced-order pore, e.g. Taylor–Hood P2–P1) is a larger solver
  change and is **not warranted** by the practical impact measured here (settlements correct; pore
  smooths under normal dt).
- Defensive UI (DONE, audit P2c): the phase editor warns when a Consolidation / fully-coupled phase's
  step Δt = duration/steps is below the Vermeer–Verruijt Δt_crit = h²γw/(η k)(1/Eoed + n/Kw) (see
  `consolidation-formulation.md` §4), which is exactly the regime that reveals the checkerboard. The
  formula is `katai::core::consolidation_critical_dt` (η = 40 tri6 / 80 tri15, h = mean element size),
  the phase-level check is `katai::app::consolidation_step_warning`. Calibrated against this study
  (`test_consolidation` `test_critical_dt`): Δt_crit ≈ 4.65e‑3 day for the 0.5 m tri6 mesh here, which
  brackets the measured onset (checkerboard severe at dt=1e‑8, gone by dt≈1).
