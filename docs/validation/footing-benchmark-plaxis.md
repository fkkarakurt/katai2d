# PLAXIS 2D Tutorial Lesson 1 — Circular Footing on Sand (CAPSTONE)

Official PLAXIS 2D documentation example, reproduced from scratch and compared to the documented
PLAXIS result. The capstone integrated validation of KATAI's core capabilities.

**Source:** PLAXIS 2D 2025.1 Tutorial Manual, Lesson 1 "Settlement of a circular footing on sand",
Case A (rigid footing). Reproduced by `tests/study_footing_plaxis.cpp` (EXCLUDE_FROM_ALL).

## Problem
Axisymmetric. Sand layer 4 m thick over a rigid base; model radius 5 m; circular footing radius 1 m.
- **Soil (Mohr-Coulomb, drained):** γ_unsat=17, γ_sat=20 kN/m³; E'=13 MPa, ν'=0.3, c'=1 kPa, φ'=30°, ψ=0°.
- **Water table** at y=2 m. **K0 procedure** (K0 = 1−sin30 = 0.5) for the initial stress.
- **Phase 1 (rigid footing):** a prescribed uniform vertical settlement of 0.05 m over r ≤ 1 m.

**KATAI model:** axisymmetric tri6 + Mohr-Coulomb; two-layer K0 effective-stress seed (dry above the
water table, buoyant γ'=10 below); staged-release baseline constant_force = the K0 internal force
(holds the geostatic state, residual 0 at u=0); the ramped prescribed footing settlement develops the
reaction. Total footing force = (final − initial) vertical nodal force at the footing nodes × 2π.

## Result

| Criterion | KATAI | PLAXIS | |
|---|---|---|---|
| **Total footing reaction at 0.05 m** | **606.9 kN** | **588 kN** | **+3.2%** |

| Performance | value |
|---|---|
| Mesh | 640 elements, 1353 nodes, 2551 DOFs (axisym tri6) |
| Convergence | converged, load factor 1.0 |
| Iterations / time | 1174 total Newton iters, 26.4 s |

## Assessment
- **Accuracy +3.2% (< 5%)** against the official PLAXIS tutorial value, integrating axisymmetry +
  Mohr-Coulomb plasticity + K0 initial stress (with water table) + prescribed displacement in one
  boundary-value problem driven to near the bearing capacity. The residual is within the
  mesh-dependence band (PLAXIS quotes "about 588 kN"); a footing-refined mesh would tighten it.
- **Performance:** 26 s for a 2551-DOF axisymmetric MC analysis pushed to 0.05 m (near collapse).
  The high iteration count is the non-associated (ψ=0) Mohr-Coulomb cut-back near the limit load;
  arc-length control and a consistent non-associated tangent (gap-analysis Phase B) would cut this
  substantially. (An associated run, ψ=φ, converges far faster but is not the PLAXIS input.)

**Conclusion:** KATAI reproduces an official PLAXIS 2D documentation example from scratch to +3.2%,
exercising its core capabilities together — the capstone confirmation that the backend is accurate
and PLAXIS-consistent for a real, integrated geotechnical analysis.
