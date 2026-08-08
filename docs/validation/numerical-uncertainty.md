# Numerical uncertainty: how far is a KATAI number from the exact solution of its own equations?

Every benchmark in this repository has, until now, been reported at one mesh density and one
solver tolerance — "+1.4% against Davis & Booker on the file's own 0.5 m tri15 mesh". That
sentence hides a question a reviewer is entitled to ask: **how much of that 1.4% is the model,
and how much is the mesh?** A deviation that shrinks under refinement is discretisation; one
that does not is the formulation, the domain size, or the reference itself. Until the two are
separated, neither can be defended.

This document records how KATAI separates them, where the procedure comes from, and what it has
measured so far.

## 1. Where the procedure comes from — and where it does not

The geotechnical finite element programs this one is measured against do not define a
discretisation-error estimator. PLAXIS, MIDAS GTS NX and GEO5 all provide mesh refinement
*controls* (global and local coarseness factors, refinement around structures and load lines)
and qualitative guidance — finer is more accurate, balance accuracy against run time — but none
publishes a procedure for putting a number on the discretisation error of a result. PLAXIS's
*Scientific Manual* chapter 9 does define a family of convergence criteria, and it is worth
being precise about what they govern: they decide when an **iterative** (Newton) loop may stop,
which is a different question from how far the converged answer of a given mesh sits from the
exact solution of the same equations. KATAI's parity work against chapter 9 is a separate matter
(see the convergence-criteria note).

The quantitative procedure therefore comes from the verification literature, where it is
standardised:

| Source | What is taken from it |
|---|---|
| P. J. Roache (1994), "Perspective: A Method for Uniform Reporting of Grid Refinement Studies", *ASME J. Fluids Eng.* **116**(3):405–413 | The Grid Convergence Index and its factor of safety: **1.25** when three or more grids *establish* the observed order, **3.0** when the order is *assumed* rather than measured. |
| I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, P. E. Raad (2008), "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications", *ASME J. Fluids Eng.* **130**(7):078001 | The step-by-step procedure adopted as that journal's editorial policy: the representative cell size, the fixed-point form of the observed order that admits **unequal** refinement ratios, the extrapolation, the error measures, and the recommendation that the refinement factor be at least **1.3**. |
| ASME V&V 20-2009 (R2016), *Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer*; ASME V&V 10-2006, *Guide for Verification and Validation in Computational Solid Mechanics* | The framework and its vocabulary. **Code verification** asks whether the code solves the equations it claims to; **solution verification** estimates the numerical accuracy of one particular calculation. This document is solution verification. V&V 10 is the companion whose scope — computational solid mechanics — actually contains this program. |

## 2. The equations, as implemented

Implemented in `kernel/math/include/katai/math/grid_convergence.hpp`. Three solutions of the
same quantity φ on three systematically refined meshes, subscript 1 = finest:

- **Representative cell size**, from the mesh that was actually built, not from the requested
  element size: `h = sqrt(total area / element count)` (Celik et al., step 1). Refinement
  factors `r21 = h2/h1`, `r32 = h3/h2`.
- **Observed order** `p`, by fixed-point iteration:
  `p = |ln|ε32/ε21| + q(p)| / ln(r21)`, `q(p) = ln((r21^p − s)/(r32^p − s))`,
  `s = sign(ε32/ε21)`, with `ε21 = φ2 − φ1`, `ε32 = φ3 − φ2`.
- **Richardson extrapolation**: `φ_ext = (r21^p φ1 − φ2)/(r21^p − 1)`.
- **Errors**: approximate `e_a = |(φ1 − φ2)/φ1|`; extrapolated `e_ext = |(φ_ext − φ1)/φ_ext|`.
- **Grid Convergence Index**: `GCI_fine = Fs · e_a / (r21^p − 1)`.
- **Convergence condition** from `R = ε21/ε32`: monotonic convergence (0 < R < 1), monotonic
  divergence (R ≥ 1), oscillatory convergence (−1 < R < 0), oscillatory divergence (R ≤ −1).

These were re-derived before they were coded, because the secondary literature restating Celik
et al. disagrees with itself: one widely-read version prints the extrapolated value as
`r^p (φ1 − φ2)/(r^p − 1)`, which is not Richardson extrapolation. Writing
`φ_i = φ_exact + C h_i^p` and eliminating `C` gives the forms above, and
`tests/test_grid_convergence.cpp` (**KV-NUM-004**) pins the implementation against triplets
manufactured from that closed form: the observed order and the extrapolated value come back to
1e-14 relative, at equal and at unequal refinement ratios, and the identity
`GCI_fine = Fs × (true relative error of the fine mesh)` holds exactly for power-law data.

## 3. What is a policy, and not a citation

Two decisions are KATAI's own. They are written out because a reader must be able to disagree
with them:

1. **When the observed order is not believed.** Roache provides both factors of safety; Celik
   et al. warn that an order far from the formal one means the triplet is not in the asymptotic
   range. Neither states a numeric acceptance window. KATAI accepts the observed order only
   when the triplet converges monotonically, `0.5 ≤ p ≤ 3.0`, and the asymptotic-range check is
   within 10% of 1. Otherwise the band is recomputed with an **assumed order p = 2 and
   Fs = 3.0**, and the extrapolated value is not quoted. The lower bound is not arbitrary: the
   amplification `1/(r^p − 1)` diverges as `p → 0`, so a small observed order turns mesh-to-mesh
   irregularity into a large, confident-looking correction — which is exactly what happened in
   the study below.
2. **A floor the evidence sets.** Outside the asymptotic range the reported band is widened, if
   necessary, to the spread across the three meshes. There the three values are not a hierarchy
   of errors crowned by the finest; they are three answers whose differences are dominated by
   the mesher's irregularity, and `e_a` — which sees only the finest pair — can happen to be the
   smallest of them. A band narrower than the range the family actually produced would claim a
   reproducibility the data denies. **Inside** the asymptotic range the same floor would be
   wrong, and it is not applied: the coarse mesh is expected to be further off, and its error is
   not the fine mesh's uncertainty.

## 4. First case measured: KV-FND-008, uniform strip load on an elastic half-plane

`tests/test_mesh_convergence.cpp` (**KV-NUM-005**), gated in the suite. The checked-in case file
is solved at three densities — only `mesh.elem_size` changes; geometry, material, load, element
order and the mesher's own local refinement rules are the file's — and the quantities are read
at exact probe points with the element shape functions rather than at the nearest node, because
the nearest node moves with every mesh and that jitter would be indistinguishable from
discretisation error.

| Requested elem_size | Elements | Nodes | Measured h |
|---|---|---|---|
| 0.700 m | 4 827 | 9 860 | 0.4071 m |
| 1.050 m | 2 184 | 4 497 | 0.6052 m |
| 1.575 m | 1 012 | 2 127 | 0.8891 m |

Refinement factors r21 = 1.487, r32 = 1.469 — both above the recommended 1.3.

| Quantity | Fine mesh | Observed p | Reported band | Basis | Closed form | Model deviation |
|---|---|---|---|---|---|---|
| σ_z at 2 m depth | −81.9428 kPa | 0.238 (monotonic) | **±0.496%** | assumed order, Fs = 3 | −81.8310 kPa | −0.137% — **inside** the band |
| σ_z at 4 m depth | −55.2331 kPa | 3.014 (oscillatory) | **±0.109%** | assumed order, Fs = 3, widened to the spread | −54.9815 kPa | −0.458% — **outside** the band |
| Settlement under the strip centre | −0.0212054 m | 6.98 (oscillatory) | **±0.107%** | assumed order, Fs = 3, widened to the spread | (none: a strip on a half-plane has unbounded surface settlement) | — |

Three things are worth stating plainly, and the first two are findings about the *study*, not
about the runs:

- **These triplets are not in the asymptotic range.** The observed orders — 0.24, 3.01, 6.98 —
  are not orders any element in this program can deliver. The reason is structural: KATAI's
  Ruppert mesher rebuilds the mesh from scratch at every density, so successive meshes are
  **not nested**, and at practical densities the mesh-to-mesh irregularity is the same size as
  the h^p trend it is meant to reveal. The GCI literature derives the method for uniform,
  systematic refinement and says as much; this is that caveat, measured on this code.
- **Richardson extrapolation would have made the 2 m answer worse.** With p = 0.238 the
  extrapolated σ_z is −80.28 kPa, i.e. +1.89% from the closed form, against −0.14% for the fine
  mesh the extrapolation was supposed to improve. The policy in §3 is what stopped that number
  from being published, and this row is why the policy exists.
- **At 4 m depth the model deviation is larger than the numerical band.** 0.458% against
  ±0.109%: the remaining deviation is *not* discretisation. The obvious candidate is the finite
  40 × 20 m domain standing in for a half-plane, which is a modelling choice, not a mesh
  property — and separating those two is the entire purpose of this exercise.

## 5. What this does not yet cover

The register is deliberately explicit about its own gaps, since an absent row must never read as
a passed one:

- **Fifteen of the sixteen corpus cases have no sweep yet.** KV-FND-008 is the first.
- **Tolerance independence is not measured.** Each case should also be re-run at solver
  tolerances 1e-2 / 1e-4 / 1e-6, the more so because the hardening-soil family currently runs at
  a 1% force residual — the reported quantity must move by less than the declared band, or the
  band is not a band.
- **Two open questions in the record are still open**: the slope factor of safety falling to
  0.872 on the finest mesh (ψ = 0, non-associated flow), and the −39% deep sheet-pile wall row,
  whose domain-size diagnosis is done but whose closure — a graded mesh over a doubled domain —
  is not.
- **Nested refinement is not available.** The most direct fix for §4's first finding is to
  refine a mesh rather than rebuild it, so that successive meshes share nodes. That is a mesher
  capability this program does not have, and until it does, the honest reading of these bands is
  the one given: an observational spread with an assumed order, not an asymptotic estimate.
