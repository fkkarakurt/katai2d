# Initial stress and the K0 procedure — formulation (implementation specification)

Every geotechnical analysis is built on an **initial stress state** (geostatic). In PLAXIS the first
phase is produced by one of two methods: the **K0 procedure** (horizontal surface and layers) or
**gravity loading** (inclined surface). The working order is the same as elsewhere: lock the
mathematics against a reference first, then verify the smallest piece to round-off (see
`mohr-coulomb-formulation.md`, `effective-stress-formulation.md`).

Sources: Jaky (1944) for K0; Terzaghi & Peck for geostatic stress; PLAXIS Reference and Scientific
Manual (Bentley) for the K0 procedure and gravity loading. [[literature-review]].

## 0. Convention (consistent with the solver)

- **Tension-positive.** Compressive stress is negative. Compression increases with depth (as z
  decreases), so the vertical stress becomes more negative.
- Vertical effective stress (geostatic, horizontal surface, unsaturated):
  ```
  σ'_v(z) = −γ·(z_surf − z)            (γ = unit weight, z_surf = surface level)
  ```
  Below the water table the **buoyant** unit weight γ' = γ_sat − γ_w applies (effective stress; see
  `effective-stress-formulation.md`).

## 1. K0 — the coefficient of lateral earth pressure

The coefficient of lateral earth pressure at rest:
```
K0 = σ'_h / σ'_v
```
- **Jaky (1944):** for normally consolidated soil, `K0 = 1 − sin(φ')`. (φ' = 30° gives K0 = 0.5.)
- Overconsolidated: `K0 = (1 − sinφ')·OCR^(sinφ')` (Mayne & Kulhawy 1982), or a user value.
- K0 is **supplied by the user**, as in the PLAXIS K0 procedure; Jaky is the suggested default.

## 2. The K0 procedure — the initial stress field

With a horizontal surface and horizontal layers the initial stresses are **assigned directly**; no
equation is solved:
```
σ'_v  = −γ'·(z_surf − z)
σ'_h  = K0 · σ'_v        (both σ_xx and the out-of-plane σ_zz: both horizontal directions)
σ_xy  = 0
```
The plane-strain Gauss-point state is `stress = [K0·σ'_v, σ'_v, 0]` with `stress_zz = K0·σ'_v`.

**Key observation — equilibrium holds for any K0.** In horizontal layering this field is in *exact*
equilibrium with self-weight, independently of the value of K0:
```
Vertical:    ∂σ_yy/∂y + b_y = γ' + (−γ') = 0     (σ_yy = −γ'(z_surf−z), b_y = −γ')
Horizontal:  ∂σ_xx/∂x + ∂σ_xy/∂y = 0 + 0 = 0     (σ_xx depends on z only, σ_xy = 0)
```
This is why PLAXIS restricts the K0 procedure to a **horizontal surface**: on inclined or
non-layered geometry the x-derivative of σ_xx is not zero, so equilibrium is violated and plastic
points appear — and **gravity loading** is used instead.

## 3. Gravity loading — the "natural" elastic K0

If no pre-stress is assigned and gravity alone is applied to an elastic solve, the lateral stress
ratio obtained under laterally confined (oedometer) conditions is:
```
K0_elastic = ν' / (1 − ν')
```
(ε_xx = ε_zz = 0 gives σ'_h = (ν/(1−ν))σ'_v.) This differs from the **free** K0 of the K0 procedure:
gravity loading locks K0 to ν', whereas the K0 procedure allows a realistic K0 such as Jaky's to be
assigned. In KATAI, `test_stress` already verifies K0 = ν/(1−ν) to round-off.

## 4. Seeding the solver (the staged-construction hook)

The initial stress is handed to the nonlinear solver as the **committed initial Gauss state**. On the
first step Δε = 0, so σ_trial = σ_init and f_int = ∫Bᵀσ_init. If the external load is gravity and
σ_init is the horizontally layered K0 field, the residual is f_grav − ∫Bᵀσ_init ≈ 0, giving **zero
displacement** — a pre-stressed body is already in equilibrium under its own weight. This is exactly
the PLAXIS rule that the first phase after K0 must produce zero displacement, and it is the mechanism
staged construction requires in order to carry σ from phase to phase.

## 5. Verification (level 0 to 1)

- **K0 procedure equilibrium, to round-off:** a homogeneous column with σ_init set to the K0 field
  (K0 from Jaky), gravity applied, gives **displacement ≈ 0** and a committed stress equal to the K0
  field, for any K0. This proves both that the solver carries the pre-stress correctly and that the
  field is equilibrium-consistent.
- **Contrast with gravity loading:** no pre-stress, gravity solved, gives σ_h/σ_v = ν/(1−ν). This
  ties the two PLAXIS methods to each other.
- **Later, staged construction:** elements activated and deactivated; on excavation the internal
  force of the removed element is transferred to its neighbours; equilibrium is re-established in
  every phase and displacements can be reset.

## 6. K0 on an inclined surface — the equilibrium (nil) step

Per the PLAXIS Reference Manual, the K0 procedure is correct **only when the surface, the layer
boundaries and the water table are all horizontal**. Otherwise the column-overburden field leaves
**genuine out-of-balance forces**: principal directions must rotate and shear must develop on the
slope face. The PLAXIS remedy is a **plastic nil-step** — an otherwise empty calculation step that
resolves the imbalance — and the equilibrium reached can then be compared with the gravity-loading
result.

The KATAI implementation, currently in `build_problem.hpp`:

- In the K0 phase the base force `B = f_int(σ_K0) (+ interface σ_n0)` is held while external loads are
  ramped, so `residual(0) = 0` on every mesh: the seed *is* the answer, and in the horizontal case
  **u = 0 exactly**.
- **When the geometry is not horizontal** (`k0_nonlevel`: inclined surface, inclined water table, or a
  non-horizontal boundary between regions differing in γ or K0), the imbalance `d = f_body −
  f_int(seed)` is **ramped together with the external loads** to `target(1) = f`, the true
  equilibrium. The nil-step is thereby embedded in the K0 phase, and the displacement field shows the
  redistribution.
- **The trigger is geometric, not a threshold on |d|.** Quadrature residues — the r-weighted cubic
  integrand in axisymmetry, or a water-table kink crossing the interior of an element — can be the
  same order of magnitude as the genuine imbalance of a gentle slope. On horizontal geometry the seed
  is the intended answer, so a magnitude test would misclassify both cases.
- **Exact column integration at strata breaks** (`K0LayeredOptions::strata_breaks`): the integral is
  taken segment by segment at layer-boundary and water-table levels, so σ'_v is **exact** for
  piecewise-constant γ'. This preserves `f_int(σ_K0) ≡ f_gravity` to round-off in the horizontal
  multi-layer and water cases.
- **Union of boundary conditions at corners:** a boundary node takes the **union** of the constraints
  of every polygon edge it lies on. The previous equidistant tie-break, in which a free edge could
  override a fixed one, was a genuine hole in the support at a node shared by a layer interface and a
  side boundary; it is fixed.

Verification (`test_k0_slope`, through the GUI path): horizontal single-layer, two-layer and
horizontal water table give max|u| = 0 exactly with no nil-step; the Griffiths & Lane slope (stable,
c = 20) triggers a nil-step, develops max|σ_xy| ≈ 38 kPa, and σ'_v below the crest is within 10% of
both the overburden value and the gravity-loading value; the unstable slope (c = 1) produces an honest
collapse report at load factor 0.83; an inclined water table triggers the nil-step.

Sources: Jaky (1944) *J. Soc. Hungarian Arch. Eng.*; Mayne & Kulhawy (1982); PLAXIS Reference and
Scientific Manual (Bentley); ISSMGE HTC L5 "Practical considerations" (non-horizontal K0 requiring a
high-accuracy nil step); [[literature-review]].
