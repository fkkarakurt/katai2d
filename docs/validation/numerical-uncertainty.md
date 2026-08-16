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
   when the triplet converges monotonically, `0.5 ≤ p ≤ 4.0`, and the asymptotic-range check is
   within 10% of 1. Otherwise the band is recomputed with an **assumed order p = 2 and
   Fs = 3.0**, and the extrapolated value is not quoted. Neither bound is arbitrary. Below 0.5
   the amplification `1/(r^p − 1)` diverges as `p → 0`, so a small observed order turns
   mesh-to-mesh irregularity into a large, confident-looking correction — which is exactly what
   happened in the first version of the study below. Above 4 nothing in a 6-node quadratic
   triangle can be responsible: its displacement converges at O(h³) in L2 with nodal values
   capable of more, and its recovered stress at O(h²) with recovery superconvergent to about
   O(h³) (Zienkiewicz and Zhu 1992). A study on the 15-node quartic element must raise both.
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
is solved on a **nested** family: the mesher is visited once, and the two finer meshes are made
by splitting every triangle at its edge midpoints (`katai::mesh::refine_uniform`, verified in
**KV-NUM-006**). Quantities are read at exact probe points with the element shape functions
rather than at the nearest node, because the nearest node moves with every mesh and that jitter
would be indistinguishable from discretisation error.

| Level | Elements | Nodes | Measured h |
|---|---|---|---|
| coarse (mesher, elem_size 1.4 m) | 1 190 | 2 485 | 0.8199 m |
| refined once | 4 760 | 9 729 | 0.4100 m |
| refined twice | 19 040 | 38 497 | 0.2050 m |

Refinement factor **exactly 2** at both steps, with every angle of the coarse mesh preserved.

| Quantity | Fine mesh | Observed p | Reported band | Basis | Closed form | Model deviation |
|---|---|---|---|---|---|---|
| σ_z at 2 m depth | −81.8861 kPa | 3.65, monotonic | **±0.003%** | observed order, Fs = 1.25 | −81.8310 kPa | −0.067% — **outside** the band |
| σ_z at 4 m depth | −55.2003 kPa | 2.32, oscillatory | **±1.066%** | assumed order, Fs = 3, widened to the spread | −54.9815 kPa | −0.398% — inside the band |
| Settlement under the strip centre | −0.0212052 m | 3.90, monotonic | **±0.001%** | observed order, Fs = 1.25 | (none: a strip on a half-plane has unbounded surface settlement) | — |

At 2 m depth the triplet is in the asymptotic range, and the Richardson-extrapolated stress —
what this formulation would give on an infinitely fine mesh — is **−81.8839 kPa, 0.06% from the
closed form**, closer than any of the three runs.

Three things follow, and the first is the result this whole exercise was built to produce:

- **The 2 m deviation from Boussinesq is not the mesh.** The discretisation band there is
  ±0.003%; the deviation is 0.067%, more than twenty times larger. Refining further will not
  remove it. The obvious candidate is the finite 40 × 20 m domain standing in for a half-plane —
  a modelling choice, not a mesh property — and that is now a measured statement rather than a
  hypothesis.
- **Nested refinement is what made this possible, and its absence was measurable.** The same
  study run by asking the mesher for three smaller element sizes gave observed orders of 0.24,
  3.01 and 6.98 — orders no element here can deliver — because three separate visits to the
  mesher produce three unrelated meshes whose differences carry the mesher's irregularity as
  well as the discretisation. With nesting the same probes give 3.65 and 3.90, both consistent
  with a quadratic element whose displacement converges at O(h³) and whose recovered stress is
  superconvergent.
- **Being coarse is not the same as being in the asymptotic range.** A first nested attempt from
  a 2.8 m base (299 → 1 196 → 4 784 elements) was still oscillatory: the coarsest level simply
  did not resolve a 4 m strip. The family had to start fine enough for the trend to exist before
  the extrapolation could see it. The 4 m probe is still oscillatory even in the final family,
  and its band is therefore the conservative one — an honest ±1.066% rather than a confident
  ±0.067%.

## 5. The number that matters most: a factor of safety that follows the mesh

Strength reduction is the everyday computation of practical geotechnics, and it is the one where
a mesh-convergence study changes what may be claimed. `tests/test_slope_fos_convergence.cpp`
(**KV-SLP-003**) runs the Griffiths & Lane (1999) Example 1 benchmark on nested families.

**First, what the reference actually says.** Their Table 2, read from the authors' own copy:

| Trial FOS | E′δ_max/γH² | Iterations |
|---|---|---|
| 0.80 | 0.379 | 2 |
| 1.00 | 0.381 | 10 |
| 1.20 | 0.422 | 20 |
| 1.30 | 0.453 | 41 |
| 1.35 | 0.544 | 792 |
| **1.40** | 1.476 | **1000 — did not converge** |

Their factor of safety is the first trial value at which their viscoplastic algorithm fails to
converge within an iteration ceiling of 1000: "the non-convergence option is taken as being a
suitable indicator of failure". So the published 1.4 is a **bracket, (1.35, 1.40]**, quantised
by the trial increments they chose and dependent on their ceiling. Bishop & Morgenstern's charts
give 1.380. Any comparison quoting "±x% against Griffiths & Lane" without that sentence is
claiming a resolution the reference does not have.

**What KATAI measures.** Two nested families of the same case, 2026-08-08:

| Family | Elements | h [m] | FoS | Character |
|---|---|---|---|---|
| base 2.2 m (gated) | 142 / 568 / 2 272 | 1.245 / 0.622 / 0.311 | 1.4210 / 1.3867 / **1.3645** | differences shrink → reads as monotonic *convergence*, band ±3.7% |
| base 1.6 m | 267 / 1 068 / 4 272 | 0.908 / 0.454 / 0.227 | 1.4217 / 1.3867 / **1.3093** | differences grow → monotonic **divergence**, band ±8.6%, −7.9% over a fourfold refinement |

At the density the case file itself asks for, the factor of safety is 1.3867 — inside the
published bracket and +0.49% from the charts. But it does not stay there: **the finer the mesh,
the lower the factor of safety**, and the coarse family's apparent convergence is an artefact of
not having refined far enough.

**This is the method, not a defect, and the literature says so.** With ψ = 0 the flow rule is
non-associated, failure localises into a shear band, and without a regularisation the band width
is set by the elements — so refining narrows the band and lowers the computed factor. In the
words of the group that works closest to PLAXIS: *"The result obtained from a phi/c reduction is
influenced by the mesh size, element type and convergence tolerances"* (N. Torggler, *Numerical
Studies of Embedded Beam Row in Safety Analysis*, TU Graz, Institute of Soil Mechanics and
Foundation Engineering, §3.1.1, citing F. Tschuchnigg, H. F. Schweiger and S. W. Sloan (2015), "Slope stability
analysis by means of finite element limit analysis and finite element strength reduction
techniques. Part I: Numerical studies considering non-associated plasticity", *Computers and
Geotechnics* 70:169–177). The same source reports that for zero dilatancy the factor of safety
"cannot be uniquely determined" because the failure mechanism itself varies during the reduction.

**What this program does about it.** Since v0.7 a Safety phase whose material has ψ < φ raises
`K2D-A005` on every run: the factor of safety is stated as mesh-dependent, at the mesh it was
computed on, with a recommendation to confirm it by refinement. A number that moves with the
mesh must never be handed over as if it did not — that is the same rule as WP-1's, applied to a
result instead of an input.

**And below the default stopping rule it is not the solver.** The same case at the same mesh,
solved at tolerated force residuals of 1e-3 (what the strength-reduction search uses when nothing
is asked for), 1e-4, 1e-6 and 1e-8, returns **1.421021 every time — bit for bit**
(`tests/test_tolerance_independence.cpp`, **KV-NUM-007**). Five orders of magnitude of stopping
rule move the factor of safety by 0.0000%, while a fourfold mesh refinement moves it by 4–8%. The
comparison against Griffiths & Lane and against the charts is therefore a statement about the
model and the mesh, not about when a Newton loop was allowed to stop — which is exactly the
confusion a reader is entitled to suspect, given that the reference's own factor of safety *is*
defined by a convergence criterion.

**Above it, the solver is exactly the answer, and it lies the unsafe way.** The strength-reduction
search asks one question of each trial — did this reduced strength reach equilibrium? — and reads
"the solver stopped" as "yes". A loose stopping rule therefore does not add scatter, it adds
**bias**, always upward:

| Tolerated residual | Factor of safety | vs the default rule |
|---|---|---|
| 3e-1 | 2.999683 | +111% — the search's own cap: *nothing ever failed* |
| 1e-1 | 2.069116 | **+45.6%** |
| 1e-2 | 1.449585 | +2.0% |
| **1e-3** (the search's own) | **1.421021** | — |
| 1e-4 / 1e-6 / 1e-8 | 1.421021 | 0.0000%, bit-identical |

A slope reported 45% safer than it is would be precisely the number this project exists to
refuse, so since 2026-08-09 a Safety phase asked for anything looser than 1e-3 raises
**`K2D-A006`** with these figures in the message, and the validator refuses a tolerated error
of 1 or more outright.

**How this was found is part of the record.** Until the same day, the strength-reduction search
hard-coded its trial tolerance, so the sweep above set a number nothing read. It returned three
identical factors of safety, and those were reported here as independence. A sweep over a
silently ignored input is indistinguishable, in its output, from a genuine insensitivity — which
is why "the numbers did not move" is evidence only once the input is known to arrive. The
controls were threaded into the search (and, in `.k2d` v7, into the file) and the study re-run;
the conclusion above is the re-measured one, and it carries the loose-side edge the original
could never have found.

What quantises our factor instead is the strength-reduction **search**: the safety strategy
bisects the reduction factor between 0.4 and 3.0 for 12 iterations, so the reported number
carries a resolution of 2.6/2¹² = 6.3 × 10⁻⁴ — about 0.05% at 1.38, finer than the mesh
dependence by two orders of magnitude, and finer than the 0.05 trial increments the reference
stepped through.

**What remains.** A regularisation (a non-local or gradient formulation, or a fixed band width),
an adaptive refinement driven by an error estimator, and a comparison against finite element
limit analysis — which brackets the collapse load from above and below and does not localise —
are the routes the literature takes, and none of them is implemented here.

## 6. The hardening family's 1% residual, measured

The driver gives the Hardening Soil and Soft Soil families a tolerated force residual of **1e-2**
— a hundred times looser than the 1e-6 it gives Mohr-Coulomb — because their tangent is a
continuum rather than a consistent one, so the Newton iteration converges linearly instead of
quadratically. That is a defensible engineering choice, and it was an undeclared one: nothing
measured what it costs.

`tests/test_input_corpus.cpp` gained the first HS **boundary-value** case for exactly this
purpose (**KV-CST-002**): a laterally confined, weightless HS column whose vertical stress steps
from 50 to 200 kPa. The vertical stress is then the surcharge itself, uniform top to bottom, and
the settlement has a closed form — the integral of the HS oedometric stiffness law,
`E_oed = E_oed^ref ((c cos φ + σ₁ sin φ)/(c cos φ + p_ref sin φ))^m`, which with c = 0 reduces to
`−ε₁ = (p_ref^m/E_oed^ref)·[σ₁^(1−m)]/(1−m)`. Until now the HS family was verified only at the
material point; this verifies the path from the file through the mesher, the cap return mapping
and the load stepping.

| Tolerated residual | Settlement of the loading step | vs the closed form |
|---|---|---|
| **1e-2** (the driver's default) | 0.018934 m | **+0.413%** |
| 1e-4 | 0.019045 m | +1.002% |
| 1e-6 | 0.019044 m | +0.998% |

Three things to read off it. The default costs **0.59%** on this quantity — bounded, and no
longer unknown. **1e-4 is already converged** *at this step count*: two more orders change the
answer by 0.005%. And the default's *smaller* deviation from the closed form is **cancellation,
not accuracy**: the tolerance-converged answer is +1.0% at these 40 increments, and the looser
run happens to sit nearer the analytic value on the way there. Reporting it as "more accurate"
would be exactly the kind of luck this document exists to strip out.

⚠️ **That +1.0% was left unaccounted for here, and §7 now accounts for it.** It is not one error
but three, and only one of them is a discretisation: the mesh contributes nothing at all on this
case, the load path is worth about 0.46% of it at 40 increments, and the ~0.54% remainder is a
model deviation rather than a numerical one. The row above is therefore a tolerance sweep at a
fixed step count, not a statement that the case is converged.

## 7. The nonlinear family: where the error actually lives, and why this case gets no band

Everything above bands a MESH. For a path-dependent model that is the wrong axis to start on,
and KV-CST-002 is where it shows. Its converged deviation from the closed form was recorded in
§6 as "+1.0%" with no account of what the 1% was made of. **KV-NUM-009** takes the same case
apart on three axes independently, and the answer is that only one of them is a discretisation
at all.

| axis swept | range swept | what it is worth |
|---|---|---|
| **mesh** | 0.5 → 0.125 m, **85 → 1105 nodes** | **5e-15 relative** — round-off |
| **iteration tolerance** | 1e-2 → 1e-8 | 1e-6 and 1e-8 agree to six figures; the HS default of **1e-2 is not converged** and makes a step sweep non-monotone |
| **load increments** | 10 → 160, tolerance converged | **+3.47% → +0.54%**, i.e. 2.9 percentage points |
| what is left | — | **≈ +0.54%**, and it is not numerical |

**The mesh contributes nothing, and that is a fact about the case rather than a limitation of
the sweep.** The column is weightless, so σ₁ is the surcharge and the strain field is uniform; a
uniform field lies exactly in the element space, so refinement has nothing to improve. This is
the estimator's `Exact` branch meeting a real problem. A mesh band published for this case would
have been a fiction dressed as rigour.

**The load path dominates — and it is not a Richardson parameter.** The estimator needs
`φ(d) = φ_exact + C·d^p`. A path-dependent integration with its own adaptive substepping does not
supply one, and the sweep says so in its own numbers: the observed order computed from the three
overlapping triplets of a **single monotone-looking sweep** comes out **2.204, 0.469 and 1.132** —
a factor of 4.7 apart. Triplets that disagree with themselves by that much are not in an
asymptotic range, and 1/(r^p − 1) turns a small wobble at p = 0.47 into a large and
confident-looking correction. Refined further the sweep stops improving at all: past about 160
increments it settles onto a **±0.02% noise floor** and starts to oscillate. So this case reports
its measured **spread** and refuses a GCI. That refusal is the result, not a gap in it.

**What survives is a model deviation, and its signature identifies it.** After an exact mesh, a
converged tolerance and a 16× refined path, ~0.5% remains. A pure offset in the fitted stiffness
would show the same relative deviation on every stress range; this one does not:

| stress range | deviation from the closed form |
|---|---|
| 50 → 100 kPa | **−0.2425%** |
| 100 → 200 kPa | **+0.6419%** |
| 200 → 400 kPa | **+1.0952%** |

Near zero at the reference pressure and growing away from it **with a sign change** — the mark of
a cap whose α and β were calibrated at p_ref (`hs_calibrate_cap`, as PLAXIS derives them by
simulating an oedometer). It is a difference between the model and the closed form, not between
the computation and the model. Reporting it as numerical uncertainty would be wrong in both
directions: it would inflate the number and blame the wrong thing.

**And there is a ceiling on path refinement, declared because a study that cannot be repeated is
not a study.** The tolerated error is an ABSOLUTE force residual, so shrinking the increment does
not shrink what each increment has to achieve. Refining the SEATING phase of this weightless
column — whose confining stress starts near zero, where the HS stiffness is smallest — to 160
increments at 1e-6 does not converge at all rather than converging better. The sweep above
therefore pins the seating phase at 40 and moves only the staged phase, and KV-NUM-009 pins the
failure itself, so that the day it stops failing the sentence gets rewritten.

## 8. The one case that tests the estimator back: KV-STR-003's peak moment

Every sweep above runs where the discretisation error is unknown — which is the point of an
estimator, and also the reason an estimator is hard to trust. The beam-bending case is the
exception, and it was worth banding for that reason alone.

Its distributed-load beam has a peak moment whose error is known **in closed form**: the moment
field is a parabola and the element's curvature is linear, so the peak overshoots by exactly
`q h²/12` — reproduced to five figures at three densities. The order of that error is therefore
not something to be discovered from the data; it is analytically **2**, decided before any run.
Feeding those three moments to the same `grid_convergence_band()` that bands everything else:

| h [m] | peak moment [kNm] | overshoot vs 50.0 |
|---|---|---|
| 0.5 | 52.08333 | +2.08333 |
| 0.25 (the file's own mesh) | 50.52083 | +0.52083 |
| 0.125 | 50.13021 | +0.13021 |

- observed order **p = 2.0000**, monotonic convergence, inside the asymptotic range;
- Richardson value **50.00000 kNm** — the manual's published number, recovered from three meshes
  **none of which produces it**;
- reported band on the file's own mesh: **± 0.3247 %** (GCI at the observed order, Fs = 1.25).

The order it recovers is the order the algebra already knew, and the value it extrapolates to is
the value the reference publishes. That is a check of the estimator, not of the beam — and it
costs nothing, because the three runs were already there as a witness that the overshoot was a
bias rather than scatter. It also puts a number on the warning attached to that witness: a peak
moment read off a coarse run is not merely "a bit high", it is high by an amount this procedure
will quantify for any mesh triplet the user cares to produce.

## 9. What this does not yet cover

The register is deliberately explicit about its own gaps, since an absent row must never read as
a passed one:

- **Most of the corpus has no sweep yet.** The corpus is now **26 files** behind **55 declared
  cases**; four carry a band (KV-FND-008 via KV-NUM-005, KV-SLP-002, the Giroud rigid footing,
  and KV-STR-003 as of 2026-08-14), and one has been swept and **deliberately given none**
  (KV-CST-002 via KV-NUM-009, §7). The arithmetic of the gap is not the obstacle it looks like:
  a three-level nested sweep on a linear elastic case costs about **2.6 seconds**, so what is
  missing here is written oracles and probes, not machine time. The nonlinear families cost more
  — KV-NUM-009's three sweeps take about 170 s together — but §7 shows that the expense is not
  the real difficulty either: **the difficulty is that a mesh triplet is often not the axis the
  error is on**, and finding the axis has to precede banding it.
- **No band exists for a path-dependent quantity, and §7 explains why rather than promising one.**
  Richardson extrapolation requires a single smooth discretisation parameter. Load-step
  refinement is not one here, and the honest output for such a case is a measured spread plus the
  reason. What would change that is a constitutive integration whose local error is controlled to
  a set tolerance rather than by an adaptive substep count with a ceiling, and a RELATIVE
  convergence criterion so that refining the path does not eventually make every increment
  unsatisfiable (the ceiling measured in §7). Both are engine changes, not reporting changes.
- **Tolerance independence is measured for two families, not all.** KV-NUM-007 covers the slope
  factor of safety (strength-reduction trials at 1e-3, spread 0.0000% below it and unsafe-sided
  above) and the Hardening Soil oedometer KV-CST-002, whose tolerance §7 now carries out to 1e-8:
  1e-6 and 1e-8 agree to six figures, so that family's converged tolerance is known and the
  shipped default of 1e-2 is measured to be short of it. The soft-soil family, the consolidation
  and the dynamic paths have not been swept. The controls are now per phase in
  the `.k2d` (v7: `tol`, `loadsteps`, `maxiter`) and reach the solver by the same route as the
  driver seam, which **KV-NUM-008** pins as an identity — a control that is read, validated and
  echoed but not applied would leave a sweep proving nothing, which is not a hypothetical.
- **One open question in the record is now measured, one is still open.** The slope factor of
  safety falling on fine meshes is §5 above: characterised, cited and declared to the user. The
  −39% deep sheet-pile wall row still has only its domain-size diagnosis, not its closure.
- **The 4 m probe is still not in the asymptotic range**, and no amount of refinement policy
  will change that: it needs either a finer family or an understanding of why that depth behaves
  differently from 2 m. Its band is the conservative one until then.
- **Only uniform refinement exists.** `refine_uniform` refines everywhere; a graded or adaptive
  refinement — finer where the error is, which is what the deep-wall case needs — is not
  implemented, and neither is an error estimator to drive it.
