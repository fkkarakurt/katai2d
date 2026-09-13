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

Mesh refinement *controls* (global and local coarseness factors, refinement around structures and
load lines) and qualitative guidance — finer is more accurate, balance accuracy against run time —
do not put a number on the discretisation error of a result, and neither do convergence criteria.
It is worth being precise about what the latter govern: they decide when an **iterative** (Newton)
loop may stop, which is a different question from how far the converged answer of a given mesh
sits from the exact solution of the same equations. KATAI's convergence criteria are a separate
matter (§11).

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
words of one study of the method: *"The result obtained from a phi/c reduction is
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
| **1e-2** (the driver's default) | 0.018680 m | **−0.932%** |
| 1e-4 | 0.018639 m | −1.152% |
| 1e-6 | 0.018647 m | −1.109% |

Three things to read off it. The default costs **0.179%** on this quantity — bounded, and no
longer unknown. **1e-4 is converged** *at this step count*: two more orders change the answer by
0.043%. And the default's *smaller* deviation from the closed form is **cancellation, not
accuracy**: the tolerance-converged answer is −1.1% at these 40 increments, and the looser run
happens to sit nearer the analytic value on the way there. Reporting it as "more accurate" would
be exactly the kind of luck this document exists to strip out.

⚠️ **These numbers were re-measured on 2026-08-20 and they moved** (the row read +0.413% / +1.002%
/ +0.998% before). The reason is in §7: the stress-point integrator was holding the
stress-dependent moduli fixed at the state each *increment* began from, which is a first-order
error in the increment size. It is fixed, and the whole of §7 is rewritten around what is left.
One consequence belongs here rather than there: below about 1e-4 it is **no longer the stopping
rule that sets the floor** — the constitutive integration tolerance does, which is why the two
tight runs agree to 0.043% instead of to six figures.

⚠️ **That −1.1% is not one error but two**, and §7 separates them: the mesh and the load path
contribute essentially nothing on this case, and what remains is the model's own distance from the
idealised power law — measured there without any finite element method at all.

## 7. The nonlinear family: where the error actually lives, and why this case gets no band

Everything above bands a MESH. For a path-dependent model that is the wrong axis to start on,
and KV-CST-002 is where it shows. **KV-NUM-009** takes the same case apart on the axes that are
actually there, and the answer is that none of the discretisations is what is left over.

> ### ⚠️ What this section used to say, and why it was wrong
>
> Until 2026-08-20 this section reported the load path as worth **2.9 percentage points** over a
> 16× refinement, and concluded that the residual after it — about 0.5% — was a **model** deviation
> identified by its *signature*: it grew with the stress level and changed sign near `p_ref`,
> which is where the cap is calibrated.
>
> The signature was real and the conclusion drawn from it was not. Both were produced by a defect
> in the stress-point integrator: `E_ur`, `E_i`, `q_a` and `q_f` are functions of the current minor
> principal stress, and they were being evaluated **once, at the state each increment started
> from**, then held for every substep of that increment. That is a first-order error in the load
> step, and no substep tolerance can see it — an error-controlled integrator reported "converged"
> while the answer still moved with the increment count. Measured at the stress point, on one
> oedometer path at a fixed integration tolerance of 1e-6, the final σ₁ over a 32× refinement of
> the outer step ran 773.27 → 791.47 → 800.46 → 805.02 → 807.27 → 808.48 kPa: every halving halved
> what was left, extrapolating to ≈809.7. **The "converged" answer was 4.5% low at 20 steps and
> still 0.6% low at 160.** Reading the moduli at each substep instead puts the whole sweep inside
> **6e-5 relative**, on the number the old sequence extrapolates to.
>
> An argument from a shape is only as good as the shape. This section no longer makes one.

| axis swept | range swept | what it is worth |
|---|---|---|
| **mesh** | 0.5 → 0.125 m, **85 → 1105 nodes** | **5.6e-16 relative** — round-off |
| **iteration tolerance** | 1e-2 → 1e-6 | 0.222% in total; 1e-4 and 1e-6 agree to **0.043%** |
| **load increments** | 10 → 160, tolerance converged | **0.020%** (it was 2.9 percentage points before the integrator was fixed) |
| **constitutive integration tolerance** | STOL 1e-3 → 1e-6 | −1.536% → −0.982%; the default is **1e-5** |
| what is left | — | **≈ −1.0%**, and §7a shows it belongs to the model |

**The mesh contributes nothing, and that is a fact about the case rather than a limitation of
the sweep.** The column is weightless, so σ₁ is the surcharge and the strain field is uniform; a
uniform field lies exactly in the element space, so refinement has nothing to improve. This is
the estimator's `Exact` branch meeting a real problem. A mesh band published for this case would
have been a fiction dressed as rigour.

**The load path still gets no GCI, on a better reason than before.** The estimator needs
`φ(d) = φ_exact + C·d^p`, and a sequence has to converge before an order can be read off it. This
one no longer does: over 10 → 160 increments the five settlements move by 0.020% in total, and of
the three overlapping triplets one reverses sign, so the axis now sits **inside its own noise**. There is nothing
to extrapolate, and quoting a Richardson correction built on noise would invent a number. The
measured spread is reported instead — and it is now roughly a fiftieth of the deviation that
remains, which is what makes that deviation attributable to something else.

**A fourth axis is now declared, because it exists.** The constitutive integration carries its own
tolerance (`STOL`), and it is a numerical control like the tolerated error and the step count, not
a property of the soil. The default is **1e-5**, chosen by measurement rather than taste: a looser
integration is not cheaper, because once the integration noise rises above the residual the phase
is asking for, the equilibrium iteration grinds against it.

| STOL | deviation | Newton iterations (seating / staged) | wall clock |
|---|---|---|---|
| 1e-3 | −1.536% | 1887 / 4883 | 91 s |
| 1e-4 | −1.310% | 1907 / 1427 | 76 s |
| **1e-5** | **−0.986%** | **432 / 1186** | **31 s** |
| 1e-6 | −0.982% | 387 / 3379 | 117 s |

### 7a. What is left, measured without the finite element method

The remaining ≈1% is claimed as a model deviation, and the claim is made by taking the finite
element method out of the measurement rather than by reading a shape.

The same material — built through the same registry entry the driver builds it with, so its
calibrated cap parameters α and β are the run's own — is walked as a one-dimensional oedometer at
a single **stress point**: no mesh, no load path, no equilibrium iteration, no linear solver. If
the boundary-value run and the constitutive routine land in the same place, then whatever separates
both of them from the closed form is the model's distance from the idealised power law.

| stress range | boundary-value run | stress point, no FE | apart by |
|---|---|---|---|
| 50 → 100 kPa | **−1.3328%** | −1.3490% | 0.016 pp |
| 100 → 200 kPa | **−0.9861%** | −0.9518% | 0.034 pp |
| 200 → 400 kPa | **−0.6081%** | −0.9130% | 0.305 pp |

So **at most about a third of a percentage point** of the deviation can be the finite element
method, and the rest is the model. Integrated more tightly (STOL 1e-6) the stress point reads
−1.398% / −1.003% / −0.781% / −0.605% over 50→100, 100→200, 200→400 and 400→800 kPa: **one-signed,
largest just below `p_ref`, and closing as the stress rises**. The model is stiffer than
`E_oed^ref (σ/p_ref)^m` everywhere, by about one and a half per cent near the reference pressure
and by half a per cent an octave above it. That is a different statement from the one this section
used to make, it is smaller than the old sweep's own uncertainty, and it does not need a finite
element run to reproduce.

### 7b′. The backends now agree

The case's 100 → 200 kPa run was the one that split between the two linear-solver backends —
**+0.641936%** on PARDISO against **+0.458298%** on Eigen, a 0.18 pp disagreement that this record
published as evidence. Re-measured with the integration under a declared tolerance:

| stress range | PARDISO (MKL) | Eigen (portable) |
|---|---|---|
| 50 → 100 kPa | 7.706387545820947e-03 m | 7.706387545820932e-03 m |
| 100 → 200 kPa | 1.093676800855408e-02 m | 1.093676800855408e-02 m |
| 200 → 400 kPa | 1.552598037192273e-02 m | 1.552598037192282e-02 m |

Fifteen significant figures on two of them and bit-for-bit on the third. A disagreement that
survived four orders of tolerance was not in the linear algebra and was not in the equilibrium
iteration; it was in what each iterate's stress path did to a subdivision that had no tolerance of
its own.

**And there is still a ceiling on path refinement, declared because a study that cannot be repeated
is not a study.** The tolerated error is an ABSOLUTE force residual, so shrinking the increment does
not shrink what each increment has to achieve. Refining the SEATING phase of this weightless
column — whose confining stress starts near zero, where the HS stiffness is smallest — to 160
increments at 1e-6 still does not converge rather than converging better. The sweep above therefore
pins the seating phase at 40 and moves only the staged phase, and KV-NUM-009 pins the failure
itself, so that the day it stops failing the sentence gets rewritten. (It did stop failing, briefly,
at an integration tolerance of 1e-4, and started again at the shipped 1e-5 — which is itself a
measurement of what the ceiling is made of.)

### 7b. The other half of the same argument: KV-NUM-010, where the axis *is* clean

A refusal is only worth something if the same procedure produces a band when it should. Otherwise
"no band" and "the machinery does not work" are the same observation. **KV-NUM-010** is the
control: Terzaghi consolidation, swept by exactly the estimator that bands meshes, on the axis its
error actually lives on.

**The order was decided before the run, as it was for the beam.** `consolidation.hpp` integrates
**fully implicitly (α = 1)** — backward Euler, global error O(Δt). The prediction is therefore
**p = 1**, and it is not something to be discovered from the data. (That this works twice does not
make it a rule; §7c is the case where the same reasoning predicts 2 and the measurement returns 3.) The default `OrderPolicy`
window of [0.5, 4.0] is not used either: that bracket exists because it is what a tri6 *mesh* can
deliver, and a time axis is entitled to its own, so the sweep is run with [0.75, 1.5] and a
fallback order of 1 rather than 2.

| time steps (Tv = 0.5) | U | vs Terzaghi |
|---|---|---|
| 30 | 0.758035991 | −0.7742% |
| 60 | 0.760967928 | −0.3904% |
| 120 (the file's own) | 0.762449266 | −0.1965% |
| 240 | 0.763193873 | −0.0990% |

- observed order **p = 0.9850** from 120/60/30 and **0.9924** from 240/120/60 — the two agree to
  **0.007**. *That* is what an asymptotic range looks like, and it is precisely what §7's load
  path could not produce (2.204 / 0.469 / 1.132 from one sweep);
- Richardson **U(Δt→0) = 0.763961996** against the series' 0.763950331 — **+0.0015%**, the closed
  form recovered from three time steps alone;
- band on the file's own 120 steps: **± 0.2480%**, and the true error there is **0.1965%**, so the
  band **contains** it — checkable here rather than merely trusted, because the exact answer is
  known.

And the mesh, for the second corpus case running: **1.7e-8** over an 11× node count. Two cases in
a row whose error lives where no mesh sweep would ever have found it. The pair is the point:
§7 refuses a band because the axis will not carry one, §7b publishes one because it will, and the
same estimator and the same reasoning produced both.

### 7c. Where "the algebra fixes the order" stops: KV-NUM-011

Twice now the order has been known before the run — 2 for the beam's `q h²/12`, 1 for backward
Euler. That is a tempting rule, and **KV-NUM-011 is the case that bounds it.** `dynamics.hpp`
integrates with **γ = ½, β = ¼** — average acceleration, no numerical damping, unambiguously
**second order** — so the prediction was p = 2. The measurement says **3**.

**Two axes, and the first one is not a discretisation at all.** A resonant amplitude is *built
up*, not imposed: the envelope approaches steady state like `1 − exp(−ξωt)`, so the file's 20
cycles is 0.9981 of the way there. Sweeping Δt without knowing that would charge the shortfall to
the integrator.

| duration | \|u_surf\| | vs closed form |
|---|---|---|
| 10 cycles | 0.049374659 | −4.3172% |
| **20 cycles (shipped)** | 0.051494210 | **−0.2098%** |
| 40 cycles | 0.051588308 | −0.0274% |

At the shipped settings the buildup is worth about as much as the time step. With the duration
taken out to 40 cycles, the time axis alone:

| time steps (40 cycles) | \|u_surf\| | vs closed form |
|---|---|---|
| 800 | 0.049972794 | −3.1581% |
| 1600 | 0.051401970 | −0.3885% |
| 3200 | 0.051563299 | −0.0759% |

- **observed order p = 3.147**, where the scheme alone justifies 2;
- 3.147 and 2.998 at 160 cycles, and **3.097 / 2.795 on a grid deliberately INCOMMENSURATE with
  the period** (T/Δt = 18.43, 36.86, 73.72, 147.45) — which **eliminates sampling phase-lock** as
  the explanation. The order is about 3, reproducibly, on two unrelated grids;
- the axis is nonetheless clean: Richardson lands on the closed form to **−0.036%**, and the two
  independent grids agree with each other to 1e-5;
- band on this triplet: **± 0.0498%**.

**What the quantity is explains it.** The published number is the *peak* of a resonant response,
not the response at an instant, and a derived quantity does not have to inherit its scheme's
order. So the rule is: **the order is a property of the quantity as much as of the integrator, and
it may only be asserted in advance when the quantity's relation to the discretisation is itself
known** — as it is for the beam's algebraic overshoot and for Terzaghi's U, and as it is not here.
The test asserts `p > 2.3` rather than `p ≈ 3.147`, because the finding is that the prediction
fails, and pinning the exact value would go quiet the day the mechanism changed.

And the mesh, for the **third** corpus case running: 459 → 6389 nodes moves \|u\| by 4.5e-6.

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
- Richardson value **50.00000 kNm** — the closed-form peak moment q l²/8, recovered from three meshes
  **none of which produces it**;
- reported band on the file's own mesh: **± 0.3247 %** (GCI at the observed order, Fs = 1.25).

The order it recovers is the order the algebra already knew, and the value it extrapolates to is
the closed-form value. That is a check of the estimator, not of the beam — and it
costs nothing, because the three runs were already there as a witness that the overshoot was a
bias rather than scatter. It also puts a number on the warning attached to that witness: a peak
moment read off a coarse run is not merely "a bit high", it is high by an amount this procedure
will quantify for any mesh triplet the user cares to produce.

## 9. What this does not yet cover

The register is deliberately explicit about its own gaps, since an absent row must never read as
a passed one:

- **Most of the corpus has no sweep yet.** The corpus is now **26 files** behind **58 declared
  cases**; six carry a band (KV-FND-008 via KV-NUM-005, KV-SLP-002, the Giroud rigid footing,
  KV-STR-003 as of 2026-08-14, KV-CON-002 via KV-NUM-010 and KV-DYN-002 via KV-NUM-011, the last
  two on their TIME axes), and one has been swept and **deliberately given none** (KV-CST-002 via
  KV-NUM-009, §7). The arithmetic of the gap is not the obstacle it looks like:
  a three-level nested sweep on a linear elastic case costs about **2.6 seconds**, so what is
  missing here is written oracles and probes, not machine time. The nonlinear families cost more
  — KV-NUM-009's three sweeps take about 170 s together — but §7 shows that the expense is not
  the real difficulty either: **the difficulty is that a mesh triplet is often not the axis the
  error is on**, and finding the axis has to precede banding it.
- **The consolidation and dynamic paths are now swept; the soft-soil one is not.** KV-NUM-010
  bands the time axis of KV-CON-002 and KV-NUM-011 that of KV-DYN-002, both validated against a
  known answer. Soft Soil Creep has a time axis of the same kind and has not been touched. Nor
  has the *coupled* (flow-deformation) path, which has both.
- **Why the observed order of KV-NUM-011 is 3 rather than the scheme's 2 is measured, not
  explained.** Phase-lock between the sampling grid and the period was tested and eliminated; the
  remaining account — that a peak of a resonant response is a derived quantity with an order of
  its own — is consistent with everything measured but is not derived. The band does not depend
  on the account: it is computed from the observed order, and it is validated by landing on the
  closed form. But an explanation would let the order be predicted for the next such quantity
  instead of discovered.
- **No band exists for a path-dependent quantity, and §7 explains why rather than promising one.**
  Richardson extrapolation requires a single smooth discretisation parameter. Load-step
  refinement is not one here, and the honest output for such a case is a measured spread plus the
  reason. What would change that is a constitutive integration whose local error is controlled to
  a set tolerance rather than by an adaptive substep count with a ceiling, and a RELATIVE
  convergence criterion so that refining the path does not eventually make every increment
  unsatisfiable (the ceiling measured in §7). Both are engine changes, not reporting changes.
- **Tolerance independence is measured for two families, not all.** KV-NUM-007 covers the slope
  factor of safety (strength-reduction trials at 1e-3, spread 0.0000% below it and unsafe-sided
  above) and the Hardening Soil oedometer KV-CST-002, whose tolerance §6 carries out to 1e-6:
  1e-4 and 1e-6 agree to 0.043%, so that family's converged tolerance is known and the shipped
  default of 1e-2 is measured to cost 0.179% on this quantity. Below about 1e-4 it is the
  CONSTITUTIVE integration tolerance rather than the equilibrium one that limits the answer (§7),
  which is why those two runs agree to a twentieth of a per cent instead of to six figures. The soft-soil family, the consolidation
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
## 10. The uncertainty that was not discretisation, and was not the linear algebra either (KV-STR-004)

Everything above bands the distance between a computed number and the exact solution of its own
equations, with the mesh or the time step as the axis. This section is about a different axis, and
it was found the way the useful ones usually are: a test that had been green locally for a month
was red on every CI run.

The case is KV-STR-004, the axial capacity of a pile row. Its mesh-independence probe halved the
element size under the fixture's 1500 kN/m head load. On the MKL composition the refined run
converges and the pile carries **600.0000 kN/m**, its declared capacity, exactly as at the coarse
density. On the `portable` composition — the same source, the same file, the same mesh, the
vendored Eigen solver instead of PARDISO — the run reported a collapse mechanism after
equilibrating **10%** of the load.

Both could not be right. What follows is the measurement that decided it, because the first two
explanations offered here were both wrong.

### What it is not

**It is not the linear algebra.** The obvious reading — PARDISO performs iterative refinement on
its solutions and the Eigen backend does not, so the Eigen tangent solve must be failing near a
singular stiffness — was written into this section, and two guarded refinement steps were added to
`eigen_direct_solver.cpp` to test it. They did not recover the run. Instrumenting the Newton loop
afterwards showed why: through the whole failing run **the linear solver refused nothing**. Every
increment that was abandoned had been abandoned with a perfectly good factorization behind it. Two
refusals do appear, but only at the very end, after the run had already spent its cut-backs.

**It is not the collapse detection.** This section previously suspected "a discontinuity in the
cut-back and mechanism-detection logic near the collapse load". There is no discontinuity, and the
load is not near a collapse: the same file, on the same mesh, converges to full load on both
backends once one setting is changed.

### What it is

The fixture's soil is Mohr-Coulomb with c′ = 1000 kPa — far above anything the run mobilises — and
it is **weightless**, so every stress in the model comes from the head load. It also carries the
tension cut-off, at σ_t = 0. A weightless body loaded at one point has tensile principal stresses
over large regions, so the cut-off is active nearly everywhere the load reaches, and the material
there is a *no-tension* one: zero stiffness in the capped direction. That is what the run is
actually solving, and it is why it crawls. Measured on the coarse mesh, same file, one field
changed:

| tension cut-off | Newton iterations, phase 2 | per increment | wall clock |
|---|---|---|---|
| on (as shipped) | **561** | 28–50 | 5.94 s |
| off | **57** | **2** | 0.61 s |

Two iterations per increment is what a consistent tangent is supposed to give. Twenty-eight to
fifty is a scheme fighting a degenerate stiffness: the Newton direction is periodically enormous —
`|δ|` measured at ten to fifty times `|du|` — and the backtracking line search has to cut α to
10⁻³ or 10⁻⁴ to find any descent at all, which is an iteration spent for almost no progress.

On the refined mesh that crawl runs into the per-increment iteration limit, which was **80**:

| composition | limit | outcome | max&#124;u&#124; | iterations | cut-backs | wall clock |
|---|---|---|---|---|---|---|
| MKL / PARDISO | 80 | converges | 0.1822857 m | — | 4 | 66.4 s |
| MKL / PARDISO | 200 | converges | **0.1822167 m** | 1278 | 0 | **37.2 s** |
| portable / Eigen | 80 | **"collapse at 10%"** | — | — | 6 | not timed (the run fails) |
| portable / Eigen | 200 | converges | **0.1822167 m** | 1369 | 0 | 220 s |

The increments on this mesh need up to **96** iterations, and their out-of-balance force falls
monotonically the whole way — nothing about them is a mechanism. PARDISO's path needed 88 and
survived a budget of 80 by cutting back four times; Eigen's path needed 96 and did not. That is the
entire disagreement: **a patience setting decided it, and which side of the setting a run lands on
turns on the linear solver's rounding.** With the budget raised, the two backends agree to seven
significant figures on a quantity neither of them could produce before.

Note the third column of the table as well. Raising the budget made the MKL run **1.8× faster**,
because a truncated increment does not merely fail — it triggers a cut-back, and re-solving a
halved increment from scratch costs more than the iterations the increment was denied.

### What was changed

The driver now derives the per-increment iteration limit from the material class, exactly as it
already derived the load-step count and the tolerated error: **200** where the tangent is
nonsymmetric — nonlinear soil, interfaces, embedded beams — and 80 otherwise, which is a limit no
linear run approaches. A file's own `maxiter`, and the caller's numerical controls, still win over
it.

Raising it is free where a mechanism genuinely forms. A collapsing run does not exhaust its
iteration budget; it ends when the tangent goes singular along the mechanism and the solver refuses
to return a direction. Measured on the same fixture with the soil weakened to c′ = 5 kPa and given
weight: **the same 41% limit load and the same 122 iterations at a limit of 80 and at 200.**

And the phase now says which of the three things happened, because they are not interchangeable and
one sentence used to be printed for all of them:

- **no descent** at any step size, or **the linear solver refusing** on a singular tangent — a
  mechanism. The equilibrated fraction is the incremental limit load, and the message says so.
- **the iteration budget running out** while the residual was still falling — a setting. The
  message names the limit that stopped the run and explicitly does *not* publish a limit load.

`test_non_convergence_names_its_reason` (KV-STR-005) pins both directions on one fixture: past
capacity it must publish the limit load; at half capacity with the limit set to 1 it must refuse to,
because the same model carries that load when the budget is adequate.

### The crawl itself: measured, understood, and deliberately not traded away

Raising the budget stopped a setting from deciding the answer, but it did not explain why the
increments wanted 50 to 96 iterations for a problem whose soil is elastic apart from a cap. That was
measured next, and the measurement named it exactly:

| line-search step accepted | measured |
|---|---|
| full step, `α = 1` | 172 of 541 iterations |
| halved once or twice | 349 of 541 |
| worst increment | 50 iterations |
| the same fixture at 8.4 residual assemblies per iteration | KV-FND-010, where `α` is typically 0.0005 |

The full Newton step is rejected in about two thirds of the iterations and the search halves it, so
the run advances at roughly half a Newton step at a time. The step is not wrong; the **test** is. A
strictly decreasing Armijo condition demands that the residual norm fall at every iteration, and on
a non-smooth problem it does not: each step moves a few stress points across the yield surface — a
tension cut-off is the extreme case, since the material has no stiffness across the cap and gets all
of it back the moment the point unloads — so the norm can rise for one iteration while the iterate
is on its way to equilibrium.

Two remedies were implemented and measured on this tree. Neither is in the program, and the second
is the more interesting refusal.

**Levenberg-Marquardt damping — refuted by measurement.** The obvious reading of a degenerate tangent
is to damp it: solve (K + μ·diag|K|)δ = r, raise μ when a step is rejected, relax it when one is
taken. It made every case worse, and the pile fixture stopped converging at all. The trace says why:
μ climbed to its ceiling within a few iterations, the step shrank to 1.2·10⁻⁶ against a displacement
of 3.6·10⁻³, and the iterate froze. Damping shortens a step that is too long; it cannot help a step
that is the right length and is being judged by the wrong yardstick.

**A non-monotone line search — it works, and it is not shipped.** Judging a trial step against the
worse of the last two residuals rather than only the last (Grippo, Lampariello and Lucidi 1986) does
exactly what the diagnosis predicts: the full step is then taken in **319 of 482 iterations**, the
worst increment closes in **37** instead of 50, the fixture set runs **8–27% faster**, and four of
the five cases print the same numbers to every figure they print.

What it costs is measurable in precisely the place this record exists to measure. The Hardening Soil
family's default stopping rule is 1%, so its answers depend on the path the iterates take, and a
looser line search takes a different one:

| | as shipped | with the window |
|---|---|---|
| KV-NUM-007, oedometer, 1e-4 against 1e-6 | 0.05% | 0.05% at a window of two, **0.36%** at five |
| KV-NUM-007, oedometer, 1e-2 against 1e-6 | 0.72% | 0.72% at two, **1.12%** at five (the case publishes better than 1%) |
| KV-NUM-009, observed orders of the load-step sweep | 2.204 / 0.469 / 1.132 | 2.751 / **−0.600** / 1.655 |
| KV-NUM-009, residual over 50-100, 100-200, 200-400 kPa | −0.2425% / **+0.6419%** / **+1.0952%** | −0.2425% / **+1.1839%** / **+1.0665%** |

(The 16× spread of that sweep is 2.915% either way; it is the shape of the sequence that moves, not
its extent.)

The last row is the one that settles it. That case identifies the ~0.5% left after an exact mesh and
a converged tolerance as a **model** deviation rather than a numerical one, and the evidence is its
signature: near zero at the reference pressure and growing away from it. With the window the growth
is no longer monotone, so the signature the conclusion rests on is no longer what the run shows.

None of these are wrong answers. They are the same equilibria reached along a different path, each
inside the tolerance its phase declares. But they are the numbers this record publishes, and moving
them is a decision about the record — not a side effect to be taken along with a speed-up. So the
crawl stays: it is a cost the program pays on non-smooth problems, its cause is understood, and the
price of removing it is a re-measurement of the Hardening Soil numerics register with its
conclusions re-derived from the new data.

⚠️ **Re-read this section against 2026-08-20.** Every number in the table above was measured with
the stress-point integrator that §7 has since retracted — the one holding the stress-dependent
moduli fixed across an increment. Two consequences, and neither of them reverses the decision by
itself:

* **The last row's argument is gone.** It turned on the window destroying a *signature* (growth with
  stress level, sign change at `p_ref`) that the conclusion rested on. That signature was the
  defect, not the model, so a line search can no longer be charged with destroying it.
* **The rest of the table is not evidence any more either.** The load-step sweep it compares
  against is now flat to 0.020% and its "observed orders" are noise; the tolerance rows were
  measured where the equilibrium tolerance still set the floor, and it no longer does below 1e-4.

The line search therefore returns to the queue as an **open question with no current measurement**,
not as a decided trade. What it costs has to be measured again against the integrator that ships
now, on a record whose conclusions no longer depend on the shape of that particular sweep.

### What this still leaves open

Two smaller things were measured on the way and are recorded here rather than fixed silently: a
line search that finds no descent in twelve halvings still applies its last, smallest step, which
was measured turning a residual of 2.57 into 534 on one iteration; and the relative residual in the
linear solver's refusal message is formatted with `std::to_string`, which prints six decimals, so a
refusal at 3·10⁻⁶ reads "0.000003" and one at 3·10⁻⁹ reads "0.000000".

**What this costs the project's claims, stated plainly.** The README and the project website say
that the `portable` preset reproduces every published number with the vendored Eigen solver alone.
That was true of every number this record publishes, and it was *not* true of the 1500 kN/m
refined-mesh run — which was never a published number, but was asserted by a test, which is the
same promise in a different place. It is true of that run now, on the default settings, on both
compositions.

---

## 11. The criterion that was never checked, and what it changes (KV-NUM-013)

Every band in the sections above rests on the same sentence: *the run converged*. Until now that
sentence meant one thing in this tree — a global force residual `||r|| <= rtol·max(||f_ext||, ||f_c||, 1)`
— and it was reported as though it meant the whole question. It does not. A complete convergence
check is a **family**: a CSP-normalised global force error, a moment residual on rotational
freedoms, and, separately, the local error at every soil stress point. Checking one member alone
leaves the independent questions the others ask unasked.

This section is the first measurement of the other members on this tree.

### 11.1 The two stresses of a stress point

The concept under the local criteria is that a stress point carries **two**
stresses during an iteration:

* the **constitutive** stress, what the material law returns for the strain the point was given,
  σ_c,j = σ_0 + D_e (Δε_j − Δε_p,j);
* the **equilibrium** stress, what the finite-element linearisation says the point carries,
  σ_eq,j = σ_c,j−1 + D_e δε_j, built from the previous iterate and this iteration's correction.

Their difference, normalised by max(τ_max, c', 1 kPa) for a plastic point (the plastic-point
criterion) or by max(τ_max, c', p_ref/200) for a point whose *elastic* stiffness depends on stress
(the non-linear elastic criterion), is
the local error. It is an error and not a measure of nonlinearity: at the exact solution δε = 0
and the two stresses coincide identically. The linear-elastic control measures it at **0.0**, and
the first iterate of every increment measures it at 2e-16.

### 11.2 Two orientations of the stiffness parameter, and which way it was resolved

CSP — the Current Stiffness Parameter — normalises the global force error, and it can be written
**two different ways up**:

| | written as | value while elastic | value at failure |
|---|---|---|---|
| energy ratio, elastic over total | elastic energy increment / total energy increment | 1 | **grows without bound** |
| work ratio, actual over elastic | integral of Δε·Δσ over integral of Δε D_e Δε | 1 | **towards 0** |

Only the work-ratio orientation is consistent with the behaviour the parameter exists to have —
equal to unity when the solution is fully elastic, approaching zero at failure — and with the uses
built on it: arc-length control engages below CSP 0.5, and collapse is reported below CSP 0.015.
Implementing the energy-ratio form would make the global force criterion **loosen** as a mechanism
forms, which is precisely where a load fraction is about to be published as a bearing capacity.

The work ratio is what is implemented, deliberately, and the measurement agrees with it. On one strip
footing, load rising:

| footing load [kPa] | CSP | load factor |
|---|---|---|
| linear elastic control | **1.0** (to round-off) | 1.000 |
| 100 | 0.43743 | 1.000 |
| 300 | 0.09765 | 1.000 |
| 600 | 0.00012 | 0.656 |
| 900 | 0.00009 | 0.434 |

The parameter falls monotonically to below the collapse value of 0.015 on the two loads the
model cannot carry. `test_convergence_criteria` pins that direction, so the inverted reading cannot
be reintroduced quietly.

### 11.3 The separation: the global criterion improves by 1e5 while the local one improves by 3

Mohr-Coulomb strip footing, q = 300 kPa, tolerated error 1e-6, iterate by iterate:

| iterate | global force error | worst local error | inaccurate plastic points |
|---|---|---|---|
| 1 | 1.30e-2 | 6.54e-2 | 335 / 335 |
| 2 | 1.29e-2 | 7.09e+0 | 344 / 344 |
| 3 | 3.39e-3 | 7.83e+0 | 347 / 347 |
| 4 | 8.63e-5 | 5.60e-1 | 347 / 347 |
| **5 (accepted)** | **3.71e-7** | **2.33e-2** | **337 / 347** |

The run is declared converged with its global error a factor of 2.7 *below* what was asked, while
337 of its 347 plastic points carry a local error four decades above it. The separation is not a
constant: it is 5x at the first iterate and **63 000x** at the last. Whatever the local criterion
is measuring, the force residual cannot see it.

Tightening the tolerance shows the local error is genuinely converging, not stalled — and where its
floor is:

| tolerated error | iterations | global force error | worst local error | settlement max abs u [m] |
|---|---|---|---|---|
| 1e-2 | 74 | 3.63e-3 | 8.10e+0 | 0.314188623 |
| 1e-4 | 114 | 8.63e-5 | 5.60e-1 | 0.314678453 |
| 1e-6 | 128 | 3.71e-7 | 2.33e-2 | 0.314684142 |
| 1e-8 | 141 | 7.56e-12 | 7.69e-5 | 0.314684146 |
| 1e-10 | 149 | 7.56e-12 | 7.69e-5 | 0.314684146 |

The local error follows the global one down until the solve reaches its round-off floor, where it
settles at 7.69e-5 — 0.008% of the local shear strength. **The displacement, on this case, is
already converged to eight figures at 1e-6.** So on a Mohr-Coulomb footing the local criteria are
strictly informative: they report a stress-point statement that the displacement answer does not
depend on.

### 11.4 Where it is not merely informative: the Hardening Soil oedometer at its shipped tolerance

KV-CST-002, run three ways. 1e-2 is what this tree ships for the Hardening Soil family. 1e-6 is
four decades tighter and is what this record uses
when it wants the answer rather than a run.

| run | iterations | wall clock | settlement max abs u [m] | vs the tight answer |
|---|---|---|---|---|
| tolerance 1e-2, global criterion only (**what shipped until 2026-08-24**) | 114 | 1.3 s | 0.018680445 | **+0.179%** |
| tolerance 1e-2, local criteria binding (**what ships now**) | 194 | 3.3 s | 0.018647102 | +0.0005% |
| tolerance 1e-6, global criterion only | 399 | 7.6 s | 0.018647012 | — |

The shipped stopping rule stops **0.18% away** from its own converged answer, and the local
criteria are what detect it. Binding them at the *same* tolerance lands on the converged settlement
in **194 iterations instead of 399** — that is, the local criteria are a cheaper route to the right
answer than tightening the global tolerance is. This was not the expected result: the prediction
written down before the run was that enforcement would be expensive and might prevent convergence
near a limit load. It did neither. On the two footing loads the model cannot carry, enforcement
changed the outcome not at all (load factor 0.656 and 0.434 either way, 464 to 478 and 506 to 519
iterations).

### 11.5 The structural half: interfaces and the pile toe

Two more members complete the family, and both are the same construction applied to a traction
instead of a stress tensor. A slipping interface point carries an equilibrium shear
tau_eq,j = tau_c,j-1 + k_s * delta-u_s,j and a constitutive shear tau_c,j from the Coulomb return;
the interface criterion normalises their difference by the point's own capacity. An embedded beam's
SKIN coupling springs are counted in the SAME tally — the criterion draws no distinction between a
soil-structure interface and the special interface a pile skin is. The foot is not a point at all:
the foot criterion forms ONE ratio over every foot in the model, tolerated at FIVE times the
tolerated error.

Measured on the two corpus cases that have the elements, with only the global criterion binding:

| case | global force error | what the structural criteria say |
|---|---|---|
| KV-STR-002 sliding block (1e-6) | 2.03e-7, met | **14 of 47 slipping interface points inaccurate**, worst 2.05e-6 |
| KV-STR-004 axial pile (1e-6) | 1.29e-9, met | skin + interface 0 of 32; **foot force error 5.6e-9 of 5.0e-6 allowed** |

The sliding block is the clean interface case: its soil is elastic, so nothing in the soil
contributes a plastic point and the interface criterion is the only local one with anything to
say. It fails while the global criterion is met by a factor of five. The pile's foot criterion is
satisfied comfortably, which is worth stating as plainly as a failure: a criterion that only ever
reports trouble is not measuring anything.

### 11.6 The decision: the local criteria now bind, and what that moved

Everything above was measured before anything was decided. The decision is taken: as of
2026-08-24 `enforce_local_criteria` is **ON** by default, in the engine, the phase strategy and
the driver. `KATAI_CONV_NOLOCAL` turns it off for a whole run, which is how every comparison in
this section is reproduced.

It was taken on the measurement rather than on the principle. The principle would have said
"check the whole family"; the measurement says something stronger — that this is not a
stricter rule bought with iterations, but a **cheaper route to the same answer**. On KV-CST-002 it
lands on the converged settlement in 194 iterations where tightening the global tolerance by four
decades costs 399.

**What moved, across the whole suite: three cases out of 154.** Every other published number is
unchanged, which is itself the useful measurement — the criteria bind where the old rule was
stopping early and nowhere else.

**1. KV-NUM-007, the slope factor of safety: a whole class of unsafe-sided answer is gone.**

| tolerated residual | factor of safety, before | after |
|---|---|---|
| 1e-1 | 2.069116 (**+45.6%**) | 1.429578 (**+0.6%**) |
| 1e-2 | 1.449585 (+2.0%) | 1.421021 (+0.0%) |
| 1e-3 (the search's own) | 1.421021 | 1.421021 |

This was never really about the tolerance. A strength-reduction search asks "did this reduced
strength still reach equilibrium?", and a run that stops on the global force residual alone can
stop with its stress points nowhere near the strengths it has just reduced them to — which the
search reads as a yes. That gap is exactly what the local criteria close. A slope reported 45%
safer than it is was the kind of number this record exists to refuse, and the stopping rule can no
longer produce it.

**2. KV-CST-010, Soft Soil Creep: a published 3% band that was measured on an under-converged
run.** Requiring the criteria moved every duration of this case AWAY from the idealised creep law,
which looked at first like the criteria being wrong for a viscous point. Sweeping the tolerance
with them OFF settled it — the model walks to its own answer and stays there:

| duration | idealised law | shipped tolerance, global only | swept to 1e-8, global only | with the local criteria |
|---|---|---|---|---|
| 1 d | 2.772589e-3 | 2.715531e-3 (−2.06%) | **2.936937e-3 (+5.93%)** | 2.941283e-3 (+6.08%) |
| 10 d | 9.591581e-3 | 9.747126e-3 (+1.62%) | **9.892608e-3 (+3.14%)** | 9.922439e-3 (+3.45%) |
| 100 d | 1.8460e-2 | 1.8598e-2 (+0.74%) | **1.8763e-2 (+1.64%)** | 1.8794e-2 (+1.80%) |

The old run was not accurate; it was **lucky**. It stopped short, and where it stopped happened to
sit inside 3% of a law the model does not exactly obey. The converged deviation is +5.93% falling
to +1.64%, and the local criteria find it at the shipped tolerance in 81 iterations against the
sweep's 189. The case's band is now 7% across the sweep, and what carries the physics is no longer
a flat band but the SHAPE: the deviation must FALL with duration, because the idealised law drops
the elastic and consolidation parts and those matter most where there is least creep. A flat band
cannot distinguish a model that obeys the law from one that happens to land near it once.

**3. KV-STR-002, the sliding block: a plateau that is no longer bit-identical.** Doubling the
imposed slip used to give the same failure force to the last bit; the two runs now stop one
iterate apart and agree to 1.4e-9 relative — nine significant figures. The band moved from 1e-9 to
1e-6, which still leaves six decades between "this is a plateau" and "this is a stiffness
reading", the two hypotheses it exists to separate.

**Two guard sentences decayed and were rewritten, both in KV-NUM-008.** One proved that the file's
tolerated error is READ by showing that the answer moved when it changed; it moves by 0.0009% now,
so the proof is taken where the control lands instead — the run reports the tolerance it ran under,
it demonstrably MET it, and reaching it cost iterations. The other pinned the iteration limit this
case needs at 3; the threshold is now FOUND by the test rather than written down, because that is
the second time a number pinned there moved because the solver got better.

### 11.7 The two criteria that had never been read, and what they turned out to measure

**First, a correction.** The paragraph this section used to open with said the moment residual
"participates in the gate". It did not. The gate called `local_ok()`, and `local_ok()` is the
three local counts plus the pile toe — the moment residual was computed, stored, printed, and
never consulted. That is the third guard sentence in this document to decay the same way, and the
pattern is now specific enough to name: **a sentence describing what the code checks decays as
soon as the check moves, because nothing fails when the sentence stops being true.**

Both halves are closed. It is in the gate now (`Convergence::enforced_ok`, 2026-08-25), and it has
been read through the only two kinds of case in this tree that carry rotational freedom.

**Case 1 — an ELASTIC plate: the criterion cannot bind, by construction.** The coupled wall +
excavation fixture (`test_excavation_wall`: Mohr-Coulomb soil under K0, a 12 m bonded sheet-pile
wall, staged excavation to 9 m), swept over the tolerance it is solved at:

| tolerated | wall tip u_x [m] | iterations | moment residual | force residual (CSP-normalised) |
|---|---|---|---|---|
| 1e-1 | −5.087516e-3 | 111 | 8.51e-15 | 3.13e-2 |
| 1e-2 | −6.459047e-3 | 206 | 1.54e-14 | 7.03e-7 |
| 1e-3 | −6.457048e-3 | 328 | 1.25e-14 | 1.95e-11 |
| 1e-6 | −6.457049e-3 | 407 | 1.66e-14 | 4.44e-14 |
| 1e-10 | −6.457049e-3 | 468 | 1.25e-14 | 2.74e-14 |

The moment residual does not move at all — nine decades of tolerance, one decade of scatter around
1e-14 — and the reason is structural rather than numerical. An elastic plate's internal force is
`K_p u`, so the rotational rows are LINEAR in the plate degrees of freedom: whatever the soil is
doing, the linear solve satisfies them exactly and the next assembly reads back round-off. The
criterion is therefore incapable of failing on an elastic plate, **including the 1e-1 run, whose
tip deflection is 21.2% away from the converged one.** It reports "ok" next to an answer that is a
fifth wrong, and that is not a defect: it is answering a question about the rotational rows, and
the run is caught by the global force gate instead. What would be a defect is a record that lets
the "ok" be read as covering the answer, which is what this table exists to prevent.

**Case 2 — a PLASTIC plate: the criterion is alive.** The M-N hinge BVP (`test_plate_plastic`: a
plate on an elastic bed, central load past its plastic moment) makes the rotational rows nonlinear,
and then the residual behaves like every other residual — it follows the tolerance down:

| tolerated | moment residual | force residual | ratio |
|---|---|---|---|
| 1e-2 | 1.530e-4 | 2.433e-3 | 0.063 |
| 1e-3 | 4.522e-5 | 8.298e-4 | 0.054 |
| 1e-4 | 8.355e-13 | 1.014e-10 | 0.008 |

Over 15 (load, tolerance) pairs — 1.1× to 5× the plate's limit load, tolerances 1e-1 to 1e-3 — the
ratio never exceeded **0.063**. On everything this tree runs today the rotational rows converge at
least as fast as the translational ones, so **binding the moment criterion costs nothing**: no case in the suite
changed by a single iteration when it was switched on.

It is bound anyway, for the case the tree does not have yet. The global gate is `||r||` — a
Euclidean norm over a vector whose entries are forces AND moments — and adding kN to kNm in
quadrature is not a norm of anything. Where the wall is what fails and the soil around it is still
elastic, that sum is dominated by a soil which is not being asked for much, and the moment criterion is the only
dimensionally honest question left about the rows that matter.

**What binding it actually cost: exactly one case, and that case found a missing floor.** The
claim above — that binding the moment criterion is free — was made on the two fixtures and then tested against
the whole suite, which is the only place such a claim can be tested. One test failed:
`test_no_silent_drop`'s plate standing along the whole of a line that is pushed down. Every node
of that plate takes the same settlement, so it translates without curving, and a plate that does
not curve carries no moment at all:

> **Correction (2026-09-13).** This paragraph first said the plate was *undriven* and that the run
> said so itself (`K2D-A003`). The numbers below were right and the explanation was not: the
> structural elements had read a driven node's motion since 2026-08-13, and `K2D-A003` was a stale
> warning, now retired. The plate moves; it does not bend.

| the same plate | reference sum \|M_nodal\| | moment error | outcome |
|---|---|---|---|
| standing on a pushed line (translated, not bent) | 1.11e-12 kNm/m | 2.06e-1 | **phase refused at load factor 0** |
| under a 100 kPa strip load | 4.71e+02 kNm/m | 4.5e-14 | converged in 2 iterations |

Fourteen orders of magnitude of denominator. The ratio that refused the phase was one round-off
divided by another. And the defect it exposed is not the binding: **every other criterion in this
family has a floor and the moment criterion did not**. The plastic-point criterion divides by
`max(tau_max, c, 1 kPa)`, the foot criterion by `max(|F_c|, 1% of |F_max|, 1 kN)`, and the reason
for them is that a point carrying almost nothing must not report an enormous relative error on a
difference that is numerically nothing. A structure carrying almost no moment is that situation in
different units. The floor is 1 kNm/m — ten decades above the translated plate's reference and two
below the loaded one — and it is recorded here as what it is: not part of the moment criterion as
first formulated, adopted from the device its siblings use, and forced by a measurement rather than
chosen by taste.

`has_moment` changed with it. It used to mean "there is a rotational freedom **and** its reference
is non-zero", which quietly excused the one case that needed saying; it now means what it says, and
an unloaded structure reports a small number instead of vanishing from the report.

**The same shape of answer on a different criterion: the non-linear elastic one.** This document has said three times
that the non-linear elastic criterion "has never been the binding count", and blamed the case list:
the Hardening Soil oedometer is plastic everywhere the moment it loads, and the Mohr-Coulomb cases
have no stress-dependent modulus at all. Run on the case that is nothing but such points — an
HSsmall **unloading**, where every point leaves its yield surfaces and re-enters on `Eur`:

| tolerated | non-linear elastic points | inaccurate | worst non-linear elastic error |
|---|---|---|---|
| 1e-2 | 96 | 0 | 1.29e-15 |
| 1e-4 | 96 | 0 | 1.29e-15 |
| 1e-6 | 96 | 0 | 1.29e-15 |

The points are counted — the criterion is wired, and this is the first case that proves it — and
the error is not small, it is **zero, and it does not move**. The cause is not case selection: it is
the integration scheme. `hs_frozen_Eur` evaluates the unloading modulus at the **committed** state
and holds it for the whole increment, which is what makes the tangent consistent and the line
search safe. A point that does not yield therefore walks the increment with a **constant** elastic
operator, and the equilibrium stress is built as `sigma_c,j-1 + D^e delta-eps` — exactly what
the constitutive routine then returns. The two stresses of §11.1 coincide identically. The
non-linear elastic criterion measures their difference, so in this tree it measures zero **by construction**, and it is a code
that updates the modulus inside the iteration where it would have something to say.

That makes it a live check of an identity rather than a dead criterion, and the test says so: if
the modulus ever starts iterating, or the elastic operator handed to the probe stops being the one
the return actually used, this count is the first thing that moves.

**A side finding that outlives this subsection: CSP is blind to structural plasticity.** Across
every plastic-plate run above — hinge fully formed, `|M|` saturated at `Mp`, loads to 5× the limit
load — the Current Stiffness Parameter reads **1.00000**, i.e. "fully elastic". The energy
integrals it is built from are accumulated in the soil Gauss loop and nowhere else, so a mechanism
that forms in a structure is invisible to it. That is not an issue for the parameter's use in this
document so far (every case it has been read on plastifies through the soil), but it bears directly
on two things that are planned on top of it: the force error's denominator is scaled by CSP, and the
arc-length switch is specified to engage below CSP 0.5. Both would sit still through a wall
collapsing.

### 11.8 What this section still does not cover

* **The global force criterion still uses the old normalisation as its gate — and the swap has
  now been measured, so the decision is a small one.** Both ratios are on the record at the same
  iterate (`Convergence::force_error` for the CSP-normalised ratio and `Convergence::global_error` for the gate —
  until 2026-08-25 only the one that does NOT gate was published, so a run could print "force error
  3e-2, tolerance 1e-1" while the quantity that let it stop was a different ratio entirely; the
  command line and the Studio report went on printing that one as "the force error" until
  2026-09-13, when both began to lead with the gate and `.res` v10 began to store it), and
  `KATAI_CONV_CSPGATE` swaps the gate over. What the swap is worth, measured on the Mohr-Coulomb
  footing walked to collapse:

  | case | CSP | CSP-normalised | the gate | CSP-normalised is stricter by |
  |---|---|---|---|---|
  | q = 100 kPa (converges) | 0.437 | 1.87e-15 | 3.31e-15 | 0.6× (looser) |
  | q = 300 kPa (converges) | 0.098 | 7.56e-12 | 8.46e-12 | 0.9× (looser) |
  | q = 600 kPa (collapses) | 0.00012 | 3.46e-3 | 2.41e-3 | **1.4×** |
  | q = 900 kPa (collapses) | 0.00009 | 2.36e-3 | 1.09e-3 | **2.2×** |

  The CSP normalisation does tighten as a mechanism forms, which is the property it is there for,
  by a factor of 1.4 to 2.2 at collapse; on every run that converges the two agree within about
  1.5× and both sit twelve or more decades below any tolerance the tree ships. Read alone, that
  table says the swap is a wash. **It is not, and the suite says so: with `KATAI_CONV_CSPGATE` set,
  5 of the 142 fast-lane tests fail** — and the clearest of them, `test_anchor_prestress`,
  equilibrates *0% of its load*.

  The reason is the same defect this subsection's other half found in the moment criterion, in the
  same place: **the denominator has no floor tied to what is being applied.** The CSP-normalised
  ratio divides by
  `||f_int|| + CSP·||f_const||`, and its only floor is `1e-6` of the phase's load scale — six
  decades down, so it bites only when literally nothing is resisting. A phase that prestresses an
  anchor against ground which has not responded yet has very little internal force, so the
  denominator is small, the ratio never reaches the tolerance, and every increment is abandoned.
  The gate in place divides by `max(||f_ext||, ||f_const||, 1)`, which contains the load being
  applied and therefore cannot collapse that way.

  So the open question is not "which normalisation do we prefer". It is that **the CSP-normalised
  ratio's denominator needs the floor discipline its siblings have before it can gate anything**,
  and until that floor is derived rather than invented, the swap is a measurement seam and not
  a candidate. The staged excavation this was predicted to separate on cannot show any of it: the
  only staged case in the corpus stays elastic (CSP = 1), where the two denominators differ by a
  factor of 1.45 at 1e-16.
* **The `portable` composition has not been re-measured** with the criteria binding. Every number
  in this section is from the MKL build.
* **KV-CST-012 does not converge at a tolerance tighter than the 1e-2 it ships with**, and the
  phase that refuses has been named. Sweeping the whole case fails from 1e-3 down (load factor
  0.041, 0.034, 0.011); tightening the phases ONE AT A TIME shows every other phase — the initial
  one, the unload, the excavation — solving at 1e-3 unchanged, and **`Surcharge on` alone
  reproducing the failure exactly**. It is not the local criteria: with them off
  (`KATAI_CONV_NOLOCAL`) the same phase fails at the same place. So the published answer for that
  case cannot presently be shown to be converged — the same class of statement as KV-CST-010's
  band, which turned out to be luck — and the phase that has to be understood is the first one.
* **Nothing in the criteria family is measured in a DYNAMIC phase.** The nonlinear Newmark solver
  has its own Newton loop and does not build the probe at all, so a shaking phase's stress points
  are held to the global force residual alone — the state this section was written to end.

