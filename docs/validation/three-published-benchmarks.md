# Three published benchmarks, end to end

This document walks three well-known published problems through KATAI 2D
exactly the way a user would run them: every material parameter, geometry
dimension and reference value is taken letter for letter from the source,
the input is a checked-in `.k2d` file, and the run is the public `katai`
command line. Nothing here is a one-off study — each case is part of the
verification corpus, re-solved from its file by `test_input_corpus` on
every build, and each appears as a row of the generated
[verification matrix](verification-matrix.md).

The three cases cover three different kinds of answer a geotechnical FE
engine must get right: a **collapse load in axisymmetry** (Cox), a
**collapse load with depth-varying strength** (Davis & Booker), and a
**factor of safety by strength reduction** (Griffiths & Lane).

| Case | Published value | KATAI 2D | Deviation |
|---|---|---|---|
| Cox (1962) — smooth rigid circular footing, limit pressure | 225.6 kPa (slip-line, exact) | 233.9 kPa | +3.7% |
| Davis & Booker (1973) — strip footing on c(z) clay, limit pressure | 7.80 kPa (analytic) | 7.91 kPa | +1.4% |
| Griffiths & Lane (1999) Example 1 — slope factor of safety | 1.380 (Bishop–Morgenstern charts) / 1.4 (their FE) | 1.384 | +0.3% / −1.1% |

Both bearing cases are also published finite-element benchmarks: the
PLAXIS 2D Validation Manual (Version 8, Bentley Systems) runs the same two
problems and publishes its own results next to the analytic values — 220.0
kPa (−2.5%) for Cox and 7.86 kPa (+0.8%) for Davis & Booker. Those numbers
are quoted below alongside KATAI's so the reader can see where an
independent, established FE implementation lands on the same problem.
Displacement-type finite elements approach a slip-line limit load from
either side depending on discretisation and formulation choices; deviations
of a few percent against an exact collapse load are the expected class of
result for both programs.

---

## 1. Bearing capacity of a smooth rigid circular footing — Cox (1962)

**Source.** Cox (1962) slip-line solution for the indentation of a
ponderable c–φ soil, as published in Section 3.1 of the PLAXIS 2D
Validation Manual, Version 8: limit pressure `p_max = 141 c = 225.6 kPa`.

**The problem, verbatim.** Axisymmetric Mohr–Coulomb soil cylinder, 5 m
radius × 4 m deep; footing radius `R = 1 m`, smooth and rigid.
`E = 2400 kPa`, `ν = 0.20`, `c = 1.6 kPa`, `φ = 30°`, `γ = 16 kN/m³`,
`K0 = 0.5`. The flow rule is **associated** (`ψ = φ = 30°`): a slip-line
solution is the associated limit load, so a non-associated run would mix a
modelling difference into a verification number. No tension cut-off — the
slip-line solution has none.

**The input.** [`tests/corpus/kv-fnd-013-cox-circular-footing.k2d`](../../tests/corpus/kv-fnd-013-cox-circular-footing.k2d)
— a K0 geostatic initial phase, then a staged phase that pushes the footing
line `(0, 4)–(1, 4)` down by a prescribed `u_y = −0.35 m` with `u_x` left
free (smooth contact). The soft soil needs that large an indentation to
reach the collapse plateau. Mesh: 0.25 m 15-node triangles, from the file's
own settings.

**The run.**

```text
> katai validate tests/corpus/kv-fnd-013-cox-circular-footing.k2d
OK: tests/corpus/kv-fnd-013-cox-circular-footing.k2d satisfies the input contract (0 warning(s))

> katai solve tests/corpus/kv-fnd-013-cox-circular-footing.k2d
phase 1/2: Initial phase
phase 2/2: Indent
phase 1/2: ok  max|u| = 0.000000e+00 m
phase 2/2: ok  max|u| = 4.845403e+00 m
solved 2 phase(s) in 172.51 s
```

The initial K0 phase displacing by exactly zero is itself an assertion: the
geostatic seed *is* the equilibrium. The large `max|u|` of the indentation
phase is the free-surface flow of the collapse mechanism at a 0.35 m punch
into soft soil — the imposed footing settlement itself is carried exactly.

**Reading the limit pressure.** The `.k2d` contract carries the reaction
output: at the fixed footing DOFs the results hold the support reactions
(the discrete internal force `Bᵀσ`, no recovery smoothing). Axisymmetric
nodal forces are per radian, so the footing force is `2π|ΣR_y|` and the
average pressure is `p = 2|ΣR_y| / R²`. From Python (the same contract the
CLI writes):

```python
import katai

pr, notes = katai.load_project("tests/corpus/kv-fnd-013-cox-circular-footing.k2d")
job = katai.run(pr)
r = job.results()[-1]
Ry = sum(r.reaction[2*n + 1] for n in range(len(r.node_x))
         if r.node_y[n] > 4 - 1e-6 and r.node_x[n] <= 1 + 1e-6)
print(f"p_max = {2*abs(Ry):.1f} kPa")   # 233.9
```

**Result.**

| | p_max [kPa] | vs analytic |
|---|---|---|
| Cox (1962), slip-line (exact) | 225.6 | — |
| PLAXIS 2D (published, same manual) | 220.0 | −2.5% |
| **KATAI 2D** (file's own 0.25 m tri15 mesh) | **233.9** | **+3.7%** |

**What the suite additionally asserts.** The K0 initial phase displaces by
exactly `0.000e+00` (the geostatic seed *is* the equilibrium); a footing
node carries exactly the imposed `u_y = −0.35 m`; and the mesh-sensitivity
lesson is recorded — a 0.5 m mesh puts only two elements across the radius
and over-predicts the collapse by about +9%, the classic coarse-mesh
bearing-capacity bias.

---

## 2. Strip footing on clay with strength increasing with depth — Davis & Booker (1973)

**Source.** Davis & Booker (1973) analytic solution for a smooth strip
footing on clay whose undrained strength grows linearly with depth, as
published in Section 3.2 of the PLAXIS 2D Validation Manual, Version 8:
`p_max = ρ [(2 + π) c₀ + B c_inc / 4] = 7.80 kPa`.

**The problem, verbatim.** Plane strain, Tresca (`φ = 0`), weightless.
`c(z) = 1 + 2z` kPa (`c₀ = 1 kPa` at the surface, `c_inc = 2 kPa/m`),
`E(z) = 299 + 498z` kPa, `ν = 0.3`. Footing half-width 1 m (`B = 2 m`),
smooth and rigid, pushed 30 mm. In the `.k2d` schema the two depth profiles
are the material's `c_inc` / `E_inc` / `y_ref` fields — this case pins
exactly those against a theoretical collapse load. No tension cut-off.

**The input.** [`tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d`](../../tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d)
— a weightless gravity initial (an exact nil by construction), then the
prescribed-displacement phase. Mesh: 0.5 m 15-node triangles. Half-model;
the average pressure under the footing is `p = |ΣR_y| / (B/2)`.

**The run.**

```text
> katai validate tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d
OK: tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d satisfies the input contract (0 warning(s))

> katai solve tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d
phase 1/2: Initial phase
phase 2/2: Indent
phase 1/2: ok  max|u| = 0.000000e+00 m
phase 2/2: ok  max|u| = 3.804946e-02 m
solved 2 phase(s) in 3.70 s
```

**Result.**

| | p_max [kPa] | vs analytic |
|---|---|---|
| Davis & Booker (1973), analytic | 7.80 | — |
| PLAXIS 2D (published, same manual) | 7.86 | +0.8% |
| **KATAI 2D** (file's own 0.5 m tri15 mesh) | **7.91** | **+1.4%** (+0.6% vs the PLAXIS number) |

**What the suite additionally asserts.** The weightless initial phase
displaces by exactly zero; a footing node carries exactly the imposed
`u_y = −30 mm`; the limit pressure sits within 5% of the analytic value
*and* within 3% of the published PLAXIS number.

---

## 3. Homogeneous slope — Griffiths & Lane (1999), Example 1

**Source.** D. V. Griffiths and P. A. Lane, *Slope stability analysis by
finite elements*, Géotechnique 49(3):387–403 (1999), Example 1 — the
homogeneous 2:1 slope with no foundation layer (D = 1). The paper reports
a factor of safety of **1.4** by finite elements (the trial at FOS = 1.4
fails to converge within the 1000-iteration ceiling — non-convergence *is*
the failure criterion) and cites **1.380** from the Bishop & Morgenstern
(1960) charts for the same slope.

**The problem, verbatim.** Slope angle 26.57° (2:1), `φ' = 20°`,
`c'/γH = 0.05`, `ψ = 0` (the paper's compromise value for slope problems),
nominal `E' = 10⁵ kN/m²`, `ν' = 0.3`. Geometry from their Fig. 1: a 1.2H
crest plateau, a 2H slope run, height H; vertical rollers on the left
boundary, full fixity at the base. No tension crack modelling — plain
Mohr–Coulomb, as in the paper. Dimensionalised here as `H = 10 m`,
`γ = 20 kN/m³`, `c' = 10 kPa`, so `c'/γH = 0.05` exactly.

**The input.** [`tests/corpus/kv-slp-002-griffiths-lane-example1.k2d`](../../tests/corpus/kv-slp-002-griffiths-lane-example1.k2d)
— the file's *initial procedure* is Safety, so loading the file and solving
it *is* the phi–c reduction. The definition of the factor of safety is the
paper's own: the number by which `c'` and `tan φ'` must be divided to bring
the slope to failure. Mesh: 1.0 m 6-node triangles.

**The run.**

```text
> katai validate tests/corpus/kv-slp-002-griffiths-lane-example1.k2d
OK: tests/corpus/kv-slp-002-griffiths-lane-example1.k2d satisfies the input contract (0 warning(s))

> katai solve tests/corpus/kv-slp-002-griffiths-lane-example1.k2d
phase 1/1: Initial phase
phase 1/1: ok  max|u| = 4.013301e-01 m  FoS = 1.384
solved 1 phase(s) in 69.12 s
```

**Result.**

| | Factor of safety | |
|---|---|---|
| Bishop & Morgenstern (1960) charts, cited in the paper | 1.380 | — |
| Griffiths & Lane (1999), finite elements | 1.4 | +1.4% vs the charts |
| **KATAI 2D** (file's own 1.0 m tri6 mesh) | **1.384** | **+0.3%** vs the charts, −1.1% vs their FE |

**What the suite additionally asserts.** The failure mechanism genuinely
displaces (`max|u| = 0.40 m` at the reduced strength — a real slip surface,
not a numerical artefact), and the factor of safety stays within 4% of the
Bishop–Morgenstern value on every build.

---

## Reproducing all of this

Any KATAI 2D install (the single-file `katai.exe` or `pip install` wheel)
reproduces the runs above with no other tooling:

```text
katai validate tests/corpus/kv-fnd-013-cox-circular-footing.k2d
katai solve    tests/corpus/kv-fnd-013-cox-circular-footing.k2d
katai validate tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d
katai solve    tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d
katai validate tests/corpus/kv-slp-002-griffiths-lane-example1.k2d
katai solve    tests/corpus/kv-slp-002-griffiths-lane-example1.k2d
```

And the guarantees are structural, not editorial: `test_input_corpus`
asserts on every build that each checked-in file is byte-identical to its
programmatic build, validates with zero errors, and reproduces the
published value inside the declared band; the reference gate regenerates
the [verification matrix](verification-matrix.md) and
[bibliography](references.bib) from the declarations inside the tests, so
this record cannot drift from the suite.
