# Four examples from the PLAXIS 2D Validation Manual, input for input

This document takes four problems from the PLAXIS 2D Validation Manual,
Version 8 (Bentley Systems) — two linear-elastic settlement problems and two
plastic collapse problems — and runs them through KATAI 2D exactly as that
manual defines them: every geometry dimension, material parameter and boundary
condition letter for letter, the input a checked-in `.k2d` file, the run the
public `katai` command line. For each problem three numbers are then placed
side by side:

1. the **closed-form reference** the manual itself verifies against
   (Giroud 1972; Gibson 1967; Cox 1962; Davis & Booker 1973),
2. the **PLAXIS 2D result as published** in that manual, quoted unchanged, and
3. the **KATAI 2D result**, from the file's own mesh settings, reproduced below
   with its full command-line transcript and the exact reading recipe.

The comparison protocol is stated up front so the numbers can be read fairly.
The PLAXIS values are published values — PLAXIS was not re-run for this
document, and the two programs do not share a mesh; each discretises the
problem with its own settings. What is identical is the *problem*: the
mechanics, the parameters, and the closed-form reference both programs answer
to. PLAXIS 2D is an independent, established finite-element implementation,
and its published validation results are quoted here so the reader can see
where such an implementation lands on the same problem. Displacement-type
finite elements approach an exact limit load from either side depending on
discretisation and formulation choices; deviations of a few percent against
an exact collapse load are the expected class of result for any such program.

None of these numbers is a one-off study. All four cases are part of the
verification corpus: `test_input_corpus` re-solves each checked-in file on
every build and asserts the values below inside declared bands, and each case
is a row of the generated [verification matrix](verification-matrix.md), with
its sources in [references.bib](references.bib).

## Summary

| Manual section | Problem | Quantity | Closed form | PLAXIS 2D (published) | KATAI 2D |
|---|---|---|---|---|---|
| 2.1 | Smooth rigid strip footing on elastic soil | footing force at s = 10 mm | 15.15 kN/m | 15.24 (+0.6%) | **15.32** (+1.1%) |
| 2.2 | Strip load on incompressible Gibson soil | centreline surface settlement | 0.050 m (half-space) | 0.047 (finite layer) | **0.0454** (−3.4% vs PLAXIS) |
| 3.1 | Bearing capacity of a smooth rigid circular footing | limit pressure | 225.6 kPa (slip-line, exact) | 220.0 (−2.5%) | **233.9** (+3.7%) |
| 3.2 | Strip footing on clay with strength increasing with depth | limit pressure | 7.80 kPa | 7.86 (+0.8%) | **7.91** (+1.4%) |

Percentages in the PLAXIS and KATAI columns are deviations from the closed
form, except for Gibson, where the finite-layer geometry makes the published
PLAXIS value the like-for-like reference (see the case).

---

## 1. Smooth rigid strip footing on elastic soil — Validation Manual 2.1

**Source.** Giroud (1972) analytic solution for a smooth rigid strip footing
on an elastic layer, as published in Section 2.1 of the PLAXIS 2D Validation
Manual: `F = 2 (1 + ν) G B s / ρ` with `ρ = 0.88`, giving `F = 15.15 kN/m`
at a settlement of `s = 10 mm`. The manual publishes `15.24 kN/m` (+0.6%)
for PLAXIS 2D.

**The problem, verbatim.** Plane strain, linear elastic, weightless.
`G = 500 kPa`, `ν = 1/3` (the `.k2d` stores the equivalent
`E = 2(1+ν)G = 1333.33 kPa`), footing half-width `B = 1 m`, smooth and
rigid: a prescribed uniform `u_y = −10 mm` with `u_x` left free. Half-model,
7 m wide × 4 m deep, the left edge the symmetry axis.

**The input.** [`tests/corpus/kv-fnd-012-giroud-rigid-footing.k2d`](../../tests/corpus/kv-fnd-012-giroud-rigid-footing.k2d)
— a weightless gravity initial phase (an exact nil by construction), then the
prescribed-displacement phase. Mesh: 0.5 m 15-node triangles, from the file's
own settings.

**The run.**

```text
> katai validate tests/corpus/kv-fnd-012-giroud-rigid-footing.k2d
OK: tests/corpus/kv-fnd-012-giroud-rigid-footing.k2d satisfies the input contract (0 warning(s))

> katai solve tests/corpus/kv-fnd-012-giroud-rigid-footing.k2d
phase 1/2: Initial phase
phase 2/2: Indent
phase 1/2: ok  max|u| = 0.000000e+00 m
phase 2/2: ok  max|u| = 1.008461e-02 m
solved 2 phase(s) in 0.35 s
```

**Reading the footing force.** The results carry the support reactions at the
prescribed DOFs (the discrete internal force `Bᵀσ`, no recovery smoothing).
The half-model carries half the footing, so `F = 2 |ΣR_y|` over the footing
nodes:

```python
import katai

pr, notes = katai.load_project("tests/corpus/kv-fnd-012-giroud-rigid-footing.k2d")
r = katai.run(pr).results()[-1]
Ry = sum(r.reaction[2*n + 1] for n in range(len(r.node_x))
         if r.node_y[n] > 4 - 1e-6 and r.node_x[n] <= 1 + 1e-6)
print(f"F = {2*abs(Ry):.2f} kN/m")   # 15.32
```

**Result.**

| | F [kN/m] | vs analytic |
|---|---|---|
| Giroud (1972), analytic | 15.15 | — |
| PLAXIS 2D (published, same manual) | 15.24 | +0.6% |
| **KATAI 2D** (file's own 0.5 m tri15 mesh) | **15.32** | **+1.1%** (+0.5% vs the PLAXIS number) |

**Interpretation.** Both programs land slightly *above* the analytic force,
and 0.5% apart from each other. The analytic solution is for an idealised
layer; a truncated finite-element domain is marginally stiffer, and the
stress concentration at the edge of a rigid punch is resolved only to mesh
resolution — both effects push the computed force up, in both programs. The
suite additionally asserts that the weightless initial phase displaces by
exactly zero and that a footing node carries exactly the imposed 10 mm.

---

## 2. Strip load on incompressible Gibson soil — Validation Manual 2.2

**Source.** Gibson (1967) closed form for a strip load on an incompressible
non-homogeneous half-space whose stiffness grows linearly with depth, as
published in Section 2.2 of the PLAXIS 2D Validation Manual: the surface
settlement under the load is uniform and equals `s = q / (2 dG/dz) = 0.050 m`
for the half-space. The manual models a **4 m finite layer** and publishes
`0.047 m` for PLAXIS 2D — below the half-space value, as a finite layer must
be.

**The problem, verbatim.** Plane strain. `E(z) = 299 z` kPa measured down
from the surface (in the `.k2d`: `E ≈ 0`, `E_inc = 299 kPa/m`, `y_ref` at
the surface), `ν = 0.495` (incompressible in the limit), strip load
`q = 10 kPa` over a 1 m half-width. Half-model, 7 m wide × 4 m deep, left
edge the symmetry axis. Weightless, so the gravity initial is an exact nil.

**The input.** [`tests/corpus/kv-fnd-011-gibson-strip-load.k2d`](../../tests/corpus/kv-fnd-011-gibson-strip-load.k2d)
— mesh 0.15 m 15-node triangles, from the file's own settings.

**The run.**

```text
> katai validate tests/corpus/kv-fnd-011-gibson-strip-load.k2d
OK: tests/corpus/kv-fnd-011-gibson-strip-load.k2d satisfies the input contract (0 warning(s))

> katai solve tests/corpus/kv-fnd-011-gibson-strip-load.k2d
phase 1/2: Initial phase
phase 2/2: Strip load
phase 1/2: ok  max|u| = 0.000000e+00 m
phase 2/2: ok  max|u| = 8.838396e-02 m
solved 2 phase(s) in 3.50 s
```

(The larger `max|u|` sits at the load edge, where the surface stiffness tends
to zero and the load terminates in a local singularity; the verified quantity
is the centreline settlement, where Gibson's solution applies.)

```python
import katai

pr, notes = katai.load_project("tests/corpus/kv-fnd-011-gibson-strip-load.k2d")
r = katai.run(pr).results()[-1]
top = max(r.node_y)
s = min(r.displacement[2*n + 1] for n in range(len(r.node_x))
        if r.node_y[n] > top - 1e-6 and r.node_x[n] < 1e-6)
print(f"s = {abs(s):.4f} m")   # 0.0454
```

**Result.**

| | s [m] | |
|---|---|---|
| Gibson (1967), half-space closed form | 0.050 | upper reference |
| PLAXIS 2D (published, 4 m finite layer) | 0.047 | −6% vs half-space |
| **KATAI 2D** (same 4 m layer, 0.15 m tri15 mesh) | **0.0454** | **−3.4% vs the PLAXIS number** |

**Interpretation.** The like-for-like reference here is the published PLAXIS
number, because both programs model the same 4 m layer while the closed form
is for a half-space; both correctly sit below it, sharing the finite-layer
bias. The case is also a discretisation stress test: at `ν = 0.495` a
low-order displacement element locks volumetrically, and a stiffness that
vanishes at the surface amplifies any such error. The 15-node (quartic)
triangles carry both without special treatment, landing within 3.4% of the
published finite-layer value.

---

## 3. Bearing capacity of a smooth rigid circular footing — Validation Manual 3.1

**Source.** Cox (1962) slip-line solution for the indentation of a
ponderable c–φ soil, as published in Section 3.1 of the PLAXIS 2D Validation
Manual: limit pressure `p_max = 141 c = 225.6 kPa`. The manual publishes
`220.0 kPa` (−2.5%) for PLAXIS 2D.

**The problem, verbatim.** Axisymmetric Mohr–Coulomb soil cylinder, 5 m
radius × 4 m deep; footing radius `R = 1 m`, smooth and rigid.
`E = 2400 kPa`, `ν = 0.20`, `c = 1.6 kPa`, `φ = 30°`, `γ = 16 kN/m³`,
`K0 = 0.5`. The flow rule is **associated** (`ψ = φ = 30°`): a slip-line
solution is the associated limit load, so a non-associated run would mix a
modelling difference into a verification number. No tension cut-off — the
slip-line solution has none.

**The input.** [`tests/corpus/kv-fnd-013-cox-circular-footing.k2d`](../../tests/corpus/kv-fnd-013-cox-circular-footing.k2d)
— a K0 geostatic initial phase, then a staged phase that pushes the footing
line down by a prescribed `u_y = −0.35 m` with `u_x` left free (smooth
contact); the soft soil needs that large an indentation to reach the collapse
plateau. Mesh: 0.25 m 15-node triangles, from the file's own settings.

**The run.**

```text
> katai validate tests/corpus/kv-fnd-013-cox-circular-footing.k2d
OK: tests/corpus/kv-fnd-013-cox-circular-footing.k2d satisfies the input contract (0 warning(s))

> katai solve tests/corpus/kv-fnd-013-cox-circular-footing.k2d
phase 1/2: Initial phase
phase 2/2: Indent
phase 1/2: ok  max|u| = 0.000000e+00 m
phase 2/2: ok  max|u| = 4.845403e+00 m
solved 2 phase(s) in 34.56 s
```

**Reading the limit pressure.** Axisymmetric nodal forces are per radian, so
the footing force is `2π |ΣR_y|` and the average pressure `p = 2 |ΣR_y| / R²`:

```python
import katai

pr, notes = katai.load_project("tests/corpus/kv-fnd-013-cox-circular-footing.k2d")
r = katai.run(pr).results()[-1]
Ry = sum(r.reaction[2*n + 1] for n in range(len(r.node_x))
         if r.node_y[n] > 4 - 1e-6 and r.node_x[n] <= 1 + 1e-6)
print(f"p_max = {2*abs(Ry):.1f} kPa")   # 233.9
```

**Result.**

| | p_max [kPa] | vs slip-line |
|---|---|---|
| Cox (1962), slip-line (exact) | 225.6 | — |
| PLAXIS 2D (published, same manual) | 220.0 | −2.5% |
| **KATAI 2D** (file's own 0.25 m tri15 mesh) | **233.9** | **+3.7%** |

**Interpretation.** The two finite-element programs bracket the exact
slip-line value from opposite sides — PLAXIS 2.5% below, KATAI 3.7% above —
which is precisely how displacement-type formulations behave around an exact
collapse load: the sign and size of the bias are properties of the
discretisation, not of the mechanics. The case's own mesh study records the
sensitivity: a 0.5 m mesh puts only two elements across the footing radius
and over-predicts by about +9% (the classic coarse-mesh bearing bias); the
file's 0.25 m mesh brings that to +3.7%. The suite additionally asserts the
K0 initial phase displaces by exactly `0.000e+00` and a footing node carries
exactly the imposed `u_y = −0.35 m`.

---

## 4. Strip footing on clay with strength increasing with depth — Validation Manual 3.2

**Source.** Davis & Booker (1973) analytic solution for a smooth strip
footing on clay whose undrained strength grows linearly with depth, as
published in Section 3.2 of the PLAXIS 2D Validation Manual:
`p_max = ρ [(2 + π) c₀ + B c_inc / 4] = 7.80 kPa`. The manual publishes
`7.86 kPa` (+0.8%) for PLAXIS 2D.

**The problem, verbatim.** Plane strain, Tresca (`φ = 0`), weightless.
`c(z) = 1 + 2z` kPa (`c₀ = 1 kPa`, `c_inc = 2 kPa/m`),
`E(z) = 299 + 498z` kPa, `ν = 0.3`. Footing half-width 1 m (`B = 2 m`),
smooth and rigid, pushed 30 mm. In the `.k2d` schema the two depth profiles
are the material's `c_inc` / `E_inc` / `y_ref` fields.

**The input.** [`tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d`](../../tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d)
— a weightless gravity initial (an exact nil), then the
prescribed-displacement phase. Mesh: 0.5 m 15-node triangles. Half-model;
`p = |ΣR_y| / (B/2)`.

**The run.**

```text
> katai validate tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d
OK: tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d satisfies the input contract (0 warning(s))

> katai solve tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d
phase 1/2: Initial phase
phase 2/2: Indent
phase 1/2: ok  max|u| = 0.000000e+00 m
phase 2/2: ok  max|u| = 3.804946e-02 m
solved 2 phase(s) in 1.83 s
```

**Result.**

| | p_max [kPa] | vs analytic |
|---|---|---|
| Davis & Booker (1973), analytic | 7.80 | — |
| PLAXIS 2D (published, same manual) | 7.86 | +0.8% |
| **KATAI 2D** (file's own 0.5 m tri15 mesh) | **7.91** | **+1.4%** (+0.6% vs the PLAXIS number) |

**Interpretation.** Both programs land within 1.5% above the analytic
collapse load, 0.6% apart from each other, on a problem that couples
plasticity to a depth-varying strength *and* stiffness profile — this case
pins the schema's `c_inc` / `E_inc` / `y_ref` machinery against a theoretical
limit load.

---

## Reading the comparison as a whole

**On the elastic problems** the programs agree with each other to well under
one percent where the geometry is like-for-like (Giroud: 0.5% apart), and
KATAI sits within 3.4% of the published PLAXIS value on the finite Gibson
layer, both codes correctly below the half-space closed form.

**On the collapse problems** both programs land within a few percent of the
exact values — on the same side for Davis & Booker, on opposite sides for
Cox. That pattern is worth stating plainly: against an exact limit load, a
few percent of discretisation-dependent bias, of either sign, is the honest
result class for a displacement finite-element method, and neither program's
number should be read as more than that.

**What "identical" means here.** Identical problem definitions and reference
values, not identical meshes or solver settings; the PLAXIS numbers are the
manual's published results, quoted unchanged. KATAI's numbers come from the
checked-in input files reproduced above, and are re-asserted (with declared
bands) by `test_input_corpus` on every build. Solve timings in the
transcripts are from an ordinary 4-core desktop and vary with hardware; the
solved numbers do not — the result contract is byte-identical across the
front ends (CLI, Python) and across the solver-backend compositions.

**Scope.** Four problems verify four specific answer classes — elastic
settlement (homogeneous and non-homogeneous incompressible), and collapse
loads (axisymmetric c–φ and heterogeneous Tresca). No broader claim is made
from them; the full [verification matrix](verification-matrix.md) records
what else is covered, case by case.
