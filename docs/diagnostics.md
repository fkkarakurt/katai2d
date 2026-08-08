# Run diagnostics

A `.k2d` project is checked twice, against two different things.

The **input contract** (`katai validate`) reads the file: units, ranges, references, the rules a
project must satisfy to be a project at all. It runs before a mesh exists, so there are questions it
cannot answer — above all, *does this object reach the soil?* A line load drawn five centimetres
above the ground surface is a perfectly valid load object. It is also a load that will never touch
the model.

**Diagnostics** are the second report. They are raised while the project is turned into an analysis,
when the mesh is there to be asked, and they carry one rule:

> An input may be used differently from the way it was written, but never in silence.

Two severities follow from that rule, and the difference between them is consequence, not rarity:

| Severity | Meaning | What the run does |
|---|---|---|
| `warning` | The analysis answers your question, but the model is not literally the drawn one — a line was clipped to the soil it could reach, a load acts at the nearest node. | continues |
| `refusal` | The object would contribute **nothing**: a surcharge that lands off the mesh, a wall the mesh never sees. | stops |

An unattached object is refused rather than warned about on purpose. If the run continued, it would
converge, report a plausible displacement field, and answer a model you did not draw — and in a
model with self-weight, the settlement would hide the missing load completely. A wrong answer that
looks right is worse than no answer.

Every code below is **stable**. Match on the code, never on the message: the prose may be improved,
the code is never reworded, reused or retired into something else.

## Codes

| Code | Severity | Raised when |
|---|---|---|
| `K2D-G001` | refusal | A **point load** lies outside every element of the mesh. |
| `K2D-G002` | warning | A **point load or anchor end** attaches to a node more than a quarter of that element away from where it is drawn. The mesher refines around such a point but does not insert it as a vertex, so a small snap is ordinary discretisation — this says how large it is. |
| `K2D-G003` | refusal | A **distributed load** line does not resolve to a chain of mesh edges: drawn above the surface, outside the model, or across a hole. |
| `K2D-G004` | warning | A **distributed load** is applied over only part of its drawn length; the message gives both lengths. Global equilibrium reflects the applied part, not the drawn one. |
| `K2D-G005` | refusal | A **plate** does not lie on the mesh. |
| `K2D-G006` | warning | A **plate, geogrid, wall or interface** is built over less than its drawn length, because the rest falls outside the soil. Its forces are those of the shorter element. |
| `K2D-G007` | refusal | A **geogrid** does not lie on the mesh. |
| `K2D-G008` | refusal | An **anchor** has neither end in the soil, or collapses onto a single mesh node. |
| `K2D-G009` | warning | A **wall or interface** line does not lie on mesh edges, so the mesh cannot be split along it: the element acts bonded to the soil, with no slip and no gap. That is stiffer, and generally less conservative, than the model drawn. |
| `K2D-G010` | refusal | A **wall's toe** — the deeper end the plate hangs from — is not a mesh node, so neither the plate nor its interfaces can be built. Usually the toe was drawn outside the soil. |
| `K2D-G011` | refusal | A **prescribed displacement** line catches no node, so it would impose nothing. |
| `K2D-G012` | refusal | An **embedded beam** has no material assigned, or a non-positive length or diameter. |

Codes in the `K2D-M…` family are about a **material parameter** the selected model does not read;
`K2D-A…` about a limit of the **analysis** that shapes how a result may be read.

| Code | Severity | Raised when |
|---|---|---|
| `K2D-M001` | warning | A **tension cut-off** is switched on for a Hardening Soil / HS-small / Soft Soil / Soft Soil Creep material. Only the Mohr–Coulomb return reads it in this build, so the run allows tension past σ_t. The schema, like PLAXIS, switches the cut-off on by default, so this is a systematic difference from the reference code in the unsafe direction. |
| `K2D-A001` | warning | A **linear Dynamic** phase reports the stress field as zero everywhere — the linear path never recovers it. Displacements, accelerations and structural forces are that phase's results; for stresses, run it nonlinear. |
| `K2D-A002` | warning | A **linear Dynamic** phase carries structural elements: geogrids stay elastic, anchors do not yield, interfaces do not slip. Their capacity is not checked during the shaking. |
| `K2D-A003` | warning | A **prescribed displacement drives a node a structure stands on**. Structural elements do not receive the imposed motion in this build, so that element's M, Q and N understate the action. |
| `K2D-A004` | note | A **structural element ends on a supported node**. The reported reactions are the soil's contribution only; the element's own end force at that support is not included. |
| `K2D-A005` | warning | A **Safety (strength reduction) phase** runs with a non-associated flow rule (ψ < φ). The factor of safety then depends on the mesh and **falls as the mesh is refined**, because failure localises into a shear band whose width is set by the elements. Measured on the Griffiths & Lane benchmark: −7.9% across a fourfold refinement. Quote the factor with the mesh it was computed on, and confirm it with a refinement study (`docs/validation/numerical-uncertainty.md` §5). |

## Where they appear

**Command line.** Under the phase that raised them; on a failed run, under every phase that ran:

```
phase 1/2: ok  max|u| = 2.145e-02 m
  warning K2D-G004  Strip q: Distributed load "Strip q" is applied over 2 m of the 7 m drawn:
  the rest of the line falls outside the soil.
```

**Python.** On every result, and on the refusal that stopped a run:

```python
job = katai.run(project)
for r in job.results():
    for d in r.diagnostics:
        print(d.code, d.message)

try:
    job = prj.run()
except katai.Refusal as e:
    if any(d.code == "K2D-G003" for d in e.diagnostics):
        ...   # a distributed load was drawn off the soil
```

## Why the list is short

It is short because it is complete for what it covers: every place in the analysis driver where a
drawn object can fail to become an engine object. That completeness is held by a source gate
(`scripts/check_silent_drop.py`), which fails the build when a build loop skips user input without
either raising a diagnostic or stating why nothing is lost — so the list grows with the code rather
than drifting behind it. The runtime behaviour is pinned by `tests/test_no_silent_drop.cpp`
(verification case `KV-DIA-001`), one fixture per code.
