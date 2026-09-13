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
| `K2D-G012` | refusal | An **embedded beam** has no material assigned, a non-positive length or diameter, or a hinged connection point that is not on the mesh at all (its head would be tied to nothing). |
| `K2D-G013` | warning | An **embedded beam** with a hinged connection is tied to the nearest mesh node, and that node is more than a quarter of the local element size from the pile head as drawn. The pile is connected where the message says, not where it was drawn; refine the mesh at the head if that distance matters. |

Codes in the `K2D-M…` family are about a **material parameter** the selected model does not read;
`K2D-A…` about a limit of the **analysis** that shapes how a result may be read.

| Code | Severity | Raised when |
|---|---|---|
| `K2D-M001` | note | A **tension cut-off** is switched on for a Hardening Soil / HS-small material, and is applied *sequentially* — to the principal stresses that model's own return produced (Material Models Manual Eq. 3-11) — rather than as one coupled multi-surface solve. The correction is compressive, so it moves away from every tensile surface; the one it can in principle disturb is the volumetric cap, which is not re-checked against the capped stress. Soft Soil and Soft Soil Creep read the cut-off on the same footing and raise nothing, having no such cap. |
| `K2D-M002` | note | An **undrained Hardening Soil / HS-small material** takes its pore-fluid stiffness from the unload/reload pair (`Eurref`, `nu_ur`) at the reference pressure. That is the model's own elasticity — the `E`/`nu` boxes it never reads no longer size the water — but it does not follow the stress-dependent E_ur(σ₃) during the run. |
| `K2D-M004` | warning | An **HS small material asks for G0 more than 20x its unload/reload shear modulus G_ur**. The model permits at most that ratio (PLAXIS MMM §7.5), so the run uses the capped G0 and says which value it used. |
| `K2D-M003` | note | A material is **Undrained (C)**: a total-stress analysis. `E`/`nu` are the undrained pair, `c` the undrained shear strength with φ = 0, no pore pressure is generated or carried, and its K0 refers to total stress. The stresses reported for it are **total** stresses and cannot be compared with the effective stresses of a Drained or Undrained (A)/(B) region in the same model. |
| `K2D-M005` | note | A phase asks for the **small-strain history to be reset** (`phases[].resetsmall`, PLAXIS's "Reset small strain", MMM §7.6). The note says how many stress points were cleared, so the phase meets the soil at G0 rather than at the stiffness the earlier phases degraded it to; stress, shear hardening and the preconsolidation pressure are carried over untouched. If the model contains **no** HS-small material, the note says so instead: the flag was read and honoured, but there was no history to clear. |
| `K2D-A001` | warning | A **linear Dynamic** phase reports the stress field as zero everywhere — the linear path never recovers it. Displacements, accelerations and structural forces are that phase's results; for stresses, run it nonlinear. |
| `K2D-A002` | warning | A **linear Dynamic** phase carries structural elements: geogrids stay elastic, anchors do not yield, interfaces do not slip. Their capacity is not checked during the shaking. |
| `K2D-A003` | warning | A **prescribed displacement drives a node a structure stands on** with a **non-zero** value. Structural elements do not receive the imposed motion in this build, so that element's M, Q and N understate the action. A component prescribed to zero is a rigid support, not a motion: it understates nothing and does not raise this warning. |
| `K2D-A004` | note | A **structural element ends on a supported node**. The reported reactions are the soil's contribution only; the element's own end force at that support is not included. |
| `K2D-A005` | warning | A **Safety (strength reduction) phase** runs with a non-associated flow rule (ψ < φ). The factor of safety then depends on the mesh and **falls as the mesh is refined**, because failure localises into a shear band whose width is set by the elements. Measured on the Griffiths & Lane benchmark: −7.9% across a fourfold refinement. Quote the factor with the mesh it was computed on, and confirm it with a refinement study (`docs/validation/numerical-uncertainty.md` §5). |
| `K2D-A006` | warning | A **Safety phase is asked for a tolerated error looser than the search's own 1e-3**. The strength-reduction search reads "this trial converged" as "the slope stands", so a loose stopping rule makes it stand at strengths it cannot carry — and the error is one-sided: the factor comes out **too high**. Measured on the Griffiths & Lane benchmark: +2.0% at 1e-2, +45.6% at 1e-1 (KV-NUM-007). |
| `K2D-A007` | note | A phase applies only **part of its staged change** (Σ Mstage < 1). The configuration the phase describes is not reached — the remainder is still carried by the soil, and the results belong to the partial stage. |
| `K2D-A008` | note | A phase **ignores undrained behaviour** for a material that declares Undrained (A) or (B): it is solved drained, so no excess pore pressure is generated in it. Strength parameters are unchanged. The answer is long-term (or state-setting), not short-term. |
| `K2D-A009` | warning | A **well is active in a consolidation or fully-coupled phase**. Those solvers take drainage boundaries and drains, not a prescribed discharge, so the well's pumping is not applied in that phase. Drains, which set the excess pore pressure to zero, are. |
| `K2D-A010` | warning | A **consolidation or fully-coupled phase runs with walls or interfaces in the model**. Those phases do not read the cross permeability (`flow_barrier`) yet, so water crosses the line as if the soil were continuous and a cut-off wall holds back less head than it would in the ground. The steady groundwater-flow calculation does read it. |
| `K2D-A011` | warning | A **flow field computed with a barrier is handed to a deformation phase**. The head is discontinuous across the barrier, but the pore load is applied from one side of it, so the differential water pressure on the wall is not in the structural forces. Read the two sides from the flow result (`head`, `head_far`) meanwhile. |
| `K2D-A012` | warning | The **constitutive integration ran out of substeps** in one or more committed increments. The Hardening Soil integrator subdivides an increment until its own error estimate is under the integration tolerance, but it may take at most a fixed number of pieces; when it runs out it returns the best it has. No equilibrium residual can see this — the stresses the check compares are both produced by the same cut-short walk — so the phase can report "converged" while its path was integrated more coarsely than asked. Raise the guard (`KATAI_HS_MAXSUB`) or loosen the phase's integration tolerance until the warning stops. |
| `K2D-A013` | warning | A **consolidation phase asked to run until the excess pore pressure reaches a target met that target at its first time step**, because the stage generated no excess pore pressure above it. The time reported is then the size of that first step, not a consolidation time. The usual cause is a phase whose load, fill or excavation was left switched off, so the stage changes nothing. |
| `K2D-A014` | warning | A **consolidation phase solved its structural elements elastically and one of them passed its capacity**. Plates, anchors and geogrids take part in the coupled solve, but with no anchor yield, no geogrid tension cut-off and no plate hinge — so a force above the capacity the engineer entered is reported as if the element could carry it, and the redistribution that yielding would have caused is missing from the whole result, the soil included. The warning names the line and its utilisation. Check the stage in a Plastic phase before using the numbers. |
| `K2D-A015` | warning | A **geogrid in a consolidation phase is in compression over part of its length**. A geogrid carries tension only — in compression it goes slack and carries nothing — but the structural branch of a consolidation phase is elastic and has no such cut, so those stations pushed back on the soil instead. The error is on the UNSAFE side: it stiffens ground the real sheet would have stopped holding. Measured on the verification fixture at 1.57% of the sheet's force with compressed ends, and 0.00000% without (KV-STR-006). |
| `K2D-A016` | warning | An **interface in a consolidation phase passed its Coulomb capacity**. The structural branch of a consolidation phase is elastic, so the joint could not slip — it carried the shear instead. A joint that cannot slip is stiffer than the real one, so the wall it holds deflects less and attracts more load than it would: the error is on the unsafe side. The warning names the joint, the peak demand/capacity and the fraction of the length over it. Check the stage in a Plastic phase, where the Coulomb return runs. |
| `K2D-G014` | refusal | An **interface takes its strength from a Hoek-Brown material**. An interface is Mohr-Coulomb whatever the surrounding model, and it is normally given c_i = R·c′ and φ_i from the soil beside it — but the Hoek-Brown model reads neither c′ nor φ′, so those boxes hold whatever was left in them and the joint would silently be given a strength unrelated to the rock it is cut into. Point the interface at a Mohr-Coulomb material (`iface_material`) whose c′ and φ′ ARE the joint strength intended — for a rock joint, the discontinuity's own friction, not the rock mass's envelope. |
| `K2D-G015` | refusal | A **material-factored design approach (EC7 DA1-C2 / DA3) was asked for on a model containing a Hoek-Brown material**. The partial factors divide c′ and tan(φ′); this model has neither, so the rock would be solved at its CHARACTERISTIC strength while the report said the design approach had been applied — a wrong verdict, not a conservative approximation. EN 1997-1 gives no partial factor for a Hoek-Brown envelope. Apply the factors to an equivalent c′/φ′ (MMM Eq 4-15/4-16) on a Mohr-Coulomb material, or use a resistance-factored approach (EC7 DA2, TBDY 2018), which does not touch the material. |
| `K2D-G016` | refusal | A **structural element is active in a Safety (strength reduction) analysis** — the initial Safety procedure or a Safety phase. The strength-reduction search solves the ground alone: a plate, anchor, geogrid or embedded beam contributes neither its stiffness nor its weight, and an interface, which splits the mesh along its line, leaves the soil on its two sides unconnected. Measured on the Griffiths & Lane slope, an active geogrid, anchor, plate carrying 150 kN/m/m and embedded beam each returned the factor of safety of the same mesh with the element deactivated, bit for bit. The error can go either way — dropping a weight that loads the slope raises the factor, dropping a member that holds it lowers it. Deactivate the structural elements in the Safety phase to obtain the factor of safety of the ground without them. An interface cannot be deactivated per phase in this build, so a model that contains one cannot run a Safety analysis yet. |
| `K2D-A017` | warning | A **Hoek-Brown material was left on the automatic K0**, which is Jaky's 1 − sin(φ′). That model has no φ′, so the formula reads the material's unused friction box and the K0 it produces is an accident rather than a measurement. The value it lands on for an untouched box (K0 = 1.00, lithostatic) is what a competent rock mass is usually assumed to be in, so the run is not refused — but a jointed or stress-relieved mass can be far lower, and that is a decision the engineer makes, not a default. Switch the automatic K0 off and enter the value the site investigation supports. |

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
