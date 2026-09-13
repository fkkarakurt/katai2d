# Changelog

All notable changes to KATAI 2D. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
MAJOR.MINOR.PATCH.

## [Unreleased]

### A Safety run answered for a model without its structures

A Safety analysis — the initial Safety procedure or a Safety phase — finds its factor of safety by
re-solving the ground under reduced strength, and that search was never handed the structural
elements. Nothing said so. Measured on the Griffiths & Lane slope, with each element active and
then deactivated on the same mesh:

| element in the Safety run | factor of safety, active | factor of safety, deactivated |
|---|---|---|
| four geogrid layers, EA = 1e5 kN/m | 1.0185791015625 | 1.0185791015625 |
| two anchors across the slip surface, EA = 1e6 kN | 1.0103271484374998 | 1.0103271484374998 |
| a slab on the crest, w = 150 kN/m/m | 1.0319091796875 | 1.0319091796875 |
| an embedded beam through the slope | 1.0103271484374998 | 1.0103271484374998 |

Bit for bit, displacement field included, on both linear-solver backends and in a chained Safety
phase as in the initial procedure. The slab's weight was verified to be applied in full in a
static phase (the vertical reactions rise by exactly w × length), so what disappeared was the
element, weight and stiffness together. The error runs in whichever direction the element acted:
dropping a weight that loads a slope raises the factor, dropping a member that holds it lowers it.
Two elements also behaved differently by backend — the plate and the embedded beam add degrees of
freedom that nothing stiffened, so the Eigen build refused every trial and reported an unstable
slope while the MKL build answered silently — and an interface, which splits the mesh along its
line, left the soil on its two sides unconnected: a joint as strong as the soil turned a slope that
stands at FoS 1.02 into "did not reach equilibrium even at the lowest strength factor" on both.

**A structural element active in a Safety run is now refused** (`K2D-G016`), by the input contract
at `initial.struct` / `phases[i].struct` and again by the engine, with no factor of safety
reported. The remedy is to deactivate the elements in the Safety phase, which runs and gives the
factor of safety of the ground without them. An interface cannot be deactivated per phase in this
build, so a model that contains one cannot run a Safety analysis yet; the test that pins the
refusal also pins that sentence, so it fails the day it stops being true.

The Python `prj.phases.safety()` docstring said "phi-c reduction of the current state". It is not:
the search starts every trial from an unstressed state, so the stresses the earlier phases left are
not its starting point. The docstring now says what the phase solves.

### An anchor on a driven line reported no force, under a warning that blamed the wrong thing

`K2D-A003` warned that a structural element standing on a node driven by a prescribed displacement
does not receive the motion, so its M, Q and N understate the action — read the soil, not the
structural diagram. The analysis stopped being that way on 2026-08-13, when every structural
element loop started reading a driven node's displacement, and the warning kept firing anyway,
including on `KV-STR-005`'s own run, whose driven geogrid force that case asserts at +0.000%.

Taking the warning apart found it half stale and half true, and the true half was not what it
said. The plate, geogrid, embedded-beam and interface force reports read the full displacement
vector, whose fixed entries carry the prescribed values. **The anchor force report did not**: it
still skipped fixed degrees of freedom, a rule its own comment attributed to the solver the
solver had since dropped. Measured before the fix, against closed forms:

| structure on driven nodes | reported | closed form |
|---|---|---|
| two-span plate, middle support settles 10 mm (force method, with the plate's shear term) | M_B = 8.988161933073 kNm/m | 8.988161933064 kNm/m |
| strut between an edge held at u_x = 0 and an edge driven by −10 mm | N = 0.0 kN | −100 kN |
| fixed-end anchor from the driven edge to a point 5 m beyond it | N = 0.0 kN | +200 kN |

The anchor report now reads what the solver reads, and **`K2D-A003` is retired**: no longer
raised, kept in the catalogue so it is never given another meaning. A fixed degree of freedom that
is a support holds zero, so every anchor not standing on a non-zero prescribed displacement
reports exactly what it did. The three closed forms are a new verification case, `KV-STR-009`.

Five places repeated the warning's account of a plate standing on a uniformly pushed line as
"undriven" — two solver comments, two tests, and `docs/validation/numerical-uncertainty.md` §11.7.
Their numbers were right: such a plate carries no moment. The reason was not — every one of its
nodes takes the same settlement, so it translates without curving. Corrected, and the published
record carries a dated correction note rather than a silent edit.

### A flow-barrier warning fired on the seams it said were not read

`K2D-A010` warned, in every consolidation and fully-coupled phase whose model contained a plate or
an interface, that those phases do not read cross permeability and water crosses the line as if the
soil were continuous. Since `KV-STR-007` that is not true of a split seam, and the warning fired
anyway — on the permeable default too, where water crossing is exactly what was declared. Measured
with a surcharge on one side of the line, the largest pore-pressure difference across it:

| the line | consolidation (Tv = 0.05) | fully coupled (Tv = 0.5) |
|---|---|---|
| interface, impermeable | 12.62 kPa | 1.91 kPa |
| wall with interfaces, impermeable | 12.57 kPa | 1.71 kPa |
| the same, fully permeable | 0 | 0 |
| plate without interfaces, impermeable | pore field bit-identical to the permeable plate's | the same |

So a split seam is read in both phases, and the one barrier that is not is a declared barrier on a
line the mesh was not split along. The warning now says exactly that and fires only there — an
active plate or interface with `flow_barrier` ≠ 0 and no seam (a plate without interfaces, or a
wall or interface that fell back to bonded, `K2D-G009`) — naming the line. Nothing pinned the old
behaviour; `KV-STR-007`'s test now pins the new one in both phases, including that the
fully-coupled phase reads the seam at all, which no test had checked.

### Two small statements that were not true

- The refusal of **Hardening Soil (and HS small) with Undrained (C)** advised "Use Undrained (A) or
  (B) with this model" — and Undrained (B) is refused for that model too, by the check right after
  it. It now advises Undrained (A), says (B) is not supported for it either, and a check in
  `test_material_registry` keeps the advice to what the model accepts.
- The README's table of the Python surface listed six material constructors and left out the
  seventh, `prj.materials.hoek_brown`, which shipped in 0.9.0.

### A correction to 0.9.0: the Hoek-Brown Safety refusal gave the wrong reason

The 0.9.0 entry *Rock, from the file to the answer* says a Safety phase on a Hoek-Brown material is
refused because no rule for reducing its shear strength was known, and the refusal message said
that only the tension cut-off's reduction is defined for this model. That was wrong. Strength reduction
is defined for this model by writing the Hoek-Brown yield function with the reduction factor inside
it (Benz, Schwab, Kauther & Vermeer 2008, *Int. J. Rock Mech. Min. Sci.* 45(2), 210–222), evaluated at each stress point's own σ′₃, so it
needs no σ′₃max. The factor that formulation gives does not correspond to the safety factor of a
Mohr-Coulomb material with equivalent strength properties — which the refusal's recommended
remedy, an equivalent Mohr-Coulomb fit, did not say.

The refusal itself stands: this build does not implement that formulation, and without it a Safety
phase would keep the rock at full strength. Its reason is corrected, in the engine and in the input
contract, and the remedy now says what it is — an approximation whose factor belongs to the fit
and to the confining range it was made over, not the Hoek-Brown factor of safety. The released
0.9.0 entry keeps the reason it gave at release; this note is the correction.

### Every reference is a primary source

A statement in this project now stands on a formula written out, a closed form or an academic
primary source, and on nothing else. Citations of other programs' manuals are gone from the source,
the messages, the Python docstrings, the validation record and this changelog. A pointer to
another manual's equation number is removed; where the equation matters it is written out or
attributed to the paper it comes from (the Hoek-Brown constants to Hoek, Carranza-Torres & Corkum
2002, the strength-factorised form to Benz, Schwab, Kauther & Vermeer 2008). `test_product_text`, which
read only the command line and the package docstrings, now reads every tracked file of the tree,
and the Studio's gate its whole tree. The staged-construction multiplier is described in this
tree's own words — the stage fraction, `mstage` — in messages, docstrings and the Studio's phase
panel (`K2D-A007` now reads "mstage = 0.5"); the file key and the Python attribute `sum_mstage`
are unchanged.

Most verification cases lost nothing, because the quantity was already asserted against its own
closed form and the removed check was a second comparison with another program's published output:
the Giroud rigid footing (`KV-FND-001`, `KV-FND-012`, against 15.15 kN/m), the Davis & Booker
strip footing from its file (`KV-FND-014`, against 7.80 kPa), the
sliding block (`KV-STR-002`, against the Coulomb closed form) and the beam-bending pair
(`KV-STR-003`, against the Timoshenko closed forms, 13.96206 and 17.43428 mm).

Two cases changed what they claim, and one record was tidied:

- **The Gibson strip load (`KV-FND-002`, `KV-FND-011`) is a plausibility band now, and says so.**
  Its only reference for the 4 m layer on a rigid base was a published finite-element number; the
  closed form is for the half-space, and it is exact only for a fully incompressible soil. Measured
  on the corpus geometry scaled in depth, the model settles −9.2% (4 m), −1.3% (8 m) and +3.4%
  (16 m) of the half-space value — a finite layer settles less, and ν = 0.495 settles more as the
  layer deepens — so the half-space number is not this model's answer either. The cases assert
  85% to 100% of it (measured 90.2% and 90.8%), which a stiffness profile read 10% too stiff fails,
  and the record states that it holds no independent reference for the finite layer.
- **The homogeneous slope (`KV-SLP-001`, `test_slope`, `test_safety_gui`, the Python example) is
  measured against the referee value 1.00** of the simple slope in Giam & Donald (1989), *Example
  problems for testing soil slope stability programs*, Monash University report 8/1989, in place of
  a 0.99 that had no single primary source. The computed factor is unchanged at 1.010
  (+1.0%), and so is the 8% band. The file carries γ = 20.2 kN/m³ against the problem's 20.0, which
  can lower the factor by at most 1%.
- **Four rows of the verification matrix no longer run on** into the section title that follows
  their declaration in `test_foundation_benchmarks`.

### Upgrading from 0.9.0

- A project whose Safety phase (or initial Safety procedure) has a structural element active is
  refused. Before this build it ran and reported the factor of safety of the same model without
  that element; deactivate the element in the Safety phase to get that number knowingly.
- Scripts that matched on `K2D-A003` will no longer see it. An anchor whose end stands on a
  non-zero prescribed displacement now reports its force; before, it reported 0.
- Two corpus files are renamed, `tests/corpus/kv-str-002-sliding-block.k2d` and
  `tests/corpus/kv-str-003-beam-bending.k2d`, and so are the examples the installer ships from them.
  The test target that runs the four footing benchmarks is now `test_foundation_benchmarks`, and
  the circular-footing study `study_circular_footing`. Contents and assertions are as described
  above.
- The validation comparison page for the four footing benchmarks is withdrawn; each case is in
  `docs/validation/verification-matrix.md` against its analytical solution.

## [0.9.0] - 2026-08-31

Two things this program could not be given, and one thing it could not say.

It could not be given **rock**: a rock mass had to be entered as a Mohr-Coulomb fit, which is a
straight line through a curve — matchable over a narrow band of confining stress and wrong outside
it in both directions. It now has the Hoek-Brown 2002 criterion, entered the way a geologist writes
it down (σci, mi, GSI, D), integrated in plane strain and axisymmetry, and verified through a
boundary value problem rather than only at a material point.

It could not be given **the ground data itself**: every geotechnical model starts from borehole
logs, and this program could only read polygons, so the logs had to be redrawn by hand — the
slowest step in setting a model up and the easiest place in the whole workflow to put a number in
the wrong row, which is the one error nothing downstream can catch. Logs are now an input, and the
polygons are generated from them.

And a run could not say **what its answer had actually satisfied**. "Converged" was a single global
force residual, and the two criteria that complete it — the current stiffness
parameter that tightens the test as a mechanism forms, and the local error at the stress points,
which asks whether the stress a point carries is the stress its own material law would return —
had never been consulted. Both are now computed, both are reported, and a result file carries which
of them were met.

Around those: a wall may now stay in the ground while the ground consolidates and keep its
interfaces while it does; axisymmetric ground may have water in it; a consolidation phase can be
asked *how long until 90%* instead of being told a duration and asked to guess; the constitutive
integration tolerance became something a file can state rather than something the reader's build
chose; and three guards that were right but silent — the mesher's step cap, the integrator's
substep cap, and a stalled Newton search reported as a collapse mechanism — were each given a wire
out of them.

### Upgrading from 0.8.1

- **Project files are now `.k2d` version 18** (from 14), through four steps that each added keys
  rather than changing existing ones: **v15** the constitutive integration tolerance
  (`phases[].substol`), **v16** borehole logs (`strata[]`, `boreholes[]`), **v17** the
  consolidation stop criterion (`phases[].cstop`, `cminp`, `cdeg`, `cfirst`, `cmaxstep`) and
  **v18** the Hoek-Brown rock model (`materials[].model` = 6 with `sigci`, `mi`, `gsi`, `hbD`,
  `sigpsi`). Every one of them is written **only when used**, so a project that touches none of
  these features produces the same bytes it did under 0.8.1 and its diff stays about what actually
  changed. This build reads every older project file. Older builds **refuse** a file written by
  this one at the version gate — which is the point of the gate: an older build cannot read
  `model = 6`, a stop criterion or a written-down integration tolerance, and would otherwise run a
  different problem from the same drawing and report it as an ordinary result.
- **Results files are now `.res` version 9** (from 6), adding which convergence criteria a run met
  (v7 soil, v8 interfaces / coupling springs / embedded-beam foot force) and the consolidation stop
  record (v9). Same rule as always: this build reads older result files, older builds refuse this
  one's. Re-open a project and re-calculate, or keep the older build for older result files.
- **`SolveResult.stopped_by` renamed one of its values.** A solve abandoned because no step along
  the Newton direction reduced the out-of-balance force used to answer `'mechanism'`; it now
  answers `'stalled line search'`, because that ending also happens on a well-supported model and
  the old name was a claim rather than a description. Code that tests `stopped_by == 'mechanism'`
  to decide whether a load factor is a capacity must be updated — read the current stiffness
  parameter beside it, as `katai.summary()` now does.
- **A Hardening Soil answer may move slightly** against 0.8.1 on the same file. The stress-point
  integrator no longer holds the stress-dependent moduli at the state each increment began from
  (measured on the oedometer: +0.41% to −0.93% against the closed form, the closed form being the
  fixed point), and an increment abandoned for a stalled search is now retried with a non-monotone
  window instead of being cut back. Both change the path, both were adopted against a measurement,
  and both are recorded below.

### The rock model meets a boundary value problem, and three things break

A model verified at the material point is not a model that works. Hoek-Brown had been checked twice
against its closed forms and once through the whole FE path from a `.k2d` file, and all
three of those tests are PLANE STRAIN and all three are element tests. Putting it in an
axisymmetric boundary value problem instead — a circular tunnel unloaded into a rock mass, against
the closed-form ground reaction curve — broke it three separate ways.

**One: the model was never integrated in axisymmetry at all.** `integrate_point_axisym` dispatches
on a closed enum, its Hoek-Brown branch had simply not been written, and the switch has no
`default:`. So a rock material in axisymmetry left `trial` and `tangent` EXACTLY as the caller
passed them — the previous iterate, handed back as though it had been integrated. Nothing failed
loudly; the equilibrium iteration merely never converged and the run reported a stalled Newton
search, which is also what a genuinely difficult problem looks like.

The branch is one paragraph. The lasting fix is `test_material_axisym_coverage`, which walks the
model REGISTRY — not a list maintained beside it — and pushes every registered model through both
integrators with the outputs **poisoned with NaN first**, so a branch that does not exist cannot
pass by handing back the caller's own values. It also asserts that the axisymmetric branch MOVES
the hoop stress, which a branch that quietly forwarded the committed value would fail. The gate was
checked by removing the new branch again and confirming it goes red.

**Two: the tangent was not strong enough for a boundary value problem.** The model handed the
solver the elastic operator rather than a consistent one. That integrates the material correctly —
every material-point test passes with it, because those tests ask what stress comes back — but it
does not iterate an ill-conditioned boundary value problem to equilibrium. The tunnel stalled as
soon as the plastic annulus passed about 4% of the radius, and load steps did not buy it back: 150
steps per stage reached exactly the same wall as 40. Newton with a wrong-but-symmetric operator
converges linearly, and linearly is not always convergent. The fix is a finite-difference
consistent tangent in both integrators, the same device the Hardening Soil and Soft Soil branches
already use, and cheap here because this return is a bisection on one scalar rather than a
substepped walk.

**Three: the retry that uses that tangent never reached the rock.** The solver's hybrid strategy —
when an increment fails with the cheap tangent, re-enter THAT increment with the stronger one — was
armed by a flag that tested for one model **by name**. So the tangent added above existed and was
never called. The evidence was itself the diagnosis: after adding it the run stalled at the
IDENTICAL load factor as before, 0.750, and at 40 steps as at 150. By the solver's own reasoning
("a limit load does not move with the load-step count, and a stall does"), a number that does not
move at all means the new code is not on the path. The flag's name was the tell — it asked "is
Hardening Soil present?" where the question is "is there a material here that can be asked for a
stronger tangent?", and a predicate that tests for an INSTANCE where it means a PROPERTY is wrong
the first time somebody adds a second instance. It is now `has_fd_tangent` and reads the question
it asks. Soft Soil and Soft Soil Creep also build a consistent tangent by finite difference and are
deliberately NOT enrolled: the retry latches for the rest of the phase, so it would change the
iteration path of runs that recover anyway, and there is no failing soft-soil case to measure that
against.

**The verification the repair earned** is `KV-CST-016`: a circular tunnel in a Hoek-Brown rock mass
under 28 MPa of hydrostatic ground, unloaded in stages, against the closed-form elasto-plastic
solution. The radial stress lands on the closed form to **0.59% inside the plastic annulus and
0.72% outside it**, and the case asserts two things on purpose — that the unloading CONVERGES
(the guard for defects two and three) and that where it converges it is right (the verification).
Either alone would have passed at some point during the repair.

The oracle is derived in the test from radial equilibrium rather than quoted, and checked by
reducing it at `a = 1/2` to the scaled forms of Carranza-Torres & Fairhurst (1999) — two
independent routes to the same expressions, so neither is a transcription of the other. The rock
mass and stress state are case D1 of Table A1.1 in Hoek, Carranza-Torres, Diederichs & Corkum
(2008); computing the rock-mass constants from GSI and D reproduces that table's printed
`m_b = 1.093`, `s = 0.0031`, `a = 0.507`.

Two of this test's own defects are recorded in it rather than quietly fixed. The relief loads
summed to less than the in-situ stress, so the tunnel started under-supported and the last stage
was taking the wall to zero rather than to the pressure it was compared against — which looked
exactly like a convergence problem. And the band on the plastic annulus was, for one revision,
**an assertion that could not fail**: the last stage stopped where the annulus reached 1.427 R
while the nearest interior node sat at 1.50 R, so the worst deviation inside it was zero by vacuum.
The test now checks that it samples the region it is asserting about.

`python/examples/rock_tunnel_ground_reaction.py` runs the same problem as a published ground
reaction curve — nineteen converged support-pressure stages with the closed form beside them.

### Rock, from the file to the answer

The rock model reached the ground. The previous increment built the Hoek-Brown criterion at the
material point; a criterion nothing can select is a header, so this one puts it in the schema
(`.k2d` v18: `sigci`, `mi`, `gsi`, `hbD`, `sigpsi`, written only by a material that uses the
model), in the material registry, in the Python surface, and in the Studio's material editor,
which asks for it in the geologist's own vocabulary and shows the derived m_b, s, a and
the resulting rock-mass strengths next to the four inputs.

**The verification is a triaxial test run by the program itself** (`KV-CST-015`): a `.k2d` file
naming Hoek-Brown, meshed, brought to an isotropic cell pressure and squeezed past yield, twice, at
1 and 8 MPa. The stress the block carries on the plateau is the Hoek-Brown envelope at the confinement
it is actually under, to **−0.0000%** at both, with the specimen 0.05% from homogeneous. Two cell
pressures rather than one because the envelope is CURVED: their ratio is **3.231** where a straight
line through the origin would give 8.000, so a build that had quietly fallen back to Mohr-Coulomb
would miss the second point even if it had been tuned to hit the first.

**The fixture is the finding.** Driven by LOAD it reproduced the envelope to −0.99% at the high
cell pressure and only −5.72% at the low one, and the obvious suspect was the load step. It was
not: halving the increment from 395 to 182 kPa moved the error from −5.65% to −5.72%, which is to
say not at all. A homogeneous specimen of a perfectly plastic material reaches the envelope
everywhere at once and has no reserve anywhere to redistribute into, so a load-controlled Newton
can only approach that limit from below and stops at whichever increment it last equilibrated.
Under displacement control the same state is an equilibrium the solver can stand on, and the answer
is read from the stress field instead of inferred from a load factor. The band went from 3%, all of
it used, to 0.1% with none of it used.

**FIVE PLACES WERE ABOUT TO DIVIDE A STRENGTH THAT IS NOT THERE.** Making a model reachable means
every surface that assumed `c'` and `phi'` now meets one that has neither, and the failure mode is
always the same shape: the surface goes on working, touches nothing, and reports as though it had.
Four of them are refused and one is warned about, each where it happens and each saying what to do
instead.

  * **Undrained (B) and (C)** enter the strength as an undrained shear strength su and put it in
    `c'`. The entry would have been taken and then ignored, and the run would have used the FULL
    drained rock envelope while the engineer believed su was in force. Undrained (A) stays: it
    changes no strength, adds the pore fluid's stiffness and lets the criterion act on effective
    stress, which is what it is written in.

  * **A Safety phase** reduces strength by dividing `c'` and raising `tan(phi')`. The rock would
    have kept full strength through every trial of the strength-reduction search, and the search
    would have walked to its cap and reported the cap itself — **"FoS > 3.0", for any rock mass
    whatever, however weak** — a wrong number in the safe-looking direction, which is the worst
    kind this program can produce. It is refused rather than approximated because no rule for
    it was known: the one thing Safety was known to do to this model is reduce the tension
    cut-off value, which says nothing about the shear strength. The conversion to an equivalent
    Mohr-Coulomb pair does exist, and reducing THAT pair would be the natural construction, but
    it is a fit over a confining range whose upper limit depends on the application and is left
    to the caller — so a factor of safety built on it would be a function of a number nobody
    entered, wearing the clothes of a measurement.

  * **A material-factored design approach** (EC7 DA1-C2, DA3) is the same defect at a worse seam:
    the partial factors divide `c'` and `tan(phi')`, so the rock would have been solved at its
    CHARACTERISTIC strength under a report saying the design approach had been applied. A design
    verification that quietly used unfactored strength is not a conservative approximation, it is
    a wrong verdict. EN 1997-1 gives no partial factor for a Hoek-Brown envelope. The
    resistance-factored approaches (EC7 DA2, TBDY 2018) never touch the material and still run —
    which the test checks, because a refusal that also blocks the working path is a different bug
    (`K2D-G015`).

  * **An interface** takes `c_i = R·c'` and `phi_i` from the material beside it. Beside rock those
    boxes hold the schema DEFAULTS — 1 kPa and 30 degrees — so the joint would have been given a
    Mohr-Coulomb strength with no relation to the rock it is cut into. The remedy already existed
    in the schema and is what the refusal names: point the interface at a Mohr-Coulomb material
    (`iface_material`) whose `c'` and `phi'` ARE the joint strength intended, since an interface
    is Mohr-Coulomb whatever the surrounding model. For a rock joint
    that is the discontinuity's own friction, not the rock mass's envelope (`K2D-G014`).

  * **The automatic K0** is Jaky's 1 − sin(phi'), and with no phi' it reads the unused box and
    lands on K0 = 1.00. That one is NOT refused: lithostatic is what a competent rock mass is
    usually assumed to be in, so the value is defensible — but it is arrived at by accident, and a
    jointed or stress-relieved mass can be far lower. The run now says so (`K2D-A017`).

And one that is neither refused nor warned but simply **fixed**: the local convergence criteria
normalise a stress difference by the material's cohesion, and `cohesion_of` would have returned
ZERO for rock. Its own comment already names that failure — "a normaliser that is quietly a factor
too small makes every point look accurate, which is the one failure mode a convergence check must
not have" — and it was about to happen on a material whose strengths are megapascals. Rock now
returns the criterion's own strength at zero confinement on a cohesion scale, sigma_c/2, the Tresca
relation the model's degeneration was already verified against.

**One tension control was implemented rather than refused.** A user-entered tensile strength may
cap the criterion's own tensile strength: where that value is lower than σ_t, the tensile capacity
is cut off at it. The tree already has the two schema fields every other model's Rankine cap
uses, so the whole implementation is one line — the cap is applied to `Constants::sigt`, and every
path that reads a tensile limit (the yield functions' second branch, the apex return region, the
tensile branch of the mobilised dilatancy) is written in terms of it, so none can miss it. It only
lowers: a value above σ_t would claim a strength the criterion has not got and is ignored. The
consequence is stated where it is visible rather than left to be discovered — this schema's default
for that field is ON with σ_t = 0, so a rock material carries no tension unless the box is
unticked, and the Studio shows the control for rock instead of hiding it, next to the σ_t the
criterion would otherwise have. Hiding it would have been the neater panel and the worse one: a
control that is absent still acts.

For the same reason the text report, the HTML report and the model summary print the four rock
inputs where they used to print `c'` and `phi'` — numbers the calculation never read.


### Ground data arrives as borehole logs, and until now it had to be redrawn as polygons

Every geotechnical model starts from the same document: a log per borehole giving the level of each
layer boundary and the water table at one position. This program could not read one. The engineer
turned the logs into polygons by hand — the slowest part of setting a model up, and the easiest
place in the whole workflow to put a number in the wrong row. A number in the wrong row *here* is
the one error nothing downstream can catch: the mesher will mesh the wrong ground, every phase will
converge on it, and the answer will be a correct solution to a model nobody meant.

**The rules are the established ones**, because the borehole workflow is the one every user of
this kind of program already knows: the layer list is **global** — every layer exists at every
borehole, and a layer absent somewhere is not a missing row but two equal levels — and a single
log makes a horizontal water surface that reaches the model boundaries, while several combine into
a non-horizontal one. Between logs the boundaries interpolate linearly; **outside** the
outermost log its levels are *held*, never continued on their slope, because extrapolating grows
ground nobody logged.

**Boreholes generate polygons; they do not replace them.** `polygons` remains the model — what the
mesher meshes, what a phase activates — and the logs are the record of where that geometry came
from, so a reviewer can see the source rather than only the drawing. Generation is therefore an
**action the user takes**, never something a run does on the way past: a geometry with two sources
of truth is a geometry that can disagree with itself. No solver path is touched by this; the
generator is a pure schema-to-schema function, testable with no mesh and no solve.

`KV-GEO-001` checks the *rule*, not the output: 48 assertions whose expected levels are computable
by hand from the logs in the test — one flat log reaching both edges, two logs interpolating
between and holding outside (20 and 18 read 20 and 18 at the edges, **not** 21 and 17), and a
three-layer section whose middle layer pinches out to exactly 0.000 m with the layer beneath rising
to meet the one above, so the ground has no gap. The refusals are asserted too, and a failed
generation leaves the project untouched rather than half-applied. GEO joins DIA as a matrix class
on DIA's own argument: it verifies what the *pre-processing* produces, and the failure it guards
against is invisible to every other check in the suite.

`.k2d` v16 adds `strata[]` and `boreholes[]`, written only when there are boreholes, so every older
model is byte-identical to what v15 produced. The Python surface exposes both types — the
schema-coverage gate is what insists on that, and it is right to.

### The last coupled family carries the structure too

A wall could stand in the ground while the pore pressure dissipated, but not while the water table
moved: the consolidation phase carried plates, anchors, geogrids and interfaces, and the
fully-coupled phase refused all of it. That is not a distinction the physics makes, and it is not
one the two solvers make either — they take the same structural stiffness and the same tied seam
pore equations, and differ only in what the fully-coupled phase adds on its own (the retention
curve and Bishop's chi, which the structure neither sees nor is seen by).

**The oracle is that sameness, measured.** The fully-coupled core reduces to the consolidation core
in the saturated limit, so the same model with the same raft, run once each way, must give one
answer: settlement **4.50097547e-02 m from both (+0.00000%)** and a raft moment of **21.621032
kNm/m from both** — and from the drained Plastic phase they both end in. That is what a reduction
should look like: the two solvers are running the same system, and what separates them is their own
iteration, a Picard fixed point against a Newton one, which here separates them by nothing
measurable (`KV-STR-008`).

The three findings the elastic structural branch owes the reader — a plate or anchor past its
capacity, a geogrid pushed into compression, an interface past its Coulomb capacity — are now
written **once** and called by both phases rather than copied into the second one. Two copies of a
statement are two statements, and one of them decays first.


### The wall may keep its interfaces while the ground consolidates

The previous increment let plates, anchors and geogrids into a consolidation phase and refused
interfaces, on the grounds that the split seam "would silently make the joint impermeable, which is
a modelling claim". That was right about the danger and wrong about the remedy: the model already
carries the interface's cross permeability, and **its default is fully
permeable** — the flow net runs through the line. The fix was to read the input that was already
there, not to keep refusing.

So a joint now carries what its own flag says. **Fully permeable** ties the two sides of the seam
to *one* pore equation, which makes continuity and the flux balance across the joint hold by
construction rather than by a constraint that could be assembled slightly wrong.  **Impermeable**
leaves the two pressures the split produced. **Semi-permeable** is refused, and the message says
what it is — a conductance q = dh/R — rather than that it is unsupported, because rounding it to
either neighbour would be silent and directional.

The flag is what the test measures, early in the dissipation where there is a difference to see:
across the seam, `|p_left − p_right|` is **0.000e+00** with a permeable joint and **12.83 kPa** (of
a 50 kPa surcharge) with an impermeable one. Measuring at the *end* of the phase, where everything
has drained, would have passed for a build that ignored the flag entirely — which is how the first
version of the test was written, and why it did not stay that way.

Mechanically the joint comes with the same limit as the rest of this phase's structural branch —
it is elastic and cannot slip — and the same treatment: measured, not declared. One below its
capacity reproduces the drained Plastic phase in both the settlement (**−0.0042%**) and the shear
the joint carries (**53.2745 vs 53.2745 kPa**). One above it raises **`K2D-A016`** naming the joint
and its exceedance, and says which way the error runs: a joint that cannot slip is stiffer than the
real one, so the wall deflects less and attracts more load (`KV-STR-007`).

Two things the work turned up on its own:

- **the excess pore pressure field was never reported.** `pore` on a result is the *hydrostatic*
  pressure the phase was set up in; the pressure a consolidation phase computes existed only as a
  maximum per time step, so an engineer could see that something was still draining but not where.
  It is now on the result as `excess_pore`, node by node, and readable from Python;
- **the drained-limit comparison is only clean on an elastic soil.** With Mohr-Coulomb the same
  pair differs by 0.15%, and it is not the joint — the drained reference does not slip either way.
  It is the soil's stress path: the coupled phase loads it undrained and the Plastic phase drained,
  and a yield surface remembers the difference. That number is measured and printed beside the
  assertion rather than folded into a band wide enough to hold it, because such a band would also
  be wide enough to hide a joint that was not in the system at all.


### Axisymmetric ground may now have water in it

An axisymmetric model with a water table was refused — *use plane strain, or remove the water
table* — which closed every circular problem that has groundwater in it: the tank, the silo, the
shaft, the circular footing, the pile load test. It now runs.

What was missing was two assemblies and **one term inside one of them**. The K0 seed was already
effective-stress and already set the hoop; the phreatic body force only needed its r weight. The
pore-pressure load needed something the plane-strain version has no place for: in plane strain the
load is ∫Bᵀ(u·m)dA with m = [1, 1, 0] and the out-of-plane direction carries no equation, but in
axisymmetry the strain has four components, the hoop is a real strain with a real equation, and
pore pressure is isotropic — so **m = [1, 1, 0, 1]**, and the hoop term lands on the *radial*
degree of freedom beside ∂N/∂r. A careful-looking copy of the plane-strain function would drop it,
and dropping it produces no error message and no obviously wrong picture.

So the test is built on the identity it breaks. A K0 state is a state of equilibrium, so the
internal force of the seeded stresses must equal the phreatic body force plus the pore load — and
that is algebra, not a solution. Measured on the assembled vectors:

| | residual | without the hoop term |
|---|---|---|
| 6-noded | 5.26e-03 | 3.04e-01 (58×) |
| 15-noded | **8.13e-13** | 3.07e-01 (4e+11×) |

The identity returns to round-off the moment the quadrature can integrate it — the r weight makes
the radial integrand a degree higher than the plane-strain one, and the 3-point rule of the
6-noded triangle cannot take a cubic exactly, which the driver already knew and is why it refuses
to read a nil-step from an assembled axisymmetric imbalance. That is what tells the tri6 residue
apart from an imbalance, and it is what makes the hoop term worth eleven orders rather than a
correction.

End to end, the same cylinder gives the buoyant effective stresses and the hydrostatic pore
pressure at every depth, and **lowering its water table by 3 m settles it 4.171531e-03 m against a
closed form of 4.173557e-03 m (−0.05%)** while moving radially by 0.08% of that. The K0 phase's own
displacement is deliberately *not* asserted: on level ground with a level table nothing is ramped,
so a zero there is arithmetic rather than evidence (`KV-CST-013`).

Structural elements in axisymmetry are still refused, and the message now says why it is not a
matter of effort: a plate in axisymmetry is a shell with a hoop membrane force and an anchor is a
ring, so they are different elements rather than the same ones integrated differently.


### The wall can stay in the ground while the ground consolidates

A consolidation phase used to be soil-only: any structural element made it refuse, so the analysis
an engineer most often wants — what the excavation does over the months while the excess pore
pressure dissipates — could not be run with the structure holding it up. **Plates, anchors and
geogrids now take part in the coupled solve.** Their stiffness enters the same system as the soil
and the water, assembled by the same function `solve_nonlinear` uses, so the wall in a consolidation
phase is the wall in a Plastic phase rather than a second implementation of one.

**The oracle is the limit the coupled solution has to walk into.** As the excess pore pressure
vanishes, the coupled problem *becomes* the drained problem — which this program already solves by
a completely different path. So the phase is run to Tv = 4 and compared with the drained Plastic
phase of the same model: settlement **−0.0025%**, plate moment **0.0000%**, anchor force
**−0.00055%**. Removing the plate moves the same settlement by **17.66%**, three orders outside
that band, which is what makes the agreement evidence rather than two soft numbers agreeing
(`KV-STR-006`).

**Two limits are reported, not declared.** The structural branch here is elastic, and what that
means differs by element:

- a plate or an anchor is elastic *until* it hinges or yields, so the limit is a capacity: a line
  past the capacity the engineer entered raises **`K2D-A014`** with its utilisation, computed in
  the units the capacity was entered in (per anchor, not per metre of wall — the conversion is a
  documented trap and the test pins it);
- a geogrid is different, and the measurement said so. Tension-only *is* its behaviour, so where
  the settlement bowl puts a sheet in compression the elastic branch has it push **back** on the
  soil instead of going slack — stiffening ground the real sheet would have stopped holding, which
  errs on the unsafe side. It raises **`K2D-A015`** naming the stations. The size is measured:
  −1.57% on a sheet with compressed ends, and **0.00000%** on the same sheet kept wholly in
  tension, which is what proves the cause is the compression cut and not the coupling.

**What is still refused now says why, and the reason is not effort.** An interface (and the
embedded wall built from one) splits the mesh, so the two sides of the joint would carry separate
pore pressures with nothing between them — silently making the joint impermeable, which is a
modelling claim, not a default. An embedded beam's skin resistance follows the effective stress
that consolidation is busy changing, so an elastic spring through it is a stronger claim than the
same spring in drained ground. Both are their own work items.


### A refused linear solve ended the process instead of the time step

Found by the test above on its first honest run, and older than it. When the coupled tangent of a
Biot solve goes singular — a soil body that reaches its strength under the load being consolidated,
or a model that is not restrained enough — the linear backend verifies its answer, sees a relative
residual of 9e-4 where 1e-6 was asked, and **refuses** rather than returning a vector that does not
satisfy the system. The static solver has caught that refusal since it began checking its answers,
and treats it as a property of the increment. Neither coupled core did: the exception escaped the
call stack and terminated the process (0xC0000409).

Both now end the time step with it and report `converged = false`, which the phase turns into the
message it already had — *the load increment may exceed the soil capacity*, which is exactly what
happened. Only `SingularSystem` is caught; a malformed request or a broken backend still propagates,
because turning one of those into "did not converge" would publish a modelling answer for a bug.
The case is pinned in `test_consolidation_stop`: a confined column dissipating more excess pore
pressure than its Mohr-Coulomb strength can carry must come back and say so.

### A consolidation phase can be asked how long, instead of being told

Until now the only thing a consolidation phase could be given was a duration. The design question is
the other way round — *how long until the excess pore pressure has gone?* — so answering it meant
guessing a span, reading the curve, and guessing again. A phase can now end when the ground gets
there instead: **at a maximum excess pore pressure**, or **at a target degree of consolidation**. The time
interval is then not used at all, and the time the target took is what the phase reports.

**The degree of consolidation is a pressure ratio, and every surface that prints it says so.** The
criterion's "degree of consolidation" is defined as the excess pore pressure left over the
maximum the stage generated, not the settlement ratio the name suggests — and the two are different
numbers, not two spellings of one. On the 1-D column where both are known in closed form, the
settlement ratio reaches 90% at Tv = 0.848 and the pressure ratio only at Tv = 1.031: **21.6% apart
in time**, and at the pressure criterion's 90% the settlement is already 93.6% done. A build that
implemented the name rather than the definition would be wrong by a fifth of its answer and would
look entirely reasonable doing it, so the definition travels with the number in the solver message,
the report (text and HTML), the Studio help, `katai.summary` and the Python result.

The march that finds the time changes its step size as it goes, and it does that by **restarting the
fixed-dt core** rather than by growing a new one — both cores are restartable by construction, and
the test asserts it (one run of 40 steps against four of 10: 6e-16 relative). Its two constants were
measured on the column whose answer is known, not chosen:

- the reported time is first order in the steps-per-doubling (8 → +5.5%, 32 → +1.5%, 128 → +0.4%),
  and 32 ships. The axis is not the mesh: a fine equal-step grid gives the same time on 99, 169 and
  453 nodes alike;
- the automatic first step is **four times** the Vermeer–Verruijt critical step, because dt_crit is
  a stability bound on the pore field and not an accuracy bound on the pressure the ratio is
  measured against — at dt_crit exactly, the first step overshoots the undrained pressure it
  generates by 12.8%, and that peak is the ratio's denominator.

Measured against the closed form the stop time is **+1.45%**; the two criteria agree with each other
to 0.00% where they ask for the same instant, as does the elastoplastic path (`KV-CON-003`).

Three endings that could have been silent are not. A march that runs out of steps **refuses**,
saying how far it got in the definition it was asked in, rather than reporting where it happened to
stop. A model with no drainage boundary is refused before the march starts, so a structural defect
is not reported as a step budget. And a stage that generated no excess pore pressure at all — a
phase whose load was left switched off — meets its target at the first time step; the run raises
**`K2D-A013`** rather than presenting the size of that step as a settlement time.

Python: `prj.phases.consolidation("...", until_degree=90)` or `until_excess_pore=1.0`, with
`first_step=` and `max_steps=`; passing both a target and `duration=`/`steps=` is an error rather
than a silent preference. File formats: **.k2d v17** (`cstop`, `cminp`, `cdeg`, `cfirst`,
`cmaxstep`, written only when a criterion is set, so older files are byte-identical) and **.res v9**
(what the phase was asked to end on, whether it got there, and the ratio's reference).


### The mesher's safety valve opened in silence, which is what a safety valve must not do

Ruppert's refinement guarantees the minimum angle it is asked for, and this tree asks for 20°,
under the ~20.7° the termination proof needs — so mesh quality here is **structural** rather than
measured element by element, which is the right design: a bound that holds by construction beats a
bound checked afterwards. The loop carries a step cap anyway, as a final safety valve, and until
now that valve had no wire out of it. `refine()` returned void. A mesh that met 20° and a mesh that
ran out of steps trying were the same object to every layer above — builder, driver, report, GUI —
and mesh quality is not something an answer is indifferent to.

`Triangulation::quality_met` / `refinement_steps` and `MeshResult::quality_met` carry it out, with
the bound that was asked. A mesh that did not reach its bound now says so in the message every
front end already shows, and says what to do about it; a mesh that did reach it says nothing extra,
because a warning that fires on the ordinary case is a warning nobody reads.

`KV-DIA-002` checks the flag **against the geometry**, not against itself: a flag that is always
true proves nothing, so the test computes the smallest interior angle over every produced triangle
from the vertex coordinates, with no help from the mesher, and asserts that report and measurement
agree (128 elements, worst angle 45.0000° against a 20° bound, 73 refinement steps).

**Honestly not shown:** the valve was not made to fire — 200 000 refinement steps is not a number
an ordinary geometry approaches — so the claim rests on the guard rather than on an observed
failure, and the test says so rather than implying a demonstration it does not contain.

This is the third instance of one pattern in this tree, which is worth naming: **the guard was
right, the silence was the defect.** `K2D-A012` fixed it for the material law, this fixes it for
the mesher, and the abandonment classification below fixes the same shape for the solver's verdict.
What they share is that the code already knew and had nowhere to put it.

### A run says when its material integration ran out of room

The Hardening Soil integrator subdivides a load increment until its own error estimate is under the
integration tolerance, and it is allowed a limited number of pieces to do it in. When it ran out it
returned the best it had — correctly, that is what a guard is for — and then **the flag saying so
was dropped between the material and the solver**. Nothing downstream could tell a run that met its
integration tolerance from one that did not, and no equilibrium check can recover the difference:
the residual is assembled from the very stresses the cut-short walk produced, so it balances
perfectly around them.

The flag now reaches the result. A phase whose committed path contains such an increment raises
**`K2D-A012`**, naming how many increments were affected and the worst number of stress points
involved, and the count is on the result for a script to read
(`convergence.saturated_increments`, `convergence.saturated_points`). Runs that met their tolerance
are unchanged and say nothing, which is the point.

### Half of what produced an answer could not be written down

A `.k2d` has carried its numerics since v7: the tolerated force residual, the load increments, the
iteration limit. All three govern the **equilibrium** iteration. The other half of what produces an
answer — how accurately each stress point is walked along its material law *inside* an increment —
could only be set through an environment variable, so a project handed to a reviewer described the
model, described the stopping rule, and silently left the constitutive integration to whatever the
reader's build happened to choose. That is not a small omission on this tree: the integration
tolerance is where a first-order error was hiding that no equilibrium residual could see, and its
1e-5 default was set by measurement rather than by taste.

`phases[].substol` (`.k2d` v15), threaded the way the other three are — the jobs-layer seam wins,
then the file, then the material class's own default, and 0 everywhere means "the class chooses".
It reaches the material routine as an **assembly** parameter rather than a stopping rule, because
that is what it is: it is read inside the Gauss loop, beside the creep time interval, not by the
iteration around it. Models without an error-controlled integrator ignore it rather than reject it.

It is asserted hostilely, because the failure mode is unusually quiet — a dropped integration
tolerance produces a perfectly convergent run with a slightly wrong stress path, which no residual
reports and no plot shows. File and seam must agree **bit for bit** (one control, two routes) and
both must differ from the default, so that "it is read at all" is proved rather than assumed:
0.018700624 m against 0.018643176 m, 0.3081% apart.

**`KATAI_HS_STOL` now wins over the file**, which is the opposite of the precedence the file has
over the class default, and deliberately so: the environment variable is a whole-run study
override, and a study exists to ask what a published number owes to a numerical choice — it cannot
ask that of the files that state the choice if the file overrides it. Same precedence and same
reason as `KATAI_CONV_NOLOCAL` over the phase's convergence setting.

### Two convergence criteria had never been consulted; both have been now

The criteria family added above was measured and reported, but two of its members had no case
behind them. Both have been read through to the end, and neither answer was the expected one.

- **The moment residual was not actually in the gate**, although the record said it
  was. It is now — and measuring what that is worth found a real defect. With an *elastic* plate the
  rotational equations are linear, so the linear solve satisfies them exactly and the criterion
  reads round-off at every tolerance, including a run whose wall deflection is 21% wrong: it is
  answering a question about the rotational equations, not about the answer. With a *plastic* plate
  it comes alive and follows the tolerance down, but never exceeded 6% of the force error over 15
  load/tolerance pairs. Binding it therefore costs nothing on anything this program runs today, and
  covers the case it does not yet have — a structure that fails while the soil around it is still
  elastic, where the force balance is scaled by ground that is not being asked for much.
- **That binding cost exactly one case, and the case was right.** A plate standing on a line that is
  pushed down is undriven — the program already says so — and an undriven plate carries no moment,
  so the criterion was dividing one round-off by another and refusing the analysis outright. Every
  other criterion in this family has a floor for exactly that situation; this one did not, because
  nothing had ever consulted it. It has one now (1 kNm/m), which is ten decades above the undriven
  case's own reference and two decades below a loaded plate's, and both ends are pinned by tests.
- **The non-linear elastic criterion is identically zero in this program, and now it is
  known why.** Run on unloading — the only place a non-yielding point can have a stress-dependent
  stiffness — all 96 points are counted and the error is 1.3e-15 at every tolerance. The unloading
  modulus is held at the state each increment begins from, so a point that does not yield walks the
  increment with a constant elastic operator and the two stresses the criterion compares are built
  from the same arithmetic. The check is therefore a check of that identity: if the modulus ever
  starts changing inside an iteration, this is the first count that moves.

**A run also now reports both global force ratios**, the CSP-normalised one and the fixed-scale one
that actually decides when a step stops. Only the first was published before, so a run could print a
force error next to a tolerance it was never compared against.

### A run now says which convergence criteria it met, not just that it "converged"

The solver checked one thing — a global force residual against a fixed scale — and reported the
result as though that one thing were the whole question. It is not: the criteria are a family,
and the members disagree.

A run now measures and reports all of these at the iterate it accepted, per phase:

- **The Current Stiffness Parameter (CSP)**, the ratio of the work an increment actually did to
  the work the same strain would have done had the response stayed elastic. It is 1 while the
  model is elastic and falls towards 0 as a mechanism forms — measured on one strip footing at
  1.00000, 0.437, 0.098, 0.00012 as the load rises past what the soil can carry.
- **A CSP-normalised global force error.** Because CSP falls as the body plastifies, this
  criterion *tightens* as a mechanism forms — which is exactly where a load fraction is about to
  be read as a bearing capacity. The old fixed normalisation does the opposite.
- **A moment residual**, wherever something in the model carries a rotational degree of freedom.
- **The local error at every soil stress point**, plastic points and stress-dependent-elastic
  points counted separately. A stress point carries two stresses during an iteration: what the
  material law returns for the strain it was given, and what the finite-element linearisation says
  it carries. They coincide only at the solution, and their difference is invisible to a global
  force residual by construction — the residual is assembled *from* those stresses.

**What the measurement found.** On a Mohr-Coulomb strip footing the global force error falls by
five orders of magnitude across the iteration while the worst local error falls by a factor of
three; at the accepted iterate the two differ by 63 000x. On the Hardening Soil oedometer at the
tolerance this program ships for that model family, the global-only stopping rule stops **0.18%
away from its own converged answer** — and requiring the local criteria at the *same* tolerance
lands on the converged settlement in 194 iterations where tightening the global tolerance by four
decades costs 399.

**They now decide whether a step has converged.** A run must satisfy the local criteria as well as
the force balance before an increment counts. The decision was taken on the measurement rather than on the principle: this is not a
stricter rule bought with iterations, it is a cheaper route to the same answer — the Hardening Soil
oedometer lands on its converged settlement in 194 iterations where tightening the global tolerance
by four decades costs 399. `KATAI_CONV_NOLOCAL` turns it off for a run, which is how every
comparison in the record is reproduced.

**Three published numbers moved, out of 154 checks, and the record re-measures all three.**

- **A slope's factor of safety is no longer inflated by a loose stopping rule.** A
  strength-reduction search asks whether a reduced strength still reached equilibrium, and a run
  that stops on the force balance alone can stop with its stress points nowhere near the strengths
  it has just reduced them to — which the search reads as a yes. At a hundredfold looser stopping
  rule the reported factor of safety used to come back **45.6% too high**; it now comes back 0.6%
  high. A slope reported 45% safer than it is was a number this program should never have been
  able to produce.
- **A creep case's published 3% agreement turns out to have been luck.** Requiring the local
  criteria moved every duration of the Soft Soil Creep column away from the idealised creep law.
  Sweeping the tolerance with the criteria switched off showed why: the model's own converged
  answer is +5.93% / +3.14% / +1.64% from that law at 1, 10 and 100 days, and the previous run had
  simply stopped short at a place that happened to sit inside 3%. The case now carries a 7% band
  and, more usefully, asserts the SHAPE — the deviation must fall as creep comes to dominate,
  which a run drifting for a numerical reason has no reason to do over three decades of time.
- **A sliding block's failure force is a plateau to nine significant figures instead of to the
  last bit.** Doubling the imposed slip moves it by 1.4e-9 relative, against the ~100% a stiffness
  reading would move.

**Two structural criteria complete the family.** A slipping interface point and an embedded beam's
skin coupling springs are measured the same way and counted together; the pile toe carries its own
out-of-balance ratio, tolerated at five times the tolerated error. On a sliding-block interface the
global criterion is met by a factor of five while 14 of 47 slipping points are not yet settled.

**Results files carry all of it.** The `.res` format moves to version 8 so that a reopened result
can still say which criteria its numbers were accepted under. Files written by earlier versions
read back with the family marked "not measured", which is the honest answer for them; older builds
refuse a version 8 file rather than mis-read it.

**Two orientations of CSP, resolved by measurement.** CSP can be written as either of two
reciprocals — elastic energy over total, or total over elastic. Only the second is consistent with
the behaviour the parameter exists to have (unity when fully elastic, approaching zero at failure)
and with what is built on it. The first form is at least 1 and grows without bound as a mechanism forms,
which would make the global criterion loosen towards collapse. The second is implemented, and a
test pins the direction so the other reading cannot return quietly.

### The Hardening Soil answer no longer depends on how many load steps it was asked for

The stress-point integrator held the model's stress-dependent moduli — `E_ur`, `E_i`, `q_a`,
`q_f`, every one of them a function of the current confining stress — fixed at the state each
load increment began from, and used them for the whole increment. That is a first-order error in
the increment size, and no substep tolerance can see it: the integrator reported that it had met
its accuracy target while the answer was still moving with the step count. Measured on one
oedometer path at a fixed integration tolerance, halving the step halved what was left, all the
way down: the "converged" answer was **4.5% low at 20 steps and still 0.6% low at 160**.

The moduli are now read at the state each substep begins from, which is what the model says they
are. The substepping itself is now error-controlled as its citation always claimed — a
modified-Euler pair measures the local error, and the subdivision is sized from it against a
declared tolerance (default 1e-5) rather than a fixed fraction of a reference stress.

What that changes, on the case the verification record publishes (a laterally confined Hardening
Soil column, `KV-CST-002`):

- **The load path stops being an error axis.** A 16× refinement of the increments moved the answer
  by 2.9 percentage points; it now moves it by **0.020%**.
- **The two linear-solver backends agree.** The same run split between PARDISO and Eigen by
  0.18 percentage points; the three stress ranges now agree to **15 significant figures**.
- **The remaining ~1% is the model's, and it can be shown without a finite element run.** The same
  calibrated material, integrated at a single stress point with no mesh, no load path and no
  equilibrium iteration, lands within 0.3 percentage points of the boundary-value answer.
- **The cap calibration is more accurate**, because it runs through this same integrator: it now
  reproduces its own `Eoed_ref` and `K0_NC` targets to better than 0.1%.

**Hardening Soil results will differ from 0.8.1 and earlier.** They differ because the earlier
ones carried an error that was invisible to every control the program offered. Two published
numbers moved in the process, and the record says so rather than quietly restating them: the
conclusion that `KV-CST-002`'s residual deviation was a cap calibrated at `p_ref` — argued from a
signature that grew with stress level and changed sign — is **retracted**, because that signature
was the frozen-modulus error. `docs/validation/numerical-uncertainty.md` §6 and §7 are rewritten
around what is measured now.

### The non-monotone line search, adopted where it is needed and nowhere else

The comment beside `NewtonOptions::line_search_window` had been describing this defect since it was
written: a monotone test reads a non-descending step as failure and halves the increment, and four
such halvings in a row abandon it — "which is how the **load path** stops being the one the file
asked for". It was a prediction, and the refinement ceiling measured a few commits earlier is its
instance: every abandonment was `four consecutive iterations without descent`, at increments
repeatedly halved, with the linear solver never once refusing and the iteration budget never
reached (21/15/31/64/59/40 out of 500). Refining further did not help, because **size was never the
obstruction**.

The mechanism was already implemented and switched off, so this is an *adoption*, and the only
question was where to switch it on. Three designs; two killed by measurement.

1. **Make it the default.** It removes the refusal — but the suite said no: binding the local
   convergence criteria at the shipped tolerance stopped reaching the answer (0.05% → 0.147%
   against the four-decades-tighter run). A non-monotone rule accepts iterates a monotone one
   rejects, so the stopping test fires further from converged, and that price is paid on every
   increment including the overwhelming majority that never stall.
2. **Escalate mid-iteration** — open the window on the fourth iteration without descent and carry
   on from the iterate that stalled. Rescued nothing. The window works by changing the whole
   trajectory, not by recovering one; by the time four steps have failed, the iterate is somewhere
   a wider gate cannot come back from.
3. **Adopted: retry the increment.** An increment abandoned for a stall is re-entered with the
   window open and its size untouched — the same shape as the hybrid tangent immediately above it
   in the same function, and for the same reason: when an increment cannot be closed the cheap way,
   try the stronger tool on *that* increment rather than paying for it everywhere.

Measured on `KV-CST-002`'s seating phase at 240 increments — refused under the monotone rule,
`-0.923430%` with the window on throughout, `-0.923431%` with retry-on-stall (eight figures) —
while the bound run's accuracy claim comes back to 0.0005%. `recent` now keeps the largest window
an increment could use and consults only its tail; with only `ls_window` entries kept an escalation
would arrive with no memory and be monotone for another five iterations, which are exactly the five
that were failing. That was measured before the line was written.

**What it costs, because it is not free.** `test_phase_numerics_009` goes from 781 s to about
2320 s and is once again the suite's longest test; the full run goes 2141 → 2335 s. The extra time
is real work — increments that used to be abandoned and cut back are now retried and converge.

### A stalled search was being published as a bearing capacity

A run that stops below full load has to say which of three things happened, because the same
fraction means different things and only one of them is a capacity. That distinction was made once
already, for the run that merely exhausts its iteration budget. It was not made for the second
case, and the second case is the common one.

When no step along the Newton direction reduces the out-of-balance force — the search stalls — the
program reported "a collapse mechanism formed (the remaining load exceeds the soil capacity)" and
published the equilibrated fraction as the incremental limit load. In four places: the message, the
Studio's Phases panel and Output window, and both reports. But that ending happens with the tangent
still non-singular — the linear solver answered every time, which is exactly what separates it from
the singular-tangent ending — so a well-supported model reaches it too. Measured on a laterally
confined weightless column, which has no mechanism available at all: the same file reported 68%,
74% and 83% of the load and on a fourth run carried all of it, decided by nothing but the load-step
count and which linear solver ran it. A capacity is not a function of either.

**The fix is a discrimination, not a retreat.** The first attempt was to stop claiming a capacity on
that ending, and it was wrong: `KV-FND-010`'s Prandtl strip footing, whose limit load is verified
against `N_c = 2 + π`, ends the very same way. What separates them is the current stiffness
parameter — **0.00010 for the footing against 0.65550 for the column**, with 1468 yielding stress
points against 96 — and the threshold is CSP < 0.5, the level below which arc-length control is
specified to engage. So the rule is now one function beside the field it reads, called by all four
surfaces, and each of them names the stiffness parameter it decided on rather than deciding
invisibly.

The asymmetry it leaves is deliberate and stated where the rule lives: a low stiffness parameter
proves the body softened into a mechanism, a high one proves nothing, because CSP is blind to
STRUCTURAL plasticity — a plate that has formed a hinge reports 1.00000. Such a run is reported as
"not established as a capacity", which is the safe side of a signal that cannot see it.

**Python surface:** `SolveResult.stopped_by` returned `'mechanism'` for this ending. The name was
the claim, so it is now `'stalled line search'`, and `katai.summary()` reads the stiffness parameter
alongside it rather than trusting the label.

### Known limits

- The integration tolerance reached its home in this release — `phases[].substol`, beside the
  tolerated error and the load-step count, so a published run now carries the accuracy it was
  computed with. What is still a build-time choice is the **default** it falls back to (1e-5 for
  Hardening Soil), set by measurement on one case rather than derived per problem.
- One increment's subdivision is capped at 200 substeps. Saturation is counted and reported, never
  absorbed; it is reached almost only on trial iterates whose stresses are discarded.
- Soft Soil, Soft Soil Creep and Mohr-Coulomb keep their own substepping rules. They do **not**
  carry this defect — both soft-soil integrators take their bulk and shear moduli from each
  substep's own trial pressure, and Mohr-Coulomb has no stress-dependent modulus to freeze.
- The Mohr-Coulomb failure bound is still a one-sided clamp on the major principal stress, which
  ratchets slightly when the failure deviator moves inside an increment (measured: a drained
  triaxial stalls 0.24% below its plateau). The consistent projection that fixes it costs a
  boundary-value problem that passes today, so it is measured, recorded and not taken.
- **Refining the load path has a ceiling, and above it whether a run converges at all is not a
  property of this program.** The claim above that the two linear-solver backends agree holds
  below that ceiling and only there. On `KV-CST-002`'s seating phase — a weightless column whose
  confining stress starts near zero, where the Hardening Soil stiffness is at its smallest — the
  measured ladder is: 40 and 80 increments converge and agree to every printed digit on both
  backends; 120 refuses on both; 160 refuses on Eigen and converges on PARDISO at 66% of the
  tolerated residual, taking twenty times as long as a rung that converges thirty times further
  inside it. A coarser path failing where a finer one passes is round-off deciding the answer,
  not refinement improving it. So the verification case asks its path-independence question at 80,
  where the answer exists in both compositions, and the ceiling itself is deliberately not
  asserted — a test that asserts a refusal is a test that fails when the program improves.
  **That ladder was measured before the retry-on-stall adopted later in this same release**, and
  the refusals on it were the abandonments that retry now re-enters: the seating phase at 240
  increments went from refused to converged, to eight figures of the always-on window's answer. The
  ladder itself has not been re-measured rung by rung, so what stands is the mechanism, not the
  numbers above 80 — and the case still asks its question at 80, where the answer existed under
  both rules. Choosing the increment size instead of being handed one is what the next release's
  automatic step control is for.

## [0.8.1] - 2026-08-19

Two things a result can be wrong about without looking wrong, and one of them had
been deciding which of the two solver backends was right.

The first is a **verdict**: a load-controlled analysis that stops short of full
load was publishing the fraction it reached as the soil's capacity, whichever of
the three ways it had stopped — including the one that is a statement about a
setting rather than about the ground. The second is a **field**: the degree of
saturation was drawn, tabulated and printed for runs that never computed it. Both
reached the calculation report, which is the document that leaves the building.

Alongside them the notation the interface and the reports write in was typeset, and
the Newton loop's remaining cost on non-smooth problems was measured and written
down rather than traded away — see *Known limits*.

### Upgrading from 0.8.0

- **Results files are now `.res` version 6** (projects are unaffected: `.k2d`
  stays at 14). This build reads every older results file; 0.8.0 and earlier
  **refuse** one written by this build, by the same forward-version guard that has
  always protected them. Re-open a project and re-calculate, or keep the older
  build for older result files.
- **The per-increment iteration limit is now derived from the material class** —
  200 where the tangent is nonsymmetric (nonlinear soil, interfaces, embedded
  beams), the phase strategy's 80 otherwise — instead of always being 80. A run
  that previously reported a collapse mechanism at a low load factor may now
  converge. To reproduce a 0.8.0 run exactly, set that phase's `maxiter` to 80 in
  the file, or "Max iterations" in the Numerical tab.
- **`τmax` is now `τmax (in-plane)`** wherever a field is named, and **`S` is no
  longer offered by runs that did not compute it** — a static analysis has no
  saturation field, and the selector no longer pretends otherwise.
- **New: `SolveResult.stopped_by`** in C++, in Python and in the `.res` file. A
  load factor below 1 cannot be read without it; see *Fixed*, first entry.

### Fixed

- **A patience setting was published as a soil capacity, and it decided which of
  the two solver backends was right.** KV-STR-004's refined-mesh probe reported
  "a collapse mechanism after 10% of the load" on the vendored Eigen solver and
  full convergence on MKL — the disagreement this project recorded as open, with
  two explanations, in `docs/validation/numerical-uncertainty.md` §10.
  Instrumenting the Newton loop refuted both: through the whole failing run the
  linear solver refused nothing, and there is no discontinuity in the collapse
  detection near the limit load. What the run hits is the per-increment iteration
  limit. The fixture is weightless with the tension cut-off on, so a point load
  makes most of the model a no-tension material — measured on the coarse mesh,
  same file, one field changed: **561 Newton iterations with the cut-off, 57
  without** (two per increment). On the refined mesh the increments need up to 96
  while their out-of-balance force falls monotonically the whole way, and the
  limit was 80: PARDISO's path needed 88 and survived by cutting back four times,
  Eigen's needed 96 and did not. With the limit derived rather than fixed, **both
  backends converge to max|u| = 0.1822167 m — seven significant figures — and the
  MKL run is 1.8× faster than it was** (37 s against 66 s), because the cut-backs
  a truncated increment triggers cost more than the iterations it was denied.
  Raising the limit is free where a mechanism genuinely forms: a collapsing run
  ends when the tangent goes singular and the solver refuses a direction, measured
  at the same limit load and the same 122 iterations under both limits.
- **One sentence was printed for all three ways an increment is abandoned.** "A
  collapse mechanism formed … this equilibrated fraction is the incremental limit
  load" is true when no descent direction exists, and true when the tangent goes
  singular; it is false when the iteration budget simply ran out, which is what
  KV-STR-004 was doing. The solver now records which of the three ended each
  abandoned increment, the phase message names it, and the iteration-budget case
  withdraws the limit-load claim and names the setting to raise. The reason
  travels with the result — `SolveResult::stopped_by`, `.res` v6 — so no front end
  has to re-derive it from the number, which is what the Studio was doing in four
  places, one of them the calculation report.
- **A field nothing had computed was drawn, tabulated and printed as a result.**
  S, the degree of saturation, is produced by unsaturated flow — a transient or a
  fully coupled phase — and by nothing else. Every other run leaves the array
  empty and the field reader answers 1.0 when asked anyway, and nothing filtered
  the ask: an ordinary static analysis offered S in the field selector, painted a
  uniformly saturated model, tabulated "min 1.000, max 1.000", and printed a
  saturation row in both reports, beside numbers the solver did produce. Above a
  phreatic surface that is not merely uncomputed; it is the one answer that is
  certainly wrong. A field is now offered only where the run produced it.
- **`tau_max` said more than it computed.** It is the radius of the Mohr circle of
  the IN-PLANE stresses; the out-of-plane σ_zz is not carried in the nodal field,
  so wherever it falls outside the in-plane pair — which plasticity can do — the
  true maximum shear is larger than the number shown.
- **The report's columns fanned out.** `printf` pads by bytes and the typeset
  names are not one byte per column; the extreme-value table is now padded by
  characters, which also fixes a Turkish material or load name having done the
  same thing since long before this release.
- **Two comments were in Turkish, in a subtree the source-language gate reports as
  clean.** They had been written double-encoded, which is what hid them: the gate
  looks for Turkish letters, and re-encoding a UTF-8 line from a legacy-codepage
  reading of its own bytes leaves none.

### Changed

- **The notation is typeset.** The interface and the reports wrote symbols the way
  a source file has to write them — `sigma'_yy`, `phi'`, `lambda*`, `kN/m3`, `>=`,
  `[deg]`. They now read `σ′yy`, `φ′`, `λ*`, `kN/m³`, `≥`, `[°]`: 225 pieces of
  user-visible text across the Studio and both report languages. The UI font loads
  Greek, the prime, super/subscripts and the mathematical operators to make it
  possible, measured against the font before it was used. Indices stay plain
  letters on purpose — Unicode has no subscript y and no subscript w, and
  subscripting most of a set while three of it cannot be is worse than
  subscripting none.
- **No other FE program is named in anything a user reads.** Twenty-two Studio
  tooltips and ten Python docstrings named one; every one of them was making a
  point that stands without the brand, and every one of them still makes it.
  Source comments and the validation record are deliberately out of scope: where a
  default came from is worth recording, and a comparison is worthless without
  naming what was compared against.

### Added

- **`SolveResult.stopped_by`** on all three surfaces, and `katai.summary` now says
  what a load factor below 1 *is* rather than only what it equals: the incremental
  limit load when a mechanism formed, and "where the iteration budget ran out, NOT
  a capacity" when one did not.
- **KV-NUM-012** (`test_non_convergence_names_its_reason`): one pile fixture
  solved past its capacity and at half of it with an iteration limit of one. The
  first must publish a limit load; the second must refuse to, because the same
  model carries that load when the budget is adequate. This is the check
  KV-STR-004 needed and did not have.
- `test_product_text` and `test_studio_product_text`: one scanner over both trees,
  reading string literals and docstrings only — it lexes past comments rather than
  splitting lines on `//`, so it sees what a user sees. Checked against a fixture
  that must fail it.
- A third detector in `check_language.py` for double-encoded UTF-8, reporting the
  line recovered rather than as stored. It found exactly the two comments above
  across 808 tracked files, with no false positives.
- Four checks pinning the saturation rule from both sides, one that measures the
  report's column alignment in characters, and a `.res` round trip for the new
  field.

### Verification

- **153/153** on the MKL composition and **151/151** on `portable` (the vendored
  Eigen solver, no proprietary component) — the first time the whole suite has
  been run on that composition rather than the subset CI builds.
- **58 declared verification cases** over 26 benchmark `.k2d` files, each with its
  oracle, locator, expected value and band stated in the test.
- The Studio's own suite (report, DXF import, product text) passes 3/3.

### Known limits

- Fifty to ninety-six Newton iterations per increment remains the cost of a
  no-tension material under a point load, and its cause is now measured: the full
  Newton step is rejected in two thirds of the iterations by a line search that
  requires the residual to fall at every single one, which a non-smooth problem
  does not do. Two remedies were implemented and measured, and **neither is in
  this release**. Levenberg-Marquardt damping of the tangent was refuted outright.
  A non-monotone line search works — 8–27% faster, the worst increment down from
  50 iterations to 37 — and was withdrawn because it moves the Hardening Soil
  numerics register: at the window that keeps KV-NUM-007 inside its published 1%,
  KV-NUM-009's residual-by-stress-range signature stops growing monotonically, and
  that signature is the evidence for its conclusion. Both measurements, and what
  adopting the change would cost, are in §10.
- A line search that finds no descent in twelve halvings still applies its last,
  smallest step — measured turning a residual of 2.57 into 534 on one iteration.
  It is recorded in §10 rather than fixed silently.

## [0.8.0] - 2026-08-17

Verification release, and the first one in which the Python surface reaches the
whole engine. Six capabilities — the interface, the plate, the geogrid, the
embedded beam, HS small, and Soft Soil with its creep — were implemented and
tested at the element, but had never been run along the path a user actually
takes: from a file, through the mesher and the driver. Run that way for the
first time, **five of the six were wrong**, and one further cross-cutting fault
(the tension cut-off) came out of the same work. Every closure arrives with a
benchmark input checked in and a stated reference for what it implements:
the verification record grows from 45 declared cases over 17 input files to
**57 cases over 26 files**, and the suite to **152 tests**.

The project file moves from version 12 to version 14. This build reads every
older file unchanged; an older build refuses these, which is the point of the
guard — each bump marks an input an older build would have dropped in silence.

### Added

- **The Python surface reaches the structural elements.** `prj.structures`
  creates walls, anchors, geogrids, pile rows and bare interfaces, one call
  each, and each returns a handle that phases activate and deactivate exactly
  like a load or a region — which is the whole of staged construction. Before
  this, a script that needed a wall had to drop to the raw bindings and assemble
  the material list, the element list and the index between them by hand, so the
  pleasant layer covered the tutorial and gave up at the first real job.
  Interfaces come from the wall that needs them (`interfaces="both"`), a wall can
  be made a groundwater screen, and a pile states how its head attaches — hinged
  by default.
- **What a structure carries is readable from a script.** `ForceStation` and
  `StructForce` are bound: the stations along an element with their arc length,
  position, N, Q, M and displacement, plus the element's name, kind, yield flag
  and extremes. A wall's bending moment and a pile's axial force existed in C++
  and in the GUI and nowhere else — which is the surface a parameter study, a
  benchmark or a thesis actually lives in.
- **Flow phases, design codes and phase water from Python.**
  `prj.phases.transient_flow(...)` and `.fully_coupled(...)`; `design=` selects
  the EC7 or TBDY 2018 partial-factor catalogue for a single phase, so a
  characteristic run and its design check live in one file; `water=` overrides
  the table or the phreatic line for one phase, which is staged dewatering.
- **`katai.summary()` and `katai.extremes()`.** A readable account of a scripted
  run: the phases named, the extremes tabulated *with the place each occurs*,
  the structural force envelopes, and every diagnostic the engine raised —
  because those change how the numbers should be read. A field that is uniform
  says so instead of inventing a location.
- **A small-strain history reset as a phase option** (`.k2d` v14
  `phases[].resetsmall`, diagnostic `K2D-M005`). A surcharge
  placed and removed to leave an overconsolidation behind also leaves a strain
  history, and in the real soil ageing erased that long before the analysis
  began. Verified by KV-CST-012 against an oracle established before the option
  existed: with the reset the run lands on KV-CST-008's fresh-K₀ answer to
  +0.0012%; without it, on a run 4.50× softer. On plain Hardening Soil, which
  has no history to reset, the flag is bit-for-bit inert and the run says so.
- **Li & Dafalias dilatancy below the phase-transformation line** for HS small,
  measured against an oracle written from the five equations and sharing no code
  with the kernel: worst difference 0.00e+00 over the branch. Recorded honestly,
  because the formulation's description disagrees with itself here — its plot of
  this function shows 1.29× the amplitude its equation gives, and a printed
  formula outranks a constant reverse-engineered from a raster plot.
- **Numerical-uncertainty cases on axes other than the mesh** (KV-NUM-009…011).
  One of them measured a prediction wrong: the Newmark scheme was argued to be
  second order and the sweep recovered three, so the rule that "the algebra
  decides the order" is now bounded rather than assumed.

### Fixed

Each of these produced a converged run, a green suite and a wrong number.

- **An interface drawn along a fixed boundary was welded shut.** Both sides of
  the split sat at identical coordinates and boundary conditions are applied by
  coordinate, so both were fully fixed: the block sheared elastically against
  its own base instead of sliding, and every check in the run reported success.
  The sliding-block benchmark returned **5,401,612 kN/m where the hand
  calculation gives 60**. Which side holds the
  support is now decided rather than copied. KV-STR-002 from the checked-in file:
  59.7202 kN/m, and deleting the interface still returns 5.4e6 — that check is
  the sentry that fails loudly if the rule is ever undone.
- **Deactivating the soil welded the beams standing in it.** The beam-bending
  benchmark is built by removing the soil cluster; every node of the
  remaining beam then touched no active element and was pinned in both
  translations. The solve converged, reported "ok" and handed back
  max|u| = 0.000000e+00 with no diagnostic. A node a plate runs through is now
  exempt; an axial-only element keeps the fixity. KV-STR-003 reproduces the
  closed-form deflections, 13.96 mm and 17.43 mm.
- **HS small rode the virgin backbone.** Masing's rule
  (γ₀.₇,reloading = 2 γ₀.₇,virgin) was not applied, degrading the stiffness twice
  as fast — measured **+5.7% / +12.9% / +34.8%** too much heave at three
  unloading sizes, a deviation that grows with strain, which is the signature of
  a wrong threshold rather than of discretisation. The ceiling on
  E₀/E_ur was missing too; it is capped now and `K2D-M004` states the G₀ the run
  actually used.
- **The tension cut-off did not reach the models it belongs to.** `K2D-M001` had
  declared since 2026-08-08 that only the Mohr-Coulomb return read it, while the
  schema switches it on by default — so every Hardening Soil, HS small, Soft
  Soil and Soft Soil Creep run in this engine allowed tension past σ_t, a
  systematic error in the unsafe direction.
  KV-CST-011 pins it, including the identity that a cut-off the tension never
  reaches is bit-identical to no cut-off at all.
- **Every embedded-beam interface spring was 2.5× too stiff.** All three springs
  must be divided by the out-of-plane spacing, as EA, EI, the weight and both
  capacities are divided; only their dimensionless stiffness factors were
  implemented. The foot also used D/2 where the equivalent radius is
  R_eq = √(12 EI/EA)/2. A pile row could not be loaded at its head, which is why
  no case had ever run.
- **Structural elements read a driven node as standing still**, so a geogrid
  under a prescribed-displacement fixture carried nothing.
- **A prescribed-displacement line set to zero raised a warning it had no
  business raising** (`K2D-A003`): nothing is understated when the imposed value
  is zero, and that is the only way this schema can support a plate at a point.
- **K₀ = 0 is accepted.** With free vertical sides it is the *only* initial state
  in equilibrium with them; negative K₀ is still refused, with a message that
  says why.
- **The installer compared a real hash against the number 55.** Without
  `-UseBasicParsing`, Windows PowerShell 5.1 returns an `application/octet-stream`
  body — which is what GitHub serves a `.sha256` attachment as — as a *byte
  array*, so splitting it on whitespace split an array of bytes and element zero
  was the first byte of the hash, as a number. The integrity check would have
  "passed" for any file whose hash began with a matching byte. All web reads now
  go through one helper that decodes explicitly, the checksum is extracted by
  matching 64 hex characters, and a checksum file that does not contain 64 hex
  characters is a hard failure — a checksum that cannot be read is not a checksum
  that passed. TLS 1.2 is forced, since 5.1 still negotiates 1.0 on some builds.
- **`CITATION.cff` carried no DOI**, so GitHub's "Cite this repository" button
  produced an entry with nothing persistent in it — the one thing a citation
  exists to provide. The archived record had a DOI all along; the file never
  named it.
- **A verification case was asserting a result that only one solver backend
  could produce.** KV-STR-004's mesh-independence probe halved the element size
  under the fixture's 1500 kN/m overload. On MKL that converges and the pile
  carries its 600.0000 kN/m capacity; on the vendored Eigen solver — same
  source, same file, same mesh — the run reports a collapse mechanism after
  equilibrating 10% of the load. Continuous integration builds the `portable`
  composition, so it had been red since 2026-08-13 while every local run was
  green. The probe now runs at 1300 kN/m, inside the range where both backends
  return the capacity exactly (measured: 1200, 1300 and 1400 all give
  600.0000 at both densities), and re-runs the coarse density at the same load
  so the comparison is like-for-like. **The backend disagreement itself is not
  closed** — it is characterised and declared in
  `docs/validation/numerical-uncertainty.md` §10, together with what was tried
  and did not work.

### Changed

- The Hardening Soil dilatancy rule had been written out four times, once per
  return mapping. It is written once now, because a rule that is right in three
  copies and stale in the fourth is a silently wrong answer on whichever stress
  path reaches the fourth.
- `docs/references/hssmall-formulation.md` is in English throughout.

### Known limits, stated rather than implied

- HS small accumulates a monotone strain-history scalar and detects no reversal
  *inside* a phase; that transformation is attributed to Benz (2006), which
  is a source to obtain rather than to paraphrase. Between phases, the
  `resetsmall` option is the remedy, and it is implemented.

## [0.7.1] - 2026-08-09

Packaging release. No engine change: the same solver, the same results, the same
verification record.

### Fixed

- The Windows executable carries its identity. It shipped with an empty version
  resource -- no product name, no version, no publisher, no copyright -- which is
  poor practice on its own and, on an unsigned and newly published binary, is
  also what a reputation heuristic reads as "anonymous": Microsoft Defender
  flagged `katai2d-0.7.0-win64.zip` as `Trojan:Win32/Wacatac.B!ml`, a generic
  machine-learning label, and removed the file. The executable now declares its
  product name, version, publisher, copyright and origin, generated from
  `version.hpp` so it cannot drift from what `katai info` prints. The rebuilt
  binary scans clean.

  This is a mitigation, not a guarantee, and the honest statement of the
  situation is in the README: an unsigned binary from a young project has no
  reputation, the published SHA-256 is how you verify what you downloaded, and
  the Python wheel -- which was never flagged -- is the alternative route.

## [0.7.0] — 2026-08-09

Capability release. Ten inputs a geotechnical model needs, and could not
express, are now in the file format; four of them closed a silently wrong
answer rather than a missing convenience. Every closure arrives with a
verification case and a stated reference for what it implements, and the
verification record grows from 27 to 45 declared cases.

The project file moves from version 8 to version 12. This build reads every
older file unchanged; an older build refuses these, which is the point of the
guard — each bump marks an input that an older build would have dropped in
silence and solved a different problem for.

### Engine — inputs that decide the answer

- **Per-material undrained stiffness** (`materials[].und_mode`, `nu_u`,
  `skempton_B`; `.k2d` v9). The pore fluid's bulk stiffness Kw/n was derived
  from a fixed nu_u = 0.495 for every undrained material — one nearly
  incompressible value applied to users who had entered something else. Either
  the equivalent undrained Poisson ratio or Skempton's B is now a per-material
  input, each converted to the same pore-fluid stiffness. The same change fixed a
  measured silent error: for the Hardening Soil family, K' was read from the
  `E`/`nu` boxes that model never reads, sizing the pore fluid by an untouched
  default — a factor of 7.6 on an ordinary data set. It now follows the
  model's own unload/reload pair, and states the remaining limitation
  (`K2D-M002`).
- **Undrained (C)** (`materials[].drainage` = 4; `.k2d` v10). A total-stress
  analysis: undrained stiffness and strength, no pore pressure generated or
  carried, K0 on total stress. Available for the Linear Elastic and
  Mohr-Coulomb models; refused by name elsewhere, and in
  consolidation and fully-coupled phases. EC7 factors its cohesion as an
  undrained strength (gamma_cu).
- **Wells and drains** (`hydros`, `phases[].hydro`; `.k2d` v11). Dewatering
  from inside the model: a well prescribes a discharge and stops at `h_min`, a
  drain holds a head (normal one-sided, or vacuum). Both switch per phase. In
  consolidation and fully-coupled phases a drain sets the excess pore pressure
  to zero and a well is reported as not applied (`K2D-A009`); in a transient
  flow phase both are refused by name.
- **Cross permeability of walls and interfaces** (`structs[].flow_barrier`,
  `hyd_res`; `.k2d` v12). A cut-off wall can block flow: the groundwater
  calculation splits its own mesh along the barrier so the two sides carry
  separate pore-pressure degrees of freedom. Impermeable, semi-permeable (hydraulic resistance d/k) or fully
  permeable — the last being the default and the previous behaviour.
- **Prescribed boundary flux** (`polygons[].edge_flux`, `Flux` in `edge_flow`;
  `.k2d` v6). Rainfall, infiltration and recharge boundaries. The kernel had
  the edge integral; no path from a file reached it, and no seepage solver took
  an external right-hand side.
- **Dilatancy cut-off** (`materials[].dilatancy_cutoff`, `e_max`; `.k2d` v5).
  A dilating soil stops dilating at its critical void ratio. Without it a dense
  sand dilates without limit and its bearing capacity comes out too high — an
  unsafe number, produced quietly. Off by default.
- **Per-phase water conditions** (`phases[].water_override`, `wx`, `wy`;
  `.k2d` v4) and **anchor prestress** (`anchors[].prestress`; `.k2d` v3).
  Between them they make an anchored, dewatered excavation expressible: the
  pit is dewatered before it is dug, and the anchors are locked off rather than
  slack. Every anchored excavation built before this was the model of a weaker
  structure.
- **Staged-construction target and "ignore undrained behaviour"**
  (`phases[].mstage`, `ignoreund`; `.k2d` v8), and **per-phase numerical
  controls** (`tol`, `loadsteps`, `maxiter`; `.k2d` v7) so that a published
  run carries the numerics it was computed with.

### Diagnostics — an input may be used differently than drawn, never in silence

- `SolveResult.diagnostics` with stable codes and a written contract
  (`docs/diagnostics.md`), surfaced by the CLI, the Python package and the
  GUI alike. Sixteen input paths that were silently discarded — a load drawn
  above the surface, a wall the mesh never sees, a prescribed displacement
  that catches no node — are now refused or reported.
- Codes added in this release: `K2D-M002`, `K2D-M003` (what a material's
  parameters are read as), `K2D-A005`–`K2D-A011` (limits a result must be
  read under, including the mesh dependence of a factor of safety and the
  one-sided hand-over of a flow field with a barrier in it).
- `scripts/check_silent_drop.py` keeps the driver's object loops honest: a
  `continue` is either a selection, or diagnosed, or declared.

### Verification and numerical uncertainty

- 45 declared cases (was 27), 17 checked-in `.k2d` benchmark inputs, 150
  automated tests.
- **A numerical uncertainty band for published numbers**: the Grid Convergence
  Index after Roache (1994) and Celik et al. (2008), in the ASME V&V 20 sense,
  with `mesh::refine_uniform` supplying the nested triplets it assumes. The
  estimator comes from the verification literature, and is itself verified
  against manufactured triplets whose answer is known.
- Measured with it: the Giroud rigid-footing benchmark converges to 15.244
  kN/m ± 0.21% — the file's own mesh reads 0.5% higher, and the difference is
  the mesh, not the physics.
- `K2D-A005`: a factor of safety computed with a non-associated flow rule
  depends on the mesh and falls as the mesh is refined (−7.9% over a fourfold
  refinement on Griffiths & Lane). The run says so.
- A verification case that measured nothing was found and fixed: the
  strength-reduction search hard-coded its own trial tolerance, so a sweep over
  the tolerance had been setting a number nothing read. Re-measured with the
  control connected, the factor of safety is bit-identical at and below 1e-3
  but +2.0% at 1e-2 and +45.6% at 1e-1 — always on the unsafe side.

### Compatibility

- Reads `.k2d` versions 1–12. Files written by this build are refused by 0.6.x,
  by design.
- The `.res` results format is unchanged (version 5).
- No breaking change to the Python surface or the CLI contract; both gained the
  new inputs.

## [0.6.1] — 2026-08-07

Patch release: the published binaries catch up with the tree.

### Engine
- Axisymmetric staged (Plastic) phases: load and prescribed-displacement
  staging now runs in r–z kinematics from every front end — the
  displacement-controlled bearing-capacity workflow works axisymmetrically.
  Previously any axisymmetric project was refused by the phase driver, so
  the capability was unreachable from the CLI and Python. What is still not
  plumbed for axisymmetry is refused honestly, never integrated as plane
  strain: activation changes (excavation / fill) and consolidation, flow,
  dynamic and safety phases.

### Verification
- Three new published benchmarks in the corpus, each a checked-in `.k2d`
  solved from the file by the suite:
  - `KV-FND-013` — bearing capacity of a smooth rigid circular footing,
    axisymmetric Mohr–Coulomb (Cox 1962 slip-line solution).
  - `KV-FND-014` — smooth strip footing on clay with strength increasing
    with depth (Davis & Booker 1973).
  - `KV-SLP-002` — Griffiths & Lane (1999) Example 1, the homogeneous 2:1
    slope, factor of safety by phi–c reduction against the published pair
    (their FE 1.4, Bishop & Morgenstern charts 1.380).
- A benchmark walkthrough document: the three cases end to end from the
  command line, with the published values beside the computed ones
  (`docs/validation/three-published-benchmarks.md`).

## [0.6.0] — 2026-07-20

The first published version: the engine, the command line, the Python
surface and the verification record.

### Engine
- Staged construction (plastic phases with activation inherited from the
  previous phase), K0 procedure and gravity initial stress.
- Constitutive models: Linear Elastic, Mohr–Coulomb with a tension cut-off,
  Hardening Soil, HS-small, Soft Soil, Soft Soil Creep; drained,
  undrained (A/B) and non-porous drainage types.
- Groundwater: steady-state and transient flow with free surfaces and
  seepage faces, Biot consolidation (elastoplastic), fully coupled
  flow–deformation.
- Dynamics: seismic time histories (harmonic, Ricker, stored accelerogram
  records) with compliant-base and free-field boundaries, Rayleigh damping,
  optional full plasticity during shaking, response spectra.
- Safety: phi–c strength reduction with an honest lower-bound flag.
- Structural elements: plates with elastoplastic Mp/Np hinges, embedded
  beam rows, node-to-node anchors, geogrids, Coulomb interfaces.
- Elements: 6- and 15-node triangles; plane strain and axisymmetric.
- Design codes: EC7 (EN 1997-1) and TBDY 2018 partial-factor catalogues.

### Front ends
- `katai` command line: `solve` / `validate` / `info` with a documented,
  test-pinned exit-code contract (0/2/3/4/5).
- Python package: the engineer-facing `katai.Project` builder over the
  published facade, plus the same `katai` command as a pip console script.
- The `.k2d` input format (version 2) with a JSON Schema, and the `.res`
  results format (version 5); one contract shared by every front end,
  pinned byte-for-byte by the suite.

### Verification
- The full CTest suite runs on every composition of the build; the
  verification corpus (13 published benchmark cases as checked-in `.k2d`
  files) and the generated verification matrix with its bibliography ship
  in the repository.
- The `portable` preset reproduces the entire suite without any
  proprietary component.
