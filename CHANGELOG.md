# Changelog

All notable changes to KATAI 2D. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
MAJOR.MINOR.PATCH.

## [Unreleased]

### A wall's joint keeps the stress it was installed with

A joint's normal stress is the stress it was installed with plus what its own opening and closing
add. The phase chain carried a joint's slip but rebuilt that installed stress in every phase from
the at-rest formula on the phase's own geometry, and the difference to the stress the phase before
had been in equilibrium with acted as a load between wall and ground: an excavation took
K0 gamma (depth dug) off the pit face at once with nothing beside it having moved, a wall started
by gravity loading took the at-rest pressure in its first phase, and a phase that changed nothing
after a wall installed on loaded ground moved the wall -- its joint's stress jumped by 48% of what
it carried. The installed stress is now carried from phase to phase (KV-STR-021). The same test
checks the seed of a wall in place from the start against the closed form: each face takes
K0 times the overburden of its own side, read a millimetre into it, which a wall in a step of the
ground surface makes visible (the low face carried the high side's 81 kPa at the toe instead of
45 before that reading).

### A Head edge above its water is a seepage face in the steady flow too

The steady groundwater flow held a Head edge at its head over its whole length, so the part of it
above the water it stands for drew the ground's water out by suction, through a face open to the
air. A dam drawn with its tailwater edge up to the crest lost its exit point: the water left at the
tailwater level (1.00 m) instead of 1.31 m up the face, and the head field differed by up to 1.8 m
from the same dam drawn with an explicit seepage face. Charny's discharge hid it (+1.7% against
+0.8% now): it holds whatever the face above the tailwater does. The part of a Head edge above its
water is now a seepage face, as in a transient flow phase and in the steady flow of a staged phase;
a model whose every Head edge stands above its water is refused with the reason (KV-FLW-010).

### A drawdown in a transient flow phase: seepage faces, wells, drains and a moving water level (`.k2d` v22)

A transient flow phase closed its seepage faces (with a warning), refused wells and drains, started
from the water table whatever the phase before it had done, and applied every head at its first
instant -- so the classical rapid drawdown of a slope could not be modelled. It now:

- resolves **seepage faces, drains and a well's h_min** as one-sided conditions in every time step,
  by the active-set rule the steady flow solve applies to the same objects (wells and drains are
  read through one shared resolution, so a well's discharge is spread along its line the same way
  in both). A face is open at zero pressure while water leaves the ground through it; a drain only
  takes water out; a well pumps its discharge until the ground falls to its h_min;
- makes the part of a **Head edge that its water stands below a seepage face**: held at the head,
  it drew the ground's water out by suction through a face open to the air;
- **starts from the head the phase before it ended with**, when that phase was a groundwater flow
  (the result carries the nodal head, suction included); otherwise from the water of the phase it
  continues;
- gives a **phase's own water level to its Head edges**: a drawdown is a phase whose level is lower;
- moves the Head edges **linearly over the duration** when asked (`phases[i].headramp`, `.k2d` v22):
  a reservoir lowered at a rate. An older build would make the change at the first instant;
- takes the **storage of a step against the state at the end of the step before**, prescribed nodes
  included, so a head that changes in time stores and releases the water it should (nothing
  changes for a constant head);
- **splits a step that does not converge** (up to 512 parts, each given four times the iterations),
  and **damps** the Picard update when it grows and caps it at a metre (measured: a 1142 m update
  where a drain asked for 6);
- is solved with the **sparse direct solver** every other phase uses. It was given a dense LU,
  refactorised every Picard iteration, which grows with the cube of the node count;
- is solved on the **ground active in the phase**: the soil an excavation had removed went on
  conducting and storing water. The **faces an excavation exposes** take the phase's water as their
  head (a seepage face above it), as the steady flow of a staged phase treats them -- a pit pumped
  down takes water from the ground around it;

A phase that takes its pore pressures from a groundwater flow (`flowpore`) right after a transient
flow phase now takes the **head that flow ended with**, not a steady flow of its own, and says so
(K2D-A038); it hands the head on, so a Safety phase after it reads the same field. A drawdown thus
reaches the effective stress and the factor of safety (KV-FLW-009: the settlement of a column whose
water table a drawdown lowered, against the closed form to 0.000%).

Verified (KV-FLW-008) against closed forms: a water boundary raised at a constant rate within 0.20%
of the rise of Carslaw and Jaeger's i2erfc solution (the change made at once lies 20% away); a
drawdown below the ground surface whose excess over the face dissipates as Terzaghi's isochrone
within 0.04%, the face at zero pressure; a well within and beyond what the ground can supply, at its
Darcy head and at h_min, exact; a drain that holds its head under the water and feeds no water into
drier ground; the same ramp through the file and the phase sequence within 0.78%; the ground under
an excavation draining as the same ground meshed alone; and a pit pumped down, run to its steady
state, within 1.3 mm (0.04% of the head drop) of the steady flow of a staged phase -- a separate
solver. K2D-A034 is retired. On the default sand retention curve the ground a drain dries four
metres above the water stopped the Picard iteration (a 6 m drop asked for a 1142 m update);
damped, capped at a metre an update and given four times the iterations once a step is split, it
converges (the chord capacity of the retention curve was tried and did not help), and a pit in the
same sand, run in 50-day steps that have to be split, reaches the steady staged flow within 0.8 mm.

### A wall along a pit's edge reads each face's own overburden

The initial normal stress of a wall's interface was read with the overburden on the wall's line
itself. A wall drawn along the edge of a pit's region has that edge under it, and the half-open
point-in-polygon rule put the line in one region or the other: after the excavation one face took
the other face's overburden, and which one depended on the side the pit was drawn on. Each face now
reads its own side, a millimetre into it. Found on the published staged excavation by bisection
(the shift of its wall followed this lookup alone): its wall now moves 23.58 / 56.73 mm in the
second and third stages (published 24.1 / 57.1; before, 23.46 / 56.47).

### Every refusal carries a code

Nineteen refusals stopped a run with a message and no code, so a script could only match the prose.
They now carry stable codes in two new families: `K2D-I…` for input a run cannot be built from (an
elastoplastic plate with no capacity, an element without a material, pore
pressures from a flow with no water, an impermeable plate the mesh is not split along, a phase
naming a material the project does not have) and `K2D-S…` for a combination this build does not
solve, refused by name (a prescribed displacement outside a Plastic phase; what axisymmetry does not
solve yet; flow pore pressures or region water in a phase that computes its own; a semi-permeable
line where the flow cannot take one; an interface activated outside a Plastic phase; a contraction
on a plate installed in the same phase). Each is triggered by a test, with the same model without
its condition (KV-DIA-003).

### The pore pressure shown is the one the phase computed

The pore-pressure field of a result showed the steady part only -- the water table or the flow the
phase was set up in -- and the total stresses were taken with it, so during a consolidation the
pressure an engineer read was the hydrostatic one, not the one still draining. It is now the steady
part plus the excess the phase generated, and the excess pore pressure is a field of its own.

### Every diagnostic code is checked to be said, and not said without its cause

Sixteen documented codes were raised by no test. Each now has the smallest model that raises it and
the same model without the condition, which must not (KV-DIA-003). The input-safety gate now fails
when a documented code is raised nowhere, a raised code is not documented, or no test triggers it;
three codes that only a numerical event or a flow field computed outside the run raises are listed
there with their reasons. On a structured mesh of right-angled 6-noded elements no point is farther
than a quarter of the element from every node, so K2D-G013 can never be raised there -- measured, and
the test is placed on an unstructured mesh.

### The phase editor reaches every phase key

A region's material in a phase, a plate's section and its line contraction, the reset of the
displacements, the load placed at the start of a consolidation, the integration tolerance, the
change of permeability `ck`, the end bearing of a wall's foot, the skin resistance at a pile's end,
the virtual thickness of an interface and the material of a free interface were in the file and the
engine but not in the editor. A new phase no longer copies what the phase before it did once (a
reset, a partial stage). A phase can be duplicated and moved up or down the chain, and removing
one asks first. A consolidation curve is drawn against time, resampled from its growing steps,
rather than against the step count, and a cross section reads any result field along a line,
interpolated with the elements' own shape functions (exact at the nodes), and copies it as CSV.

### A factor of safety with a Hoek-Brown rock

A Safety phase on a model with a Hoek-Brown material was refused: the criterion has no c' and
phi' to divide. The strength is now reduced inside the criterion (Benz, Schwab, Kauther &
Vermeer 2008): at each stress point's own sigma'_3 the envelope's tangent is a Mohr-Coulomb line,
reduced as Mohr-Coulomb is reduced (c and tan phi divided by the factor), and the reduced envelope
passes through it. The published form f/eta is implemented with its numerator rationalised, so it
does not lose digits towards the tensile apex, and the apex itself does not move under the
reduction. A user's tensile cut-off is divided like a Mohr-Coulomb one, and the dilatancy is held
within the reduced friction at zero confinement. The reduced envelope equals an independent
tangent construction to 7e-16 and a homogeneous rock block under a confining and a major stress
has the factor at which the reduced envelope passes through its state, within 0.1% (KV-CST-021).
The phase says that its factor is not that of a Mohr-Coulomb fit to the same envelope (K2D-A037).

### Axisymmetric models: excavation and fill, interfaces, and a factor of safety

An axisymmetric model refused every change of the active regions, on the ground that the
r-weighted assemblies took no element-activity mask; they did, and two calls did not pass it. A
shaft, a cylindrical pit or a tank foundation dug and backfilled now runs: the heave of a confined
column losing its top layer is the one-dimensional closed form (exact on 15-noded elements, 0.1%
on 6-noded), and the stress after a pit is dug is the stress of the same geometry loaded without
the pit, to 3e-12 kPa (KV-EXC-003).

Interfaces are built in axisymmetry and integrated per radian, each Newton-Cotes point weighing its
radius: the shaft and base of a pile modelled as a volume, the joint under a tank. A cylinder
pulled through a cohesive joint carries c_i a H per radian at full slip, to every printed digit,
and the elastic load-transfer closed form within 0.3% before it slips (KV-STR-016). Plates,
anchors, geogrids and embedded beams stay refused there: a plate in axisymmetry is a shell with a
hoop force, a different element.

A Safety phase runs in axisymmetry: every trial of the strength-reduction search is a nonlinear
solve with the phase's options, and the r-z kinematics now reach them all, with the parent's
internal force r-weighted. A circular footing on weightless Tresca clay at half its collapse
pressure has a factor of safety whose product with the pressure is the collapse pressure the same
model reaches by loading, within 0.7%, at 5.5 s_u -- the circular footing's, not the strip's 5.14
(KV-FND-017).

### A joint slips under a damped earthquake, and the ground at the model's right edge has weight

Two defects, both found by the first closed-form check of the nonlinear dynamic phase (Newmark's
sliding block, KV-DYN-004: the permanent slip of a block on a frictional joint under an
acceleration pulse):
- **The stiffness-proportional part of Rayleigh damping was built on the full initial stiffness,
  the joints' elastic springs included.** A joint spring is stiff by construction, so a joint
  slipping at a rate v met a viscous resistance beta k_s v far above its Coulomb strength: with 5%
  damping the block slipped 0.8% of the closed form. Every nonlinear dynamic run with damping and an
  interface (a wall beside its soil) had its joint slip suppressed in this way, in the unsafe
  direction. Contact springs -- interfaces and the embedded beams' skin and foot springs -- are now
  left out of the damping matrix, as is usual for zero-thickness contact elements.
- **The soil column above a node on the model's right or top boundary weighed nothing**, because a
  point on those edges is outside every region by the half-open point-in-polygon rule. An interface
  reaching the right boundary started its end station with no normal stress, and so no friction:
  8.6% too much slip on 6-noded elements, 3.9% on 15-noded ones. A boundary point now takes the
  region it bounds.

With both, the slip is Newmark's within 0.5% on both element orders. A pulse below the yield
acceleration leaves no slip.

### An undrained soil in a linear dynamic phase carries its pore fluid

The linear dynamic path assembled the skeleton's own elastic pair for an Undrained (A)/(B) soil,
so a saturated undrained clay shook as if it were drained. The shear wave was right; everything
volumetric was too soft, measured as a 9.9% larger peak displacement on a block standing on its
base. The undrained pair keeps the shear modulus and takes the undrained Poisson's ratio (the
static solver's own K' + Kw/n), and a depth profile scales with it. It now equals the nonlinear
path, which solves with the models themselves, to 1.7e-13.

### A nonlinear dynamic phase that stops early has failed, and reads its iteration controls

A time step that did not converge ended the phase with ok = true and "STOPPED early" in the
message only, so the phases after it ran as if the shaking had been completed. Such a phase is now
reported failed, with the time it reached as its load fraction and its fields up to then
viewable. The phase's own tolerated error and maximum iterations are read by the nonlinear dynamic
Newton; they were not, and a phase that set them was solved with 30 iterations and 1e-6.

### A prescribed inflow reaches a transient flow phase; a seepage face there is said to be closed

A prescribed-flux edge (a rainfall infiltration, a recharge) reached every time-dependent phase as
a closed edge, silently. A transient flow phase now applies it, integrated onto the nodes as the
steady flow solve does: a saturated column under an inflow reaches the closed-form head q L / k
above its base to every printed digit. A consolidation or fully coupled phase solves the excess pore
pressure, for which the edge is closed, and now says that the inflow is not applied (K2D-A033). A
seepage face in a transient flow phase was closed there and said so (K2D-A034); it is now solved
(see the drawdown entry above).

### A permeable wall lets water through in a transient flow phase

A wall with interfaces, or an interface on mesh edges, splits the mesh. The coupled phases tie the
two sides of a permeable seam to one pore equation; the transient flow solver was handed the split
mesh as it was, so every seam -- fully permeable included -- stopped the water (a 10.3 kPa jump
across a wall declared permeable). It now solves on a mesh in which each tied pair is one node. A
semi-permeable seam is refused there, as in the coupled phases.

### Wells and drains in a staged phase's groundwater flow are said, not dropped

The groundwater flow of a staged phase (`flowpore`) takes its boundaries from the water table; a
well or drain active in the phase was not part of it, silently. It is now said (K2D-A035), and a
well or drain active in a phase that takes its pore pressures from the water table, on which it
has no effect, is noted (K2D-A036).

### A reopened result says what the run said (`.res` v11)

The results file now stores the total displacement since the start of the staged construction
(the phase sequence's own sum, which leaves out the initial phase and the phases that commit no
state and starts again after a reset), whether a factor of safety is only a lower bound, the
design approach that makes it an over-design factor, the diagnostics, the excess pore pressure
field, the consolidation curve, the degree of saturation, the support reactions and the seismic
records. A reopened lower bound read as a definite factor, and a reopened run had lost the
warnings that qualify its numbers. Older files read as before, the new fields empty.

### Paths with letters outside ASCII; a save replaces the file whole or not at all

Paths are UTF-8 everywhere in the tree and are opened as UTF-8 on every platform. A narrow file
stream reads them in the ANSI code page on Windows, so a project under a folder or name with
Turkish letters could not be opened by the program that saved it. The command line and the Studio
run with the UTF-8 process code page. A project or results file is written beside its target and
moved over it once complete, so a save that fails part way leaves the file that was there.

### Deleting an object keeps every phase on the objects it was set for

Every phase holds what it does to the project's objects by index. Deleting a structure, a load, a
prescribed displacement or a well erased it from its list and left the phases alone, so every
later object took over its neighbour's settings in every phase -- the model solved was another
model than the one drawn. Deleting a region shifted its activation but not its dry flag, its water
or its phase material, and deleting a material left the layers, interfaces and phase changes that
named it on the wrong one. The kernel now provides the edits (`erase_struct`, `erase_load`,
`erase_disp`, `erase_hydro`, `erase_material`, `erase_polygon` for every phase setting,
`insert_polygons`, `split_edge`, `erase_vertex`) and `material_uses` to say what a deletion
changes; regenerating the ground from boreholes resets the region settings when the regions
change, and says so.

### New verification cases

KV-GEO-004 (an edit equals the project built without the object, in every phase setting),
KV-CON-005 (Soft Soil consolidating under a fill against an independent finite-difference column:
the degree of consolidation within 0.0011 at every step), KV-CON-006 (the change of permeability
against the same column, within 0.15%), KV-SLP-004 (initial stresses on a slope against the
infinite slope, within 0.03%), KV-DYN-004 (Newmark's sliding block, within 0.5%), KV-EXC-003 and
KV-STR-016 (axisymmetric excavation and interfaces), KV-FND-017 (an axisymmetric factor of safety), KV-CST-021 (Hoek-Brown strength reduction),
KV-FLW-008 (transient flow conditions), KV-FLW-009 (a drawdown's pore pressures reach the deformation:
the settlement against the closed form to 0.000%), KV-STR-019 (a tunnel lining that contracts in
elastic ground, against the thick cylinder and the ring's hoop equilibrium: displacement +0.47%,
hoop force -0.49%), KV-STR-020 (a plate that changes its section in a phase keeps the force it
carries and the new section takes what comes after: 5000 / 15000 / -5000 kN/m exactly),
KV-CST-022 (a soil region that changes its material keeps its stresses and the new material takes
the load after: the settlements of both phases exactly),
KV-DIA-003 (diagnostic codes); and declarations for tests that
had none: KV-STR-017 (the plate hinge's M-N diamond), KV-FLW-006 (transient saturated flow: the
isochrone, the erfc front, the steady limit), KV-FLW-007 (Darcy seepage), KV-STR-018 (Rankine's
thrusts), KV-FND-018 (the Prandtl-Reissner bearing capacity).

### A plate on the ground surface with an interface on each face is solved by either linear solver

Such a plate has no ground on its upper face. The mesh is split along it all the same, and the
nodes the split leaves on that side belong to no element. Nothing was done about them unless
the phase had excavated something: their joint stayed in place, and as soon as it opened under
load they had no stiffness -- a row of zeros in the system. The direct solver of the default
build passed over it by perturbing the pivot; the portable one refused it, as it should, and
the run stopped at a tenth of its load. A node that belongs to no element, that no plate moves
on, that is tied to nothing and carries no support is now fixed in every phase, and the joint
against it is switched off, as a joint beside excavated ground is. A node in no element that a
boundary condition holds is a support -- a joint drawn on a fixed boundary slides against
exactly such nodes -- and is left as it is.

The portable build is the one the continuous integration runs, and it failed there in three
tests. This was one; the test of an interface's effect was the second (see "Gravity loading
with a wall in place is gravity loading"); the third compared a restarted consolidation run
with an uninterrupted one to 1e-12, which is the round-off of one linear solver and not of the
other (1.3e-14 and 5.6e-12 measured), and now states the bound both meet.

### Hardening Soil is checked against the four results it has in closed form

Where one of the model's two mechanisms works alone there is a formula to compare with, and a
test now does (KV-CST-020): the drained triaxial curve is the hyperbola the model is built on
(five strains, to 0.001%); the steepness of the cap for m = 1 and no cohesion follows from the
two conditions of the oedometer (0.6602 calibrated against 0.6599); undrained elastic loading
keeps the mean effective stress and raises the deviator exponentially (four strains, to six
figures); and the undrained path meets the cap at the deviator the ellipse through the
pre-consolidation state gives (236.6 kPa, within 0.1). Where the two mechanisms work together
-- undrained shearing past the cap -- there is no formula, and that is where this model and
published curves of other programs differ by a few percent (KV-CST-019).

### Consolidation phases use every core, and build their matrix only to solve with it

Nearly all of the time of a consolidation phase with hardening soils is the integration of the
stresses, and that loop ran on one core; the static phases have used all of them for a long
time. It is now the same two-phase assembly: the element matrices are computed in parallel,
each element writing only its own block and its own stress points, and added to the system in
element order, so the sums do not depend on the number of threads.

The balance of an iterate needs the internal force and the coupling blocks, not the system
matrix; those blocks do not change inside a time step. Every trial of a line search used to
add the whole system to a sparse builder and discard it. The matrix is now built from the
element matrices of the last assembly, where it is solved with.

The thread pool hands a range out in small pieces to whichever thread is free. It used to cut
the range into one contiguous block per thread, and the work is far from uniform: the stress
points of freshly placed fill take two hundred substeps where their neighbours take one, such
elements are neighbours in the numbering, and they fell into one or two blocks -- two cores
busy of eight, measured. The static phases gain from this as well. Which thread takes which
piece depends on timing; the answer does not, because every index writes to its own storage.

Measured on fill placed on soft clay: 1.5 m on 15-noded elements 473 -> 129 s, 2 m 582 ->
151 s; 2 m on 6-noded elements 69 -> 25 s. The settlements moved within the tolerance of the
phase (192.3 -> 192.8 mm, 296.5 -> 293.3 mm on 15-noded elements, where 6-noded elements give
191.5 and 293.0), because the balance is summed in another order. A developer timer of the
phase is printed with KATAI_NL_PROFILE.

### A consolidation step is not made worse, and a cut step is kept within one percent

Two things kept fill placed on soft clay from finishing on 15-noded elements, where the same
model on 6-noded elements ran in a minute.

The search along a Newton step took the least bad of its lengths even when every one of them
made the balance worse. Freshly placed fill has stress points with no stress and so no
stiffness; the matrix they leave is close to singular, and one such step took the balance from
0.04 of its scale to 3.9e7. A step whose best length still doubles the balance is no longer
taken: the iteration stays where it is and the next method -- the elastic stiffness,
accelerated -- takes over.

The fill, with 1 kPa of cohesion, no dilatancy and a base that spreads under it, is at its
strength while it is placed, and there the balance of the first time steps comes down to about
1% of the forces in play on every method and at every size of step, and no further (the
non-associated flow rule has no tighter equilibrium on some meshes; Tschuchnigg, Schweiger and
Sloan 2015). The first attempt at a time step is still held to the tolerance of the phase.
Once it has been cut in two and stalls again, its best iterate is kept if it is within 1% --
the tolerated error such analyses are run with by default, and the limit an increment of a
Plastic phase that cannot be cut further is kept at. The step is counted, the remainder is
carried into the next one, and the warning (K2D-A019) now prints the largest remainder a step
was kept with instead of the limit.

Measured. Fill of 0.5 / 1.0 / 1.5 / 2.0 m on soft clay, 15-noded: 109 / 121 / 479 / 582 s with
settlements of 36.1 / 112.2 / 192.3 / 296.5 mm, against 33.8 / 114.9 / 191.6 / 293.0 mm on
6-noded elements; only the second of the four finished before. The 6-noded results are
unchanged to every printed digit. A staged embankment on soft clay with two consolidation
periods finishes on 15-noded elements for the first time: excess pore pressure under the
centre 26.0 / 0.30 / 31.4 kPa after the three stages and 52.3 days to 1 kPa, against 26.5 /
0.29 / 31.5 kPa and 52.0 days on 6-noded elements; 41 of its 222 time steps were kept this
way, the largest remainder 8.4e-3.
The factor of safety after its first lift is 1.064 on 15-noded elements, 1.4% under the
published 1.08; on 6-noded elements it is 1.129. The first lift takes 72 minutes on 15-noded
elements where it takes two on 6-noded ones, which is the next thing to bring down. Seventy-five
time steps per stage instead of twenty-five leave the time to 1 kPa at 52.0 days.

### The undrained triaxial comparison of HSsmall is made on the published curves themselves

The comparison with published element tests of four clays (KV-CST-019) took its reference
values off a figure whose strain axis ends at 3%, and compared them with the model at 3%. The
figure is a vector drawing: its curves run on to 10% axial strain, and the values that had been
read were the ends of the paths. The deviator at 3% and at 10% is now taken from the curve
coordinates and the mean effective stress at 10% from the end of the stress paths, and the
bands are narrower than the 10% they were: 8% and 5% on the deviator at 3% in compression and
in extension, 5% and 3% at 10%, 4% and 3% on the mean stress.

The model is low in the deviator: 5 to 7% in compression and 3.5 to 3.8% in extension at 3%,
3 to 4% and under 2% at 10%. That is open. The source prints two values of the normally
consolidated stress ratio for these clays and gives no pre-consolidation for the tests; its
curves, read for all seven of its clays, gain deviator faster after the cap is met than the
hardening law does with a stiffness that follows the falling confining stress, and its
shallow, more overconsolidated clays lose a quarter as much mean stress in shear as this
model. Forty-eight combinations of the four choices the documentation leaves open or states
twice were measured against the fourteen curves; none reproduces them, and none was adopted.
The entry in the verification matrix carries the numbers.

### A pile row's springs are those of the soil around each of them

The springs between a pile row and the ground are multiples of the shear modulus of the
surrounding soil. One modulus was used for the whole pile, that of the soil at its middle: a pile
through soft clay with its foot in sand stood on a foot spring made of the clay's modulus, and
its shaft in a clay that stiffens with depth had the springs of the clay at mid-depth all the
way. Each skin point now takes the modulus of the soil it lies in (with its increment with
depth, where the material has one) and the foot that of the soil it stands in.

Measured on a building on two pile rows beside a tunnel (clay with E' = 3400 kPa + 400 per
metre, pile feet in a sand with E50 = 25 MPa; published: the near corner settles 19 mm, the far
one 2 mm, tilt 0.168%): 17.7 and 5.4 mm, tilt 0.123%, before; 18.8 and 2.9 mm, tilt 0.159%,
now. The far pile row could not hold its end of the building up. In a two-layer test the
slope of the pile's axial force changes at the boundary by the ratio of the moduli (10.02 for a
ratio of 10; it was 1.19).

The foot of a pile with a base resistance also has its elastic zone -- the half disc of the
pile's equivalent radius below it -- where the mesh is fine enough to have a stress point there.

### Resetting the displacements resets the small-strain history

A phase that starts the displacement count again is a new event read on its own, and the strain
the earlier phases left in the soil is not the history it meets: the small-strain history of
Hardening Soil with small-strain stiffness is cleared with it, as "reset small strain" clears
it. Only the read-out was reset before. The phase says so (K2D-M005).

### An end-bearing zone the size of the wall's foot

The zone of excluded plasticity at the lower end of a plate with end bearing was every element
with its centroid within the plate's equivalent thickness d of the foot: nothing on a mesh
coarser than d, and whole elements -- the mesh's size, not the wall's -- on a finer one.
Measured on an anchored excavation with d = 0.35 m and a published factor of safety of 1.49,
on a 1 m mesh of 6-noded elements: 1.71 with the option and 1.55 without. The zone is now the
foot of a wall of thickness d: the stress points in the half disc of radius d / 2 beyond the
lower end, the ground the foot bears on -- the same construction as under the foot of a pile.

The zones are not kept in a Safety phase. A factor of safety is found by reducing every
strength in the model, and a zone whose plasticity is excluded has a strength that cannot be
reduced: on the same excavation with 15-noded elements the factor was 1.4875 with no stress
point in the zones and 1.536 with four at the feet of the walls and seven at the ends of the
grout bodies. The phase says that it left them out (K2D-A029). A beam with a free connection
-- a grout body, a nail -- is not a pile standing on its foot and has no zone.

### A factor of safety beside a wall says which elements it was computed on

On the same excavation, 6-noded elements give 1.63 on a 2 m mesh and 1.55 on a 1 m mesh where
15-noded elements give 1.4875 on the 2 m mesh (published 1.49). A Safety phase on 6-noded
elements with a wall or an interface in the model now says so (K2D-A031).

### A stress point in tension no longer returns NaN from the small-strain flow rule

The mobilised friction, q / (sigma_1 + sigma_3 + 2 c cot phi), was not bounded. At a point in
tension beside a strip load, a fraction of a kPa from the apex, it reached 5; the contractant
branch of the small-strain dilatancy rule then took the logarithm of a negative number, the
stress came back NaN, the tension cut-off replaced it with its own value, and only the tangent
-- and a linear solve that returned NaN -- showed it. Increments were kept at their best iterate
with a warning, and the run went on. The mobilised friction is bounded by phi, which is what it
is at failure. Plain Hardening Soil reads the same quantity and never produced the NaN; past
the failure line its mobilised dilatancy is now psi where it could fall back to zero.

### Freshly placed Hardening Soil fill

A stress-free point that has never been loaded stands at the tip of its cap, where the cap's
gradient is zero: an explicit step from there is elastic, the same step read at its own end
is plastic, and the two disagree by the same proportion at every step size. Every integration
from the tip ran to its ceiling of substeps and was reported as one that had not met its
tolerance; an increment of ten kPa from half a kPa did the same. Three cases are now told
apart: pulled, a point on the tip stays there (the cap hardens only by compaction); pushed, it
takes the uniform walk at its full count, which is what finds the self-similar path of primary
loading, without a measurement that cannot succeed; and a point near the tip asked for an
increment more than four times its own size is walked in substeps a fiftieth of what the state
then amounts to. A consolidation step whose iteration moves without getting anywhere -- the
best of eight iterations not a tenth below the best of the eight before -- is treated as
stalled, where it used to spend its whole budget on steps of a hundredth of the Newton
direction.

### The elastic-stiffness pass of a consolidation step takes its full step, accelerated

A time step whose Newton iteration has stalled goes on with the elastic stiffness, as a Plastic
phase does. There it had to lower the balance at every step, which a fixed-point iteration on
a response that is not smooth does not do: it took the least bad of fourteen shrinking steps
and stayed where it was (measured on fill placed on soft clay: the balance held at 2.3e-2 of its
scale with steps of a five-thousandth). It now takes its full step unless the balance rises
tenfold, as in a Plastic phase, and the step is Anderson-accelerated from the last six
(Anderson 1965; Walker and Ni 2011): the same placement on 6-noded elements in 71 s where it
took 155, with the same settlement. A developer check of the consolidation matrix against the
balance is there as well (KATAI_CONS_TANGENT_CHECK).

On 15-noded elements that was not enough; see "A consolidation step is not made worse, and a
cut step is kept within one percent".

### One cap calibration per soil

The cap of a Hardening Soil material is calibrated by walking an oedometer path some three
thousand times, and the material models are built for every phase: on a tunnel analysis of six
phases with two such soils the same two calibrations were done 37 times and took 45% of all the
constitutive work of the run. The result is a function of the parameters alone and is kept.

A stress point held at the tension corner that is pulled further is returned there directly
instead of by a walk that ends there, and the numerical tangent of Hardening Soil is differenced
with a step below the increment (a fixed 1.5e-8 was as large as the smallest increments).

### HSsmall phases that reverse the ground converge

In this model the stiffness of a step is set by the direction of the step's strain increment
against the strain history, so the stress is a function of the increment with slopes of order
1 / |d eps| wherever the increment is small beside the history. Three things followed, and a
tunnel lining installed after grouting -- a phase that reverses the ground everywhere -- stopped
at 1 to 8% of its load.

The tangent did not contain the slope of the stiffness, and where it was differenced the step
was a fixed 1.5e-8, as large as the increments it differenced. It now carries that slope, taken
with a step that is a fraction of the increment: at 3000 random stress points the tangent
differs from the difference of the stress by more than 10% at 23 (it was 748). An iteration
starts with G0 where there is no increment yet to give a direction, the stiffer of the two
answers a point can give.

A threshold below which the stiffness was simply G0 put a step of 0.2% in it; it is gone (the
secant goes to G0 by itself), and with it three quarters of the steps in the stress that a scan
of straight lines through increment space finds.

An increment the iteration still cannot close is solved with the stiffness read from a fixed
strain increment -- the best iterate it reached -- which is a smooth problem; the fixed
increment is then replaced by the answer and the increment solved again, twice. The phase says
that it did this and how far the two were apart at the end (K2D-A032). The tunnel analysis runs
in 3 minutes without a single cut-back where it was refused after 13. A consolidation phase has
the same device per time step: a strip load taken off a sand with small-strain stiffness and
consolidated to rest ends where the drained unloading does (2.540 against 2.539 mm heave).

### HSsmall remembers its strain history as a tensor, and hardens as its formulation says

Two parts of the small-strain model were not the model's own (Benz 2007).

The strain history was one number that only grew. A point that had been loaded and was then
unloaded stayed on the stiffness the loading had degraded it to: after a surcharge placed and
removed, an excavation heaved 3.20 mm, the plain Hardening Soil answer to 0.03%. The history is
now the deviatoric strain tensor of the formulation. Each strain increment is looked at in its
own principal directions; where it runs against the history, the history is reset there, and
the stiffness of the increment is the secant of the stiffness curve over the increment itself.
The same excavation heaves 1.005 mm against 1.018 mm for the unloading law continued from the
surcharge's relief (-1.2%), and a first unloading from rest now meets its closed form to five
figures (7.132240e-4 m against 7.132244e-4; it was 0.30% short with the tangent at the start
of the increment). "Reset small strain" does what it did.

The hardening laws did not carry the small-strain multiplier. While the small-strain stiffness
is active the plastic strains of primary loading must be held back in proportion,
d gamma_p = h_i d lambda with h_i = G_m^(1 + E_ur/E_i), or plasticity takes the response over
from the first increment and the model shows none of its small-strain stiffness in primary
loading. With it, the undrained triaxial compression of four published clay sets is at -7.2 to
-7.7% of the published deviator at 3% axial strain (read off a figure; it was -6.6 to -8.7%),
with the mean effective stress within 1.2%, and within 3% of the deviator at failure. Analyses with HSsmall materials become stiffer at small strains: walls and footings
move less in their first stages.

The small-strain multiplier belongs to the ELASTIC stiffness and to nothing else. It was applied
to E_ur itself, which is also a parameter of the shear yield surface, so the surface moved as
the multiplier fell: every increment of a monotonic loading began inside its own yield surface,
answered with the small-strain stiffness until it found the surface again, and took the
integration to its ceiling of substeps doing so. On a strip load on sand: 13 substeps per
integration against 3 for plain Hardening Soil, 35 s against 6; now 2.6 substeps and 8 s, and
a settlement of 5.48 mm where it was 5.21.

Plain Hardening Soil is not touched by any of these.

### Gravity loading with a wall in place is gravity loading

A model started by gravity loading was run by the K0 procedure instead whenever a wall with
interfaces, or a bare interface, was active in the initial phase -- the joints' starting stress
was the K0 procedure's, and the whole start was switched to match. Nothing said so. The ground
did not settle, and its lateral stress was K0 sigma_v where elastic ground under its own weight
has nu / (1 - nu) sigma_v: 0.50 against 0.25 at nu = 0.2. Gravity loading now starts stress-free
with its joints stress-free as well: 14.56 mm of settlement with the wall against 14.58 without
it, the lateral stress ratio 0.246 away from the wall. A wall left out of the initial phase and
installed in the next one, the usual way, was never affected.

A wall installed in ground that carries an excess pore pressure takes over the total stress at
its line. Its joints were seeded with the effective stress alone, and since friction acts on
the transmitted stress less that pressure, the pressure came off twice: under a 50 kPa
undrained surcharge the joint started at 18.7 kPa where the ground pressed on its line with
67.0 (18.7 effective, 48.3 in the water).

The joints' own results -- shear, normal stress, slip, and whether the Coulomb limit is reached,
station by station -- are readable from Python (`result.interface_forces`); they were computed
and shown in the application only.

The test that an interface changes the result had been passing on that defect: it compared a
model with an interface, started by one procedure, with one without, started by the other, and
read the difference as the interface's. Both are now started alike and the interface is read
where it acts: beside a weak joint under a surcharge the ground on the unloaded side settles
1.9 mm where it is dragged down 39.5 mm when bonded, the loaded side 112.9 mm against 74.8,
and the largest bending moment of a wall changes by 10 to 24% with joints on its faces. The
largest settlement of those models is away from the joint and moves by 0.4%.

### With walls and anchors, the factor of safety is the plateau of its path

A Safety phase on ground that carries structures reduced the strength in prescribed steps and
stopped where the iteration gave a step up. That is not a property of the ground. On a
dewatered excavation between anchored walls (6-noded, 2 m) it reported 1.555 at 120 mm of
movement; with the faster solver of this release, which closes steps the old one abandoned,
1.615 at 700 mm -- the same model, the same start, every accepted step converged. Followed with
the factor as an unknown (the arc-length search already used for ground without structures)
the path is one curve with a plateau: 1.53 at 95 mm, 1.585 at 235 mm, 1.62 at 800 mm, level at
1.630. Both earlier answers were points on the way up. The search with structures now follows
the path to its plateau as well, in about the time the prescribed steps took; the prescribed
search remains for a path that cannot be started. Factors of models with structures can rise by
a few percent with this. The closed forms of a block on a joint continued from a Plastic phase
-- pushed sideways, held by a prestressed anchor, the anchor removed -- come out at 1.501906,
2.002542 and 1.501907 against 1.501907, 2.002542 and 1.501907.

### A support beside undrained ground reports the water's share of its force

The reaction at a support, and the force an imposed displacement needs, were taken from the
effective stress alone. Beside drained ground that is the whole of it. Beside undrained ground
the excess pore pressure carries most of a load, and none of that was in the number: a strip
footing pushed 5 cm into undrained clay (s_u = 20 kPa, weightless) reported a footing force of
22.7 kN/m, where the same clay declared drained reported 98.6 and the closed form is 102.8; the
base of a confined column loaded undrained with 200 kN/m reported 6.9 of it. Nothing said so.
The reaction is now taken from the stress the solver balances, effective stress and excess pore
pressure together: 104.7 kN/m for the footing, 200.0 for the column. Displacements, stresses
and factors of safety were never affected -- only this reported force.

The pressure itself is reported as well. A Plastic phase on undrained ground generates an
excess pore pressure and every later phase is built on it, but its result carried the steady
pressure only; the field of the excess pressure was filled by consolidation phases alone. It
is now there after any phase whose ground carries one (`excess_pore`, node by node): 48.27 kPa
in the water and 1.73 in the skeleton under that column's 50.

Still not in a reaction, and now written at the field: the steady pore pressure (under the
water table a base reads the buoyant weight above it), the share of a load applied at the
supported node itself, and the end force of a structure or a joint that stands on a support
(K2D-A004, which now also speaks for a joint drawn on a fixed boundary).

### A joint reads the pore pressure of the ground on either side of it

Friction in a joint acts on the effective normal stress, the transmitted stress less the excess
pore pressure of the ground beside it. The ground was looked for on one side only -- the side
the joint's own nodes are on -- so a joint drawn on a fixed boundary, whose own nodes are the
boundary's, found none and slipped on the total stress. It now reads the ground on whichever
side has any, and the mean of the two between two soil bodies. A block on a rough base, held
between rigid sides, loaded with 50 kPa and pushed sideways: 36.6 kN/m of friction with no
load, 100.2 loaded drained, 38.6 loaded undrained (closed forms 37.4, 102.4, and 37.4 plus the
little of the load the water does not take); friction on the transmitted stress gave the
undrained block the drained one's 100. This is the first test of that rule
(python/tests/undrained_reactions.py).

### Soft Soil returns to the cap it is returned to, and converges in two dimensions

A strip pressed into Soft Soil clay took 963 iterations where the same strip in Mohr-Coulomb
took 119, and the out-of-balance force would not come under a few thousandths of the forces in
play. The cause was in the stress update, not in the iteration. Where two principal stresses
are equal -- every point of a ground at rest, whose two lateral stresses are -- a compression
increment is returned to the cap's surface or to its edge, and on the edge the pair is brought
to a common value after the return. The cap was solved with the pair still unequal: its
deviatoric measure weighs the two differently, so making them equal afterwards moved the stress
off the cap, 0.11 kPa inside it at a mean stress of 40 kPa. The returned stress jumped by that
wherever an increment passed from one return to the other, and the numerical tangent, taken
across such jumps, was not the derivative of the update at any point checked (10% to thirty
times off). The cap, the Mohr-Coulomb surface and the creep rate of Soft Soil Creep are now read
at the state that is kept.

The same strip: 170 iterations and 1.9 s against 963 and 26 s; at a tolerance of 1e-5, 234
against 2786; on 15-noded elements 270 against 2264; on a 0.25 m mesh 252 against 1846. The
footing force is 65.46 kN/m against 65.30 before (+0.24%), the earlier value being that of
points left inside the cap. A test scans the update across the switch for both models
(test_soft_soil, h).

### A footing pushed into sand finishes on a fine mesh, and every increment starts from the last

A rigid circular footing pushed 5 cm into sand (Mohr-Coulomb, phi' = 30, psi = 0, 15-noded
elements) needed 2572 iterations and four minutes on a 0.25 m mesh and did not finish in fifty
minutes on a 0.15 m one. Three things in the solver, none of which changes the equilibrium an
increment has to reach:

- Each load increment is started from the displacement of the one before it, scaled to its size,
  where that leaves a smaller out-of-balance force than the start it had. The old start read
  every point lying on its yield surface as elastic, so each increment began in ground eight times
  stiffer than it was and spent ten of its fifteen iterations finding the plastic zone again.
- The line search, having found a step that reduces the out-of-balance force, tries up to three
  longer ones between it and the step it turned down, and keeps the best.
- An attempt at an increment whose best iterate has not improved by a tenth in thirty iterations
  is ended there instead of at its budget of two hundred.

The same footing: 925 iterations and 67 s at 0.25 m, 14 minutes at 0.15 m. With associated flow
(psi = phi') 158 iterations against 395, and a strip footing 189 against 1615, the footing force
unchanged to every printed digit. Under a load the first of these does the same work: a strutted
excavation in Hardening Soil 624 iterations against 920, an embankment built undrained on soft
clay 1539 against 1897 with the factor of safety at 1.4010 against 1.4015.

A phase that imposes a displacement no longer ends short of it. Under a load a phase can stop
because the ground cannot carry more; a prescribed movement has no such limit. An increment of
it that can be made no smaller and still does not close is now kept at its best iterate when
that is within 1% of the load scale, and the phase says how many were and by how much
(K2D-A019). A strip footing on the same sand, which the faster start had sent down a path that
stopped at 59% of the settlement, reaches all of it at 173.9 kN/m against 174.4 before; the
largest force left out of balance in any of these runs was 0.19%.

What is not solved: with psi = 0 the out-of-balance force still settles between 1e-4 and 2e-3 of
the load scale in some increments, and the footing force then depends on the number of load
steps by up to 0.9% (595.7 kN in 20 steps, 590.2 in 80, on a 0.2 m mesh).

### A factor of safety above 3 is computed, not capped

The strength-reduction search stopped at a factor of 3 and reported anything that stood there as
"at least 3". A cohesive road embankment 3 m high stands at 3.2 to 3.8, so it had no factor of
safety at all. The search now goes on until the ground fails, up to a ceiling of 20: the
incremental search walks on, and the bisection used with structures from the unstressed state
widens its bracket, doubling, only when nothing has collapsed below 3 -- every factor below 3 is
found exactly as before. The 3 m embankment: 3.224, identical to the same slope run with c' and
tan(phi') halved and the factor doubled, as strength reduction requires. Only a model with no
mechanism at all reaches the ceiling, and that is still reported as a lower bound.

A Safety phase with a cohesive soil whose tension cut-off is on now says what that costs against a
slip-circle method (K2D-A026): the search counts no tensile strength in the zone where the soil
would be in tension, about 2c'/(gamma sqrt(Ka)) deep, where a limit-equilibrium method takes c'
along the whole slip surface. On those embankments it was about 2 of the 5-7% by which the factor
lay below Bishop's, the non-associated flow rule (K2D-A005) the rest; with both switched to the
limit-equilibrium assumptions the factor converges under mesh refinement to within 3-4% of it.

### A wall drawn in pieces is the wall drawn in one

A wall with interfaces drawn in two or more pieces -- cut at a layer boundary or at the excavation
level, or the chords of a lining -- gave a different answer from the same wall drawn in one, and
said nothing. Each piece was built as a wall of its own, and a wall hangs from its toe, the one
node it shares with the ground: the upper piece took the point where the pieces meet as its toe,
so the wall was tied to the soil on one side at every point it had been drawn through. On an 8 m
cantilever wall in front of a 4 m excavation the largest moment was 20.5 kNm/m in one piece and
18.1 in two, the head movement 11.8 mm and 9.7. A piece that starts where another ends now takes
that point as what it is, a point of the wall with an interface on each side. One piece and two:
20.52 and 20.53 kNm/m, 11.81 and 11.80 mm, on the same mesh.

### The total displacement, and a phase that starts the count again

A phase reports what it did; the settlement of a staged construction had to be added up by hand.
Every result now also carries the total since the start of the staged construction
(`displacement_total` in a script), the initial phase counting for nothing: it generates the
stresses the ground starts with. A phase can reset displacements to zero (`reset_displacements`,
`resetdisp` in the file), and the total is then counted from its start -- the settlement an
excavation or a tunnel causes to a building that stood there before. The state is not touched.

### A plate can be given another section in a phase

`plates={plate: dict(EA=..., EI=..., w=..., nu=...)}` on a phase (key `smat`): from that phase on
the plate has the new stiffness and weight -- a lining that replaces a shield, a wall that is
thickened. The plate keeps the normal force, the shear force and the bending moment it carries.
A plate is evaluated on its total deformation, so replacing its stiffness alone would multiply
those forces by the ratio of the two; the part of the deformation the new section does not
answer for is moved into the state the element keeps at its integration points, which now holds
a shear strain offset beside the axial and bending ones.

Checked on a plate on a linear-elastic block loaded in two steps with the section changed in
between (eight times the bending stiffness): in the phase of the change nothing moves and the
moments are those of the phase before to 2e-14; under the second load the moments and normal
forces are the old section's under the first load plus the new section's under the second, to
4e-13.

### An interface beside undrained ground slips on the effective normal stress

The normal stress an interface transmits is a total one: beside ground that is loaded undrained
it carries that ground's excess pore pressure with it. Its friction was computed from that
total stress. It now acts on the effective normal stress -- the transmitted stress less the
excess pore pressure of the soil element beside the joint, read at the start of the increment.

The difference is largest where the ground beside a wall is unloaded undrained: its pore water
goes into suction, and the suction was counted as a loss of contact stress, so the joint slipped
at a strength the ground beside it did not have. Measured on a strutted excavation in undrained
clay against a published solution: the wall's movement after the last stage was 162.6 mm
against 141.6 (+14.8%) and is 151.2 (+6.8%); the earlier stages move by 1% (23.5 and 56.5 mm
against 24.1 and 57.1). Drained ground has no excess pore pressure and is unchanged. The
opening of a gap is still decided on the transmitted stress, and the joint's own reported
stresses are the transmitted ones.

### Soft Soil returns a stress that is continuous in the strain increment

The Soft Soil return cut a strain increment into a whole number of equal substeps, the number
growing with the increment. An increment a hair larger than the last was therefore integrated
with one substep more, and the stress it returned jumped by the difference in integration
error. Inside an equilibrium iteration, where every iterate changes the increment a little,
that jump is noise in the internal force that no iteration can get under.

Measured on a footing pushed 0.15 m into this soil: the out-of-balance force stopped at
3-4e-3 of the forces in play, in pairs of opposite sign on neighbouring nodes all over the
mesh; the phase took 1 300 iterations at a tolerance of 1e-2 and was abandoned at 70% of its
settlement at 3e-3. The substeps now have the size the rule allows and the last one takes the
remainder, so where the count goes up by one the new substep has no length. The same footing
completes at 3e-3, at 1e-3 and at 1e-5. A drained column's settlement moves from 682.93 to
683.07 mm.

Soft Soil Creep, whose substeps are counted from the strain and the creep of the step, is cut
the same way, strain and time alike.

### The hardening soils are solved to 3e-3, not 1e-2

The tolerated error of a phase with Hardening Soil or Soft Soil in it was 1e-2, and at that
level the answer moved with the numerics by more than it is meant to be trusted to. One
retaining wall drawn in one piece and in two -- the same model, numbered differently -- gave a
largest moment of 51.56 and 50.12 kNm/m (2.8% apart); at 3e-3 it gives 51.42 and 51.52 (0.19%).
A strutted excavation's first-stage wall movement went from 6.01 mm to 6.19 against a published
6.1. The default is now 3e-3. It costs time where the iteration is slowest: an embankment built
undrained on soft clay takes 260 s instead of 160 for the same factor of safety (1.4012 against
1.4016); at 1e-3 it would take 657 s, which is why the default stops at 3e-3. `tolerance=` on a
phase still overrides it.

Not every increment can be brought to 3e-3: fresh fill at next to no stress leaves the
out-of-balance force wandering at 3-5e-3 of the forces in play. At the default tolerance such an
increment -- one that has stopped descending within 1e-2 -- is kept as it was before the default
was tightened, counted, and reported with the residual it was kept at (K2D-A019). With a
tolerance set on the phase nothing is kept above it.

### The initial shear hardening of Hardening Soil is that of the present stress

The pre-consolidation state an OCR or a POP defines is what the cap remembers; the shear
hardening a point starts with is the one that puts its current stress on the shear surface. It
used to be the larger of that and the shear hardening of the pre-consolidation state, so ground
with even a little overconsolidation began with a shear surface out of reach and answered its
first unloading stiffer than it should.

Measured on a strutted excavation in clay with a pre-overburden of 5 kPa, against a published
solution: the wall's movement after the first three stages was 6.01 / 22.61 / 55.50 mm against
6.1 / 24.1 / 57.1 (-1.5 / -6.2 / -2.8 %) and is now 6.01 / 23.68 / 57.25 (-1.5 / -1.7 / +0.3 %).
Normally consolidated ground is unchanged: its two states are the same.

The substep estimate of the same model measures its error against the stress, and no longer
against less than p_ref / 20 (it was p_ref / 200): at a point that carries next to nothing the
relative tolerance asked for the ceiling of substeps to resolve a thousandth of a kPa. A fill
on a soft base runs in two thirds of the time with the same answer.

### End bearing for a plate

`prevent_punching=True` on a plate (key `punch`). A plate has no thickness, so a wall loaded
down its own axis meets the ground under it only at a node and can punch through it. With the
option, soil plasticity is excluded in the ground within the plate's equivalent thickness
d = sqrt(12 EI / EA) of its lower end: the points there
answer with the elastic part of their own model (for Hardening Soil the unloading stiffness at
their stress) and their hardening state does not move. The zone is set phase by phase from the
plates active in it, and the phase says so (K2D-A029); a mesh with no element within d of the
foot has no zone, which is said too (K2D-A030). It is not the wall's tip resistance, and
it is not meant for sheet piles, whose tip is far thinner than d. The compression-index models
are not given an elastic branch and are integrated as usual inside the zone.

Checked on a 0.8 m wall pushed down by 400 kN/m: in elastic ground the option changes nothing
(32.251690 mm both ways); in a clay that yields the head settles 38.9 mm with it and 237 mm
without.

### The skin resistance of an embedded beam can vary along it

`Tskin_end` on a pile (key `Tskin_end`): the skin resistance then varies linearly along the beam
from `Tskin_max` at the point it is drawn from to `Tskin_end` at the point it is drawn to.
Checked on the pile-capacity fixture: with T at the head and 3T at the foot the pile carries
(2 T L + F_max)/L_spacing, the mean over its length, within the fixture's 2%.

### A structure that reaches a supported boundary is held there

Under the standard fixities -- no boundary conditions drawn -- a plate with a node on the model
boundary, where the boundary fixes at least one displacement direction, has its rotation fixed
at that node: at a rigid base that is the base, on a side it is the symmetry. On a drawn edge
the rotation is fixed where the edge holds only the direction normal to itself, which is a plane
of symmetry; any other drawn edge stays what was drawn, and a fully fixed edge under the end of
a beam is a pin. And a wall that an interface separates from the ground has displacements of its
own on the boundary, which take the fixity of the boundary they stand on.

Neither was done. A plate ending on the axis of a half model could rotate freely, so its bending
moment there was zero where symmetry makes it largest. And a lining ring separated from the
ground by its interface was not held on the axis at all: only friction kept it from drifting
sideways. Measured on a shield tunnel in a half model: with the interface the shield phase
stopped at 56% of its load with 178 mm of movement, and now completes (2.5 mm); without the
interface the shield's moments were one spike of -176 kNm/m at the springline and are now the
ovalisation they should be, +126 at the invert and -119 at the springline.

### Line contraction of a plate

`contraction={plate: 0.005}` on a phase (key `scontr`): the area loss of the closed line the
plate is part of, as a fraction, total at the end of the phase. The plate is given half of it as
a shrinkage strain -- the area of a closed line goes with the square of its length -- applied
with the phase's staged change. By itself it changes no stress; the ground around the plate
decides how much of it is realised. A later phase keeps the value, and 0 takes it back.

Checked on a plate along ground a million times softer than itself: with one end free it
shrinks by C/2 x L (50.0000 mm against 50.0000) and carries no force; held at both ends it does
not move and is in tension by EA x C/2 (5000.00 kN/m against 5000.00); taking the contraction
back undoes both.

### A consolidation phase is solved to the tolerated error of the soils in it

The equilibrium iteration of a consolidation phase was held to 1e-6 whatever the soil. A
hardening soil, whose iteration converges linearly and which every other phase solves to 1e-2,
was asked for a balance ten thousand times tighter there: it spent its forty iterations, the time
step was cut in two, and again, down to a sixty-fourth. A consolidation phase now uses a third
of the tolerated error a Plastic phase of the same soils uses (1e-3 for the hardening models),
and never a tighter one than before. The level is measured: on a Soft Soil column loaded
undrained and consolidated to rest, 1e-3 is within 0.1% of the answer at 1e-6 in a twentieth of
the time, and 1e-2 is not (7% short).

Three things in the same iteration were found on the way and changed with it:

- The force a time step leaves out of balance is carried into the next step. It was dropped, so
  a phase of some hundreds of steps lost it some hundreds of times -- nothing at 1e-6, half the
  settlement of that column at 1e-2.
- The line search shortens the step down to 2^-13 and, when no length improves the balance,
  takes the least bad one. It stopped at 1/64 and took that step whatever it did; on a
  compression-index model, whose stiffness is exponential in the volumetric strain, the balance
  of a dissipation step went from 0.16 of its scale to 1e+17.
- The volumetric strain is accumulated in a consolidation phase as in every other, so the void
  ratio the dilatancy cut-off reads no longer stays where the phase before left it.
- A balance that is not a number is never an improvement. The line search took a trial step
  as "not worse" unless it compared greater, and nothing compares greater than a NaN; one
  overflowed trial was accepted and every iterate after it was NaN.
- A step whose tangent has stopped giving a direction of descent is continued with the elastic
  stiffness -- the initial-stiffness iteration a Plastic phase falls back on in the same
  situation -- before it is judged. It used to be cut in two with the balance where it was,
  down to a sixty-fourth, and then given up. With the elastic stiffness a full step is a short
  one, so it is lengthened -- doubled while the balance keeps improving, up to sixteen times;
  without that the balance fell by 0.3% an iteration.
- A linear solve that is refused, or returns something that is not a number, is retried once
  with a fraction of the elastic operator on the tangent of the plastic points -- what a Plastic
  phase has done since it met a tangent made singular by points at a vertex. Only the tangent
  changes; the balance a step is driven to is the same.
- An iteration that has stopped moving is stopped: six iterations in a row at the shortest step
  length with no gain end it, where it used to spend the whole budget and then the same again
  on each half of the step (2 000 iterations on one time step, measured). For the soils solved
  to 1e-3 such an iterate is kept when it is within 3e-3 of the forces in play -- 0.4% from the
  answer on the column above -- counted and reported (K2D-A019), which now names the level
  applied.

Measured on a staged embankment over peat and clay (Soft Soil, undrained, with a Hardening Soil
fill): the first lift did not complete one time step in 50 minutes; the four phases now run
through, and with the change of permeability below the excess pore pressure under the centre
is within 2% of a published solution after each lift and the time to dissipate it within 7%.

### The permeability can change with the void ratio

`ck` on a material (key `ck`): in a consolidation phase log10(k/k0) = de/ck, with de the change
in void ratio since the initial state. The void ratio is read at the start of each time step.
Checked on a column compressed to a known void ratio and then consolidated: its settlement
after a time t is that of the constant-permeability column after t x k/k0, to 0.07% of the
final settlement. Entered where every soil is linear elastic it is not applied, and that is
reported (K2D-M007).

### Three small things a model said one way and solved another

- A stiffness or a cohesion that increases with depth is the reference value above the
  reference level and grows below it. It was extrapolated upwards as well, so a layer whose
  reference level lay inside it was softer and weaker above that level than the value entered.
- A plate bends against EI / (1 - nu^2): it goes on in the out-of-plane direction and cannot
  strain that way. With the nu = 0.15 of a concrete wall that is 2.3% stiffer than the EI
  entered; a plate with nu = 0 is unchanged.
- An increment with depth entered on a model that does not read it is no longer dropped in
  silence (K2D-M006).

### A hard increment is finished, or made smaller, before work is spent on it

Two changes to how a load-driven phase of a hardening soil handles an increment that does not
close.

An increment whose out-of-balance force is inside the tolerance but whose plastic stress points
have not settled -- the local criteria count the points whose stress still moved between two
iterates -- used to be iterated on with the tangent stiffness. On ground that has softened to a
few percent of its elastic stiffness a tangent step is a large one: the points went on moving
from one branch of their yield surface to another, the count never fell, and the increment was
cut back to nothing and reported as the limit of the phase. Such an increment is now finished
with the elastic stiffness, whose corrections are as small as the force still out of balance.

An increment that is simply too large was first retried three ways -- each from its beginning --
and only then halved. One still open after forty iterations is now halved at once, and the
retries are kept for the increment that cannot be made smaller.

Measured on an embankment built undrained on soft clay, both soils Hardening Soil: on a 1 m mesh
the construction phase did not finish in 100 minutes -- the increment in which the fill's core
reaches its unconfined strength took 162 iterations, then 264, then a third pass, where its half
takes 27 -- and now finishes with a factor of safety of 1.392 (1.402 on a 2 m mesh, 1.396 on
1.5 m; a published 1.4). The 2 m run takes half the time it did.

### The pore fluid of an undrained hardening soil follows the skeleton

The pore fluid of an undrained Hardening Soil material was sized at the reference pressure and
kept there. Kw/n is a multiple of the skeleton's bulk modulus -- that multiple is what makes the
undrained Poisson's ratio 0.495 -- and the skeleton's modulus follows the stress: under a few
metres of soft clay it is ten to a hundred times lower than at the reference pressure, so the
water was that much stiffer against it than the model intends. The fluid now follows the
skeleton at every stress point, as it does for Soft Soil (K2D-M002 says how it is derived). On
an embankment built undrained on soft clay the factor of safety is 1.4020 with it and 1.4019
without: the answers this was measured on barely move.

### Soft Soil can be undrained

Soft Soil and Soft Soil Creep with Undrained (A) used to be refused: the pore fluid is sized from
the skeleton's stiffness, and this model has no Young's modulus to size it from -- its bulk
modulus is (p' + c cot phi)/kappa*. The fluid now follows that law, stress point by stress point,
so the water is as much stiffer than the skeleton at 10 kPa as at 200. On a confined column of
soft clay under a surcharge of 20 kPa the vertical effective stress changes by -0.06 kPa when it
is loaded undrained (by -20.000 drained) and the column settles 3 mm instead of 683; left to
consolidate, it ends within 3% of the drained settlement. Undrained (B) stays refused for these
models: it replaces the strength by c = su with phi = 0, and the model is written with the
friction angle throughout.

### A region can change its material between phases

Ground that is improved, a fill that sets, a clay that is undrained while it is loaded and drained
in the long term: a phase can now say which material a soil region is made of in that phase
(`materials={region: material}` in a script, `pmat` in the file), and it stays so in the phases
after. The region keeps the stresses and the excess pore pressure it carries. The state that
belongs to the constitutive model -- a pre-consolidation pressure, a shear hardening, a
small-strain history -- is the soil's stress history and is kept where the new material is the
same model, so an overconsolidated clay made undrained is still overconsolidated; where the
model itself changes it starts from the stresses. The new material's unit weight applies from
that phase on. The phase says so (K2D-A028). On a confined elastic column loaded twice, the
first load meets the soft material and the second the stiff one it was changed to: 66.6667 and
25.0000 mm, q H / E_oed of each to the last digit, and changing the material alone moves nothing.

### The automatic K0 starts from the model's own K0nc, and Soft Soil follows its own law

The automatic K0 was Jaky's 1 - sin(phi') whatever the model. The hardening and soft-soil models
carry K0nc as a parameter -- it sets their cap -- and an entered value was ignored by the initial
stresses: a soil with K0nc = 0.65 started at 0.5 under a cap built for 0.65. The automatic K0 now
starts from the entered K0nc, with the overconsolidation terms on top as before. A model with no
friction angle to read (linear elastic, Hoek-Brown) takes 0.5 and, for rock, says so (K2D-A017).

Soft Soil and Soft Soil Creep:

- The elastic law is the model's, in p' + c cot(phi): K_ur = (p' + c cot phi)/kappa*. The shift
  was missing, so a soil with cohesion was softer in unloading and reloading than the model says,
  by (p' + c cot phi)/p' -- a quarter at 10 kPa for c = 1 kPa and phi = 22 degrees. With c = 0
  nothing changes. The hardening law stays in p_p itself, as the model writes it.
- The pre-consolidation pressure of an overconsolidated soil is that of the state it was
  consolidated to -- the vertical stress times the OCR (or plus the POP), the lateral stresses
  K0nc times that -- not the measure of the current state times the OCR, which gave a cap about
  a tenth too large at OCR = 2. Normally consolidated ground is unchanged.
- In a script the unloading/reloading Poisson's ratio defaults to 0.15, the model's own, instead
  of the 0.2 shared with the hardening models.

The springs of an embedded beam row are derived from the soil's own shear stiffness, as an
interface's are: beside a hardening soil they used the unused Young's modulus box.

### A consolidation phase places its load over its duration, and its sides drain

A consolidation phase with a duration put everything it activated on at its first instant: the
whole undrained pore pressure at once, and at the end of the phase what was left of it after all
of the phase's time. A fill is not placed in an instant. The load is now placed linearly over the
duration, so a phase ends with the pore pressure its loading generated less what drained while
it was being placed -- and the first time step no longer has to carry the whole stage, which is
the step an elastoplastic consolidation failed on. The pore pressure a phase inherits from the
one before it is a state, not a load, and is still taken up at the first instant: an undrained
loading followed by 400 days of consolidation settles 2.38 + 87.30 = 89.68 mm on the 10 m column
whose closed form is 89.66. `load_at_start` (`cloadstart` in the file) puts the load on at once,
for the column of the one-dimensional theory; a phase that ends on a state is unchanged.

With no flow boundary declared, the excess pore pressure used to drain through the ground
surface only, the sides sealed. The sides of a model are ground that goes on: every exposed face
now drains except those that look down. A side that is a barrier -- an axis of symmetry -- is
closed by declaring the boundaries, as before.

### An interface is as stiff as the soil beside it, and takes the soil it lies against

Three things about a wall's interfaces were not those of the ground it stood in.

Their stiffness was derived from the Young's modulus and Poisson's ratio boxes whatever the soil
model, and a hardening material does not use those boxes: every interface beside one was given
G = 5000 kPa, the boxes' defaults. It is now derived from the soil's own stiffness -- E50 for the
hardening models, the stiffness they shear with.

Their virtual thickness was a tenth of the elements along the wall, so refining the mesh at a wall
stiffened the contact. It is a tenth of the model's element size, and the tenth can be entered
per interface (`iface_vt` in the file, `iface_thickness=` in a script): the thicker, the more the
contact moves before it slips.

And all of a wall's interfaces, on both faces, were those of the soil at its midpoint: a wall
through clay into sand slid against clay over its whole length. Each interface element now takes
the soil it lies against, on its own side.

On a strutted wall in soft clay excavated in three stages the wall's movement at mid-depth went
from 19.4, 51.1 and 155.2 mm to 22.7, 55.7 and 157.1, with the next entry included.

### A load increment is solved before it is accepted

An increment whose out-of-balance force was already inside the tolerance at its first iterate was
accepted there, with no displacement. With the 1e-2 tolerance of the hardening models and a stage
load small against the weight of the model -- one fortieth of an excavation stage -- that is every
increment until enough load has gathered in the residual, and up to that much of the stage was
never applied. The displacement of such a phase came out too small. An increment is now solved
once before it can be accepted; only one that is in balance to rounding is kept as it stands.

### A footing pushed into sand reaches its settlement on a finer mesh again

A rigid footing pushed 50 mm into sand (axisymmetric, 15-noded, psi = 0) stopped at 13% of its
settlement on a 0.25 m mesh, where it had reached all of it before; the 0.4 m mesh was unaffected,
and it was the only one tested. An increment that Newton cannot close is retried with the elastic
stiffness before it is made smaller -- a remedy brought in for load-driven phases, initial phases on
sloping ground among them. Under an imposed displacement it is the wrong one: the slow elastic
iteration sent the phase down a different path. It is no longer used there. The footing: 604.5 kN
on the 0.25 m mesh (618.5 on 0.4 m, a published 588), and the 6-noded 0.15 m mesh, which stopped at
65%, now runs to the end. The 0.25 m case is in the tests.

### 6-noded elements in an axisymmetric model say which way their collapse load lies

A footing pushed into sand in an axisymmetric model gave 656 kN with 6-noded elements on a 0.4 m
mesh, against a published 588 and 619 with 15-noded ones -- 12% on the unsafe side, and nothing
said so. Plastic flow at constant volume constrains every stress point, and in axisymmetry the
hoop strain adds a constraint the quadratic triangle cannot meet; the error falls with the mesh
(638 at 0.25 m, 615 at 0.15 m) but slowly. Such a phase now carries a warning (K2D-A027) that
names the measured error and asks for 15-noded elements for a limit load.

### A staged phase can take its pore pressures from a steady groundwater flow

A Plastic or Safety phase can now take its pore pressures from a steady-state groundwater flow on
the ground active in it, instead of the water line (`flowpore` in the file, `pore_pressure="flow"`
in Python, a check box on the Studio's Water tab). Below the water line every boundary of that
ground -- the model's sides and the faces an excavation exposes -- is held at the line's head;
above it is a seepage face; the base is closed; and a wall or interface made impermeable closes its
line. Water then flows under a cut-off wall towards a pit pumped down to its floor, and the pore
pressure in front of the wall is that of the flow rather than of two separate hydrostatic levels.
This was refused before: a flow result could only drive a single phase, and with a wall it reached
the ground from one side only, so the difference of water pressure on the wall was missing from
its forces. The flow is now solved on the deformation mesh itself, each side of a wall on its own
nodes. Checked three ways: with one water line it is hydrostatic to 4e-13 kPa; with a full cut-off
down to the closed base each side is hydrostatic at its own line and the wall carries the moment
of the same two lines under the phreatic path (22.96 against 22.88 kNm/m); with a wall stopping
half way, the head rises in front of it and falls behind it. A cut-off wall ending on the model's
boundary no longer leaks through its shared end node.

### A region can have a water level of its own, or interpolated water

A Plastic or Safety phase can now give a soil region a pore pressure of its own (`pwater` /
`pwlevel` in the file, `region_water=` in Python): a hydrostatic level of the region's own -- the
pit lowered to its floor while the ground around it keeps the water table -- or a pore pressure
interpolated linearly with depth between the regions above and below it -- the ground between a
pit floor pumped dry and the toe of its wall, where the water arriving from outside meets the
pumped level. That is how a design analysis usually states the water under an excavation without
solving its flow, and a benchmark specified that way could not be rebuilt before: with the water
hydrostatic from the floor instead, the passive ground at the toe carried about twice the
effective stress it is meant to. On the ERTC7 propped-wall benchmark (Schweiger 2006), whose
reference study states its water exactly this way, the factor of safety was 12% above the
published curve with the water hydrostatic from the floor; with the interpolated water it lies
within -2.1% and +1.8% of it for embedments of 2 to 4.5 m (1.337 against 1.34 at 3 m), and the
wall moment and strut force within 2% from 3 m down. Checked against the closed forms: the interpolated pore pressure
is linear from the floor to the water table's value below to 1e-4 kPa, a region's own level is
hydrostatic from it, and a region left on the water table is the water table bit for bit.

### Soil dug out in a phase whose pore pressures come from a flow no longer weighs on the floor

A phase whose pore pressures come from a groundwater flow assembled the self-weight of every
element, the excavated ones included, so part of the weight of the soil dug out of a pit stayed on
the shared nodes of its floor: an excavated column heaved 10.01 mm where one-dimensional unloading
gives 10.70 (and the same phase on the water table 10.68). It now leaves the excavated elements
out, as the water-table path always did. Every staged phase with flow pore pressures and an
excavation is affected.

### Hydraulic heave is checked below exposed ground

A phase whose pore pressures come from a groundwater flow now checks the ground below every
upward-facing exposed surface standing in water -- an excavation floor, flooded or pumped dry -- for
hydraulic heave. Column by column, at every depth, the submerged weight of the ground above is set
against the water pressure in excess of hydrostatic: the stress form of the heave limit state,
which is the critical gradient over the gradient averaged over that depth. The least ratio and
where it was found are reported (K2D-A024), and a ratio under 1 is a warning (K2D-A025): the
ground there has no vertical effective stress left and lifts, which a converged deformation
analysis need not show -- a soil with cohesion can carry the tension over a thin layer. On a
uniform upward gradient the ratio is gamma'/(gamma_w i) to round-off at every depth, on 6- and
15-noded elements; in front of a cut-off wall with a pit pumped dry it is found beside the wall at
the depth of its toe, 4.55 against 4.56 read by hand from the heads on the wall's face. Side
friction and soil strength are not counted, and no partial factors are applied.

A node on an impermeable wall now takes the water table of its own side. A table drawn with a
step at the wall gave a node on the wall itself the mean of the two levels: the corner of a pit
floor pumped to -4 m beside a wall holding -1 m behind it was held at -2.5 m, and the exit gradient
beside the wall grew without bound as the mesh was refined (a heave ratio of 2.13, 1.08 and under 1
on elements of 1, 0.5 and 0.25 m; now 4.55, 4.25 and 4.04, converging on the column at the wall).
On the tutorial's anchored excavation with flow pore pressures both pit corners had been held 3 m
above the pit's water: its factor of safety on 6-noded elements moves from 1.615 to 1.726, and on
15-noded elements the Safety phase, which could not be started because the dewatering phase before
it collapsed at 99% of its load, now runs to 1.555 (published about 1.49). The heave ratio there
is 6.36, below the pit floor at the depth of the walls' toes.

### A script phase keeps the water of the phase before it

In a Python script, a phase that named no water fell back to the project's water table, where a
new phase in the Studio is a copy of the one before it, water included. A Safety phase after a
dewatering step therefore filled the excavation again -- and the search took that for a change
of its own and applied it at full strength, the water in the pit holding the walls back. On the
tutorial's anchored excavation (6-noded, 2 m mesh) the factor of safety was 1.847; with the
dewatered water kept it is 1.597, against the published 1.49. A phase that names no water now
keeps the previous phase's; `water="project"` returns to the project's own table.

### The factor of safety follows the strength path by arc length: faster, and no longer set by the stopping rule

The strength-reduction search prescribed the factor at every step and solved for the equilibrium
there. Near the factor of safety that equilibrium exists but lies far along the developing
mechanism -- on a c-phi slope the last 0.5% of the factor carried the displacement from 9.5 to
78 mm -- and a step asked for it either crawled for hundreds of iterations or was given up. Those
steps were most of the run time, and what the search reported depended on how tightly its steps
were solved: 1.3903 at a tolerated residual of 1e-3 and 1.3838 at 1e-4 with small steps, 1.3874
with the steps it used -- every one of them above the converged answer, the unsafe side.

The factor is now an unknown of the solve. Each step prescribes the length of the displacement
increment (an arc-length constraint) and solves for the displacement and the factor together, so
the steps stay small and well posed up to the limit, and the factor of safety is read where the
path stops rising: at its peak, or where it has risen less than the resolution over each of two
successive doublings of the displacement. On the same slope the factor is 1.38680, 1.38665 and
1.38680 at 1e-3, 1e-4 and 1e-5, in 19 s instead of 48; on 15-noded elements 1.3658 in 101 s
where the old search took 418 s for 1.3731 (1.3624 at 1e-5). The tutorial embankment's drained
factor is 1.8305 (published about 1.8) in 171 s instead of 264. A path that cannot be started
falls back to the prescribed search.

The continuation is used for ground without structures. A plate, anchor, geogrid or pile keeps
its own strength while the soil's is reduced, so once the ground has formed a mechanism the
structure takes the load on and the path keeps rising with the displacement -- there is no
plateau to read: on the tutorial's anchored excavation it walked to 3.8 m of displacement for a
factor of 1.852 where the prescribed search stops at 1.847 in a quarter of the time. With
structures the prescribed search is kept. A factorization is now also reused when the same matrix is
solved again, which the arc-length step does for its second right-hand side.

### Hardening Soil's cap hardens with the cohesion shift of the model's own law

The cap of the Hardening Soil model hardens with the pre-consolidation stress shifted by
c cot(phi), like the model's other stress-dependent laws:
d eps_v^pc / d p_c is proportional to ((p_c + c cot phi)/(p_ref + c cot phi))^-m. The shift was
missing, which made a cohesive soil's cap too soft at low stress -- by 7% at 50 kPa and 26% at
20 kPa for c = 5 kPa, phi = 29 degrees, m = 0.8. The cap parameters are calibrated on the
oedometer through the same law, so they follow it; a cohesionless soil is unchanged bit for bit.

### An increment that cannot reach its tolerance is no longer ground for hundreds of iterations

An increment whose Newton iteration stagnates is kept at its best iterate when that is within
1e-3 of the load scale (K2D-A019) -- but only after every retry had been spent on it, and on a
phase whose tolerance is out of reach that is every increment: an embankment built on clay
(Mohr-Coulomb, 15-noded) took 3478 iterations for 20 increments, each accepted at about 1e-4 in
the end, in 396 s. An increment whose best iterate is already inside that band and whose residual
has not fallen by a tenth over its last thirty iterates is accepted there, on the same terms and
counted the same way. The embankment now takes 79 s for a displacement field the same to 0.004%.
An iteration that is still converging never meets the rule, so every phase that reached its
tolerance before reaches it now: the Hoek-Brown tunnel, the Prandtl footing and the sliding block
are unchanged.

### A prescribed displacement is imposed once

A prescribed displacement was imposed again in every phase in which it stayed active: a footing
pushed 10 mm, followed by a phase that changed nothing, settled to 20 mm. It is now imposed by the
phase that activates it; a following phase in which it stays active holds its nodes where they
are, as an unchanged load adds nothing, and a phase that switches it off releases them.

### A node moved and then held stays where it was moved to

A node on a fixed boundary that a prescribed displacement moved (a tunnel contraction reaching
the symmetry axis, a footing edge on the model side) stays where it was moved to once the
displacement is switched off, and the boundary holds it. The structures read such a node as if it
had never moved: a plate through it saw it spring back, and a plate installed afterwards was
reported against the moved node with nothing taken off for its installation -- a tunnel lining
installed after the contraction showed +51 800 kN/m in its element at the axis next to -583 kN/m
in its neighbours. Structures now read the held node at its position, and a plate installed on
it starts stress-free.

### Excavated regions can be kept dry

An excavated region below the water table holds the water standing in it, and that water presses
on the ground around it -- right for a pit left to flood, wrong for a tunnel or a pit pumped dry,
where it pushed the crown up and turned the surface settlement over a shield tunnel into a heave.
A phase can now name excavated regions as dry (`dry` in the file, `dry=[...]` in Python, `wet=[...]`
to let the water back in): their faces carry no water. The floor of a pit dug 4 m below the water
table heaves by (gamma_sat - gamma_w) h H / E_oed with the water left in it and by
gamma_sat h H / E_oed kept dry, both to the closed form.

### Soft Soil stands still in its own K0 state with the water table at the surface

Soft Soil and Soft Soil Creep hold their bulk modulus at no less than that of 1 kPa of mean
stress, so that the stiffness p'/kappa* does not vanish at a stress-free point. The floor also
clamped the STRESS: a point committed below 1 kPa was handed back at 1 kPa under zero strain. With
the water table at or above the ground surface the effective stress near the top is below that,
so a K0 phase on level ground, whose seed is an equilibrium, moved the column by 1.26 mm. The floor
now applies to the stiffness only: above 1 kPa the kappa* law is exact as before, below it the mean
stress changes linearly with the floored modulus, continuous in value and slope. The column stands
still to zero, and every Soft Soil and Soft Soil Creep verification is unchanged.

### An initial phase on sloping ground no longer stops part way

The equilibrium step of an initial phase on non-level ground (and any plastic step that starts to
creep) raised a "stall floor" -- a share of the elastic stiffness in the plastic points' tangent --
while its Newton steps kept being cut. With non-associated flow the stiffened direction need not
descend either, and the floor kept climbing until the step was abandoned: a cohesive slope's K0
phase on 15-noded elements (phi = 30, psi = 0) stopped at 61% of its load with "did not converge",
and so did the same slope with a Soft Soil layer, on any tolerance. A step that finds no descent
now hands the next one back to Newton's own direction; both phases complete, and the slope that
the floor was introduced for still runs in seconds.

An increment whose Newton direction stops descending, every retry spent, is now re-entered once
with the elastic tangent at its plastic points, full steps and no demand that the residual fall
at every step -- the initial-stiffness iteration, a contraction that converges where a
non-associated Newton direction on a non-smooth surface does not descend, but whose residual may
rise for an iterate on its way. Only increments Newton has already given up on take it. On a sweep
of initial phases (Mohr-Coulomb and Soft Soil; one and three layers, level and sloping ground, with
and without water, K0 and gravity loading, 6- and 15-noded elements: 64 models), the failures fell
from 19 to 4, all of them K0 phases on sloping ground with 15-noded elements. Those four are now
solved too: a K0 phase on non-level ground whose equilibrium step does not close is solved again by
gravity loading -- the procedure for non-level ground, which ran every one of them -- and the run
says so (K2D-A022). All 64 initial phases complete; the Hardening Soil and HS small runs of the same
sweep had no failures.

### A model with no boundary conditions drawn gets the standard ones

A model with no boundary conditions was refused ("No boundary conditions"). It now gets the
conditions a soil model is built with unless there is a reason for others -- the base fully fixed,
the left and right sides fixed horizontally, every other edge free -- and the run says so
(K2D-A021). The result is the same field, bit for bit, as the same model with those conditions
drawn.

### A factor of safety for Hardening Soil and Soft Soil, in minutes and independent of the path

Started from the phase before it, a Safety phase could reduce the strength of a Hardening Soil or
Soft Soil material, and the answer depended on the route the search took: reduced together with
their hardening, the same embankment on Hardening Soil reported 1.864 along one sequence of
reduction steps and 1.676 along another, and a run took from five minutes to more than forty. In
a Safety phase these materials now take part as the Mohr-Coulomb material with their own c',
phi', psi and nu_ur and a stiffness taken from the state the phase starts from (K2D-A020). The
embankment of a published tutorial gives 1.821 on 15-noded elements against 1.8 published (and
1.823 with the same strength modelled in Mohr-Coulomb), with the Safety phase itself in about two
minutes.

### Plates that meet are joined rigidly

Each drawn plate numbered its own rotations, so two plates meeting at a node kept two independent
rotations there and every joint was a hinge, with nothing said. On a plate on elastic ground
loaded at mid-span the moment under the load was +56.65 kNm/m in one piece and -0.23 kNm/m with
the plate drawn as two pieces meeting there, its sign reversing a metre away; a tunnel lining,
which has to be drawn as the chords of its arc, carried its ring force with no bending at all. A
rotation at a node an earlier plate already turns is now the same unknown. KV-STR-015: the plate
drawn in two pieces gives the moments of the plate drawn in one, to 0.

A plate with interfaces (a wall on its own degrees of freedom) had the same joint, and worse: a
wall that is not vertical was split from the soil -- and built, with its interfaces -- up to a
metre past its drawn end, so a lining chord or a slab ending inside the soil was longer than drawn
and the force report of each piece ran into the next (a wall drawn from x = 3 to 9 m reported
forces to 9.75 m, and at a joint 59.49 kNm/m on one side against 52.37 on the other). A wall is
now split and built only where it is drawn, the ends of walls that meet are one point, and a plate
ending on a wall's end turns with it: the two-piece wall with interfaces gives the moment of the
one-piece wall to 3.8e-4.

### The ground surface drains in a consolidation phase with no flow boundaries drawn

With no flow boundaries declared, the drained boundary was the set of boundary nodes at the
highest elevation of the whole mesh, inactive regions included. With a fill not yet placed on
top, no active node drained and the ground consolidated as a sealed body; with it placed, only
its crest drained and the ground surface beside it stayed sealed. The upward-facing exposed
surface of the active ground now drains -- the ground surface, the fill's crest and slopes, an
excavation floor -- and the sides and base stay impermeable, as before.

### The initial phase is drained

The ground's own weight was carried over geological time, so an initial phase that only sets up
that state generates no excess pore pressure, whatever the drainage type. A non-level K0 phase's
equilibrium step (and a gravity phase's whole self-weight) loaded Undrained (A) soil undrained,
and every later phase inherited a pore pressure nobody had generated: a Safety phase on an
undrained slope differed from the same slope drained although nothing had been loaded yet. An
initial phase that also applies a load or a prescribed displacement keeps the soil's drainage
type for that load's response, as before.

### Free water presses on the ground it stands against

The pore pressure enters the calculation as a load on the soil it is in, and on the boundary of
that soil it leaves a term that only the pressure of free water standing against the boundary
balances. That pressure was never applied. Where the ground surface is at the water table the term
is zero and nothing was missing; anywhere else -- an excavation dug below the water table, a lake
or river bed, the upstream face of a dam -- the water's push on the ground was left out, and the
unbalanced term lifted the face by the full water pressure. An elastic column dug 4 m below the
water table heaved exactly twice the closed form (the dry pit's, although the water stood in it),
and a submerged excavation rebuilt from a published tutorial collapsed in its last stage with a
strut force four times too high. The water's pressure is now applied on the exposed boundary of the
active ground (a wall's split seam stays interior until one side is dug away), from the phase's own
water conditions; the column heaves 18.17 mm against 17.83 mm from the closed form, the difference
being the unit weight of water (9.81 against the closed form's 10).

### A prestressed anchor holds its lock-off force through the phase that installs it

A lock-off force was applied when the anchor was activated, and from that moment the anchor was a
spring: the wall's movement towards it in the same phase unloaded it. On a tie-back excavation
rebuilt from a published tutorial only 15% of each lock-off force (500 and 1000 kN) was left at the
end of its own phase, and the wall moved about four times as far as it should, with no message --
the unsafe side for the anchor and the wall alike. In the phase that installs it the anchor now
holds exactly its lock-off force, with no stiffness of its own (the jack), and is locked where the
wall has come to; from the next phase on it is a spring from that force. KV-STR-001 pins both
halves: the installation phase is the same field as no anchor plus the lock-off force as a point
load (0 m difference) with the anchor reporting exactly 500 kN on soft ground, and a later load
phase is the slack-anchor response with the anchor at 500 kN plus what the spring picks up.

### A consolidation phase carries plastic soil

An elastoplastic consolidation step was one Newton solve with no line search and no cut-back, held
to a criterion that compared the out-of-balance FORCE with the size of the flow term, which in a
dissipation step is tiny: on an embankment over Hardening Soil sand a step balanced to 7e-6 of the
forces in play oscillated for its whole iteration budget and stopped the phase. Soft Soil,
Hardening Soil and a Mohr-Coulomb fill activated in a consolidation phase all failed in the first
time step, at any number of steps, or ran for half an hour with no progress. Each balance is now
measured on its own scale (force against the forces in play, volume against the step's flow and
volume change), a line search judges both together, an iterate that stops descending below 1e-3 of
the force scale is kept and counted as a static phase keeps it, and a step that does not converge
is cut in two with its load increment, down to 1/64. The embankment's first construction phase now
runs with Hardening Soil sand, a Mohr-Coulomb fill and Soft Soil layers.

### Values entered without their switch are no longer dropped

Rinter, K0nc, K0, OCR and POP each take effect only with a switch that selects them, and the Python
interface set the value and left the switch alone: Rinter = 0.5 ran as a rigid interface and POP =
25 kPa as a normally consolidated soil, with nothing said. A value given now brings its switch with
it (OCR and POP together without saying which is refused), and a project file with a value its
switch leaves unread is warned about. The automatic K0 now reads a POP as well: it raises the
lateral stress point by point, K0nc (s'v + POP) - nu/(1 - nu) POP, where it used to reach only the
cap of the advanced models.

### A Safety phase starts from the ground the phases before it built

A Safety phase used to re-solve the whole model from the unstressed state -- the self-weight
applied again from zero, under each reduced strength -- and read its factor of safety from that
second construction rather than from the one the engineer had staged. Three things followed. A
prestressed anchor and an undrained soil had to be refused, because a lock-off force means
nothing on unstressed ground and an undrained soil loaded from zero carries its own weight in
the water. Every structure carried the settlement of the whole self-weight, including what
happened before it was installed. And the answer depended on how many steps that second
construction was given: a geosynthetic-reinforced wall that had just carried its own weight at
full strength was reported with a factor of safety of 0.80 on 8 load steps and 1.07 on 40.

A Safety phase that follows other phases now starts where they left the ground: their stresses
and pore pressures, the structures with the forces and plastic state they carried (a lock-off
force included), and their internal force held as the load while the strength is reduced from
there. Whatever the Safety phase itself changes -- an element removed, a load switched -- is first
carried at full strength, as a Plastic phase would carry it. An embankment on clay rebuilt from a
published tutorial now gives 1.82 drained and 1.42 undrained, against 1.8 and 1.4 published (the
undrained case was refused; from zero it had given 0.73). A prestressed anchor installed before
the Safety phase takes part, and on the sliding block of KV-STR-010 the factor meets the closed
form with the anchor at its capacity to +0.09%. Hardening Soil and Soft Soil, refused in a Safety
analysis because reducing their strength from zero is path-unstable, are accepted from a parent
phase. A Safety run with no phase before it, or with a structure installed in the Safety phase
itself, still starts from the unstressed state, and the refusals stay for that start.

### Cohesionless Hardening Soil no longer carries strength it does not have near the surface

The Hardening Soil model reads its stiffness at the minor principal stress, floored at a tenth of
the reference pressure so that the stiffness never falls to zero at a stress-free point. The
failure deviator read the same floored stress. So a cohesionless point shallower than that floor
carried the strength of the floor instead of its own: on a sand of 30 degrees, a deviator of 20 kPa
at zero confinement, an apparent cohesion of up to 6 kPa — on the unsafe side, and largest exactly
where a footing's bearing capacity is decided. The floor now applies to the stiffness only; the
failure deviator, the asymptote of the hyperbola and the mobilised friction angle read the stress
itself, and a stress-free point is at the apex, not on the failure plateau, so a compression takes
it into the cone by hardening. A point
sheared at 2 kPa of confinement now fails at 5.38 kPa, the Mohr-Coulomb value of its own stress,
where it used to hold 26.9 (a test pins it); a point pulled past the apex is returned to it, with a
zero tangent there.

With the true strength, shallow points reach the failure line, and two things that had never been
exercised there were wrong. A point brought onto the line other than by hardening along it — a
column loaded from zero stress, whose first elastic trial lies beyond the line — kept the plastic
shear strain it had before, far inside the surface its stress had reached; it now carries at least
the hardening the line implies. And the yield function that decides whether a trial is plastic
switched between the hyperbola and the Mohr-Coulomb excess at the failure deviator instead of taking
the larger of the two, so it jumped there. Both left a discontinuity in the response that held the
equilibrium iteration: the corpus oedometer, whose column starts stress-free, took 531 s to reach a
tolerated error of 1e-6 and now takes 10 s, and at the default tolerance 8 s instead of 52 s, with
the loading settlement closer to the closed form than before.

### Hardening Soil flows cleanly along the edge of its cap

The cap of the Hardening Soil model measures the deviator with q̃ = σ1 + (δ − 1)σ2 − δσ3, which
has an edge wherever two principal stresses are equal — and every point of a normally consolidated
ground starts on that edge, its two horizontal stresses equal. There the cap was one surface with
the averaged gradient, which does not hold the stress on the edge: loaded in compression, the
substeps alternated between the edge and a face and the stress returned jumped with the strain
increment. Measured on a strip footing on sand, 0.2% of the plastic points jumped at every
equilibrium iteration, deep ones included, and the iteration stalled. At an edge the cap is now two
faces with one hardening (Koiter), unless both shear planes of the corner are already flowing and
holding the edge — four surfaces in a three-dimensional stress space would be singular. The
preconsolidation seeded into an initial state is measured with the same q̃, so a state whose
principals are not pairwise equal (a slope, a history under an inclined surface) starts on the cap
the model reads rather than outside it.

The test that decides whether a substep yields in shear now reads the yield surface with the
stiffness of the stress it tests, as the flow and the drift correction already did. On a path that
lowers the minor stress the surface shrinks faster than the deviator falls, and the frozen test
called the start of every such substep elastic while its end point yielded.

A strip footing pushed into cohesionless sand with this model used to stop at 5 to 7% of its
settlement; it now reaches the full settlement, at a bearing pressure within 2% of the collapse
pressure the Mohr-Coulomb model of the same strength reaches on the same mesh.

### Hardening Soil dilates only where its rule says it does

The mobilised dilatancy of the Hardening Soil model follows Rowe's law above a mobilised friction of
three quarters of the friction angle's sine and is zero below it; a zero or negative dilatancy angle
is taken as it is. The threshold was missing — Rowe was applied as soon as it turned positive, i.e.
from the critical-state angle, which lies below the threshold whenever the dilatancy angle is large
(with phi = psi = 41 degrees the soil dilated from the first increment of shear) — and a negative
dilatancy angle, which the input accepts, was silently read as zero. The small-strain variant is
unaffected: below the critical state it contracts by its own rule, which a test already pins. A new
test pins the three branches.

### The line search finds its step in fewer tries

A rejected trial step used to be halved, so the step of 1/16 that a plastic mechanism typically
needs cost five evaluations of the whole internal force — on a Hardening Soil footing, three quarters
of the run time went into them. The next trial is now the minimum of the quadratic through the
residual at the start, its slope along the Newton step and the rejected trial, kept between a tenth
and a half of the last one. The acceptance test is unchanged, so no worse step can be taken, and the
worst case is the halving it replaces. The whole test suite runs in 671 s instead of 771 s.

### An iteration that creeps near a limit state is steered out of it

Near a limit state many stress points sit on the corners of their yield surfaces or change between
elastic and plastic from one iterate to the next. The tangent is exact on the side each point is
on, the full Newton step lands on the other side, and the line search can only take a sliver of it:
a 30-degree slope of sand on 15-noded elements spent 200 iterations at steps of 1e-3 before its
increment was accepted as stagnated, its force residual long under the tolerance and its local
criteria never. After three steps in a row cut to under a quarter, the plastic points' tangents now
get a share of their elastic stiffness — a tenth, doubling while the iteration still creeps, halving
with every good step — which steers the direction towards the elastic-stiffness iteration that
crosses those corners without having to see them. The residual the increment has to reach is
unchanged. The slope takes 32 s instead of 257 s with the same displacements to 1e-6; the
self-weight phase of a geosynthetic-reinforced wall on the same elements, which gave up at 78% of
its weight, now carries all of it; and a strip footing on cohesionless sand pushed to 0.3 m, which
stopped at 69% of the settlement, now stops at 77% on the same collapse pressure (283.0 against
283.5 kPa) -- past the limit load, where how far a non-associated mechanism can be followed depends
on the path the iteration takes.

An increment that runs out of iterations while its residual has not fallen by a tenth over its last
thirty iterates is now recorded as a stall rather than as a spent budget; it keeps its whole budget
first. The distinction matters because a budget is never published as a capacity: with the floor a
collapsing footing finds a hair of descent at every iterate, and without this the Prandtl footing of
the corpus refused to name its limit load.

### A phase that stops short keeps the state it last equilibrated

A phase that could not reach its full load published its load factor and nothing else — no
displacement, no stress, no reactions — and a footing pushed to collapse under an imposed
settlement reported a footing force of 0.0. The collapse load and the mechanism are read off
exactly that last equilibrated state, so the phase now carries it, and still reports that it
stopped: nothing continues from it. On a Tresca footing loaded past its capacity the reactions of
the kept state balance the equilibrated load to 1e-16, at N_c = 5.148 against 2 + π = 5.142.

### Hardening Soil hands Newton the tangent of the stress it returns where tension is cut off

A shallow Hardening Soil point pulled into tension is capped after the model's own return, but the
tangent handed to the equilibrium iteration was the one from before the cap: a principal the cap
pins to the tensile strength still carried the full elastic stiffness. Measured against a central
difference of the whole update on 20 000 random loaded states, every capped point disagreed, by up
to 1e5 relative. The cut-off's own Jacobian — closed form, since with its active set fixed the
return is affine — is now chained after the model's; the capped points agree to the round-off of
the difference, and a test holds one and two capped principals to it. The stresses are unchanged
bit for bit; only the iteration that reaches them is.

The same measurement, repeated from states a previous increment has actually reached rather than
from arbitrary ones, puts the rest of the Hardening Soil tangent where it should be: away from the
cut-off it matches the difference to within the increment's own nonlinearity, and the mismatch
shrinks with the increment. What remains is at shallow points, where one increment is large next
to the stress it starts from — a gap between the continuum and the algorithmic tangent, not a wrong
branch.

### A factor of safety in a minute instead of ten

The strength-reduction search re-solved the whole self-weight from the unstressed state for every
one of its twelve bisection trials, and half of those trials are collapses, each of which only
ends after the load stepping has been cut to its minimum with every retry exhausted. Measured on
the Griffiths-Lane slope, the collapsing trials took 50 of 58 s on 6-noded elements and most of
673 s on 15-noded ones.

* **The strength is now reduced from equilibrium.** The ground is brought to equilibrium under its
  own weight once, at full strength, and the strength is reduced step by step from each converged
  state — the step growing while it converges and halving when it does not. A step that fails is
  not taken for a collapse (a drop too large for the iteration budget fails just the same, and a
  first draft that kept such a failure as its upper bound reported 1.100 for a slope that stands to
  1.361); the factor is bracketed only when a step of 0.1% fails from the last equilibrium. The
  answers are the bisection's — Griffiths-Lane 6-noded 1.383 against 1.384 (44.6 → 13.5 s), on a
  0.5 m mesh 1.361 against 1.361 (151.8 → 60.6 s), 15-noded 1.351 against 1.351 (673 → 91 s) —
  and the refinement study KV-SLP-003 runs in 75 s where it took 535. A search with structural
  elements keeps the bisection, whose trials carry them.
* **Every reduction step is held to 1e-3.** A step accepted with a looser residual carries its
  out-of-balance force into the next, and along the path that compounded: a phase set to 1e-1
  reported the factor 26% too high. A looser setting now governs the self-weight solve only, and
  1e-1 and 1e-2 report the default's factor exactly (1.421024); K2D-A006 says so.
* `NewtonOptions::min_step_fraction` sets the smallest increment before a failing solve is a
  collapse (1/8, as before); a strength-reduction step, whose load is already applied, fails at
  its first cut.

### Hardening Soil integrates twice as fast, bit for bit

The active-set solve inside every substep built its planes, gradients and matrices on the heap —
several allocations per call, at every plastic stress point, several times per substep — and the
stiffness and strength laws re-evaluated the same trigonometry six times per call. Both are now on
the stack or hoisted, with the same operations in the same order: a set of oedometer walks and cap
calibrations runs in 3.0 s where it took 6.7, with identical output to the last digit. The
dewatering phase of the Berlin excavation (19 000 nodes, 27 000 yielding stress points) takes
233 s.

### Two solves at once in one process no longer hang or crash

The worker pool that parallelises element assembly served whichever caller came last: a second
thread starting a solve while another was running wrote its job over the first one's — the
function, the range and the count of workers still out — so workers ran pieces of the wrong job,
the count never came back to zero, and the process hung or crashed. `test_jobs`, which runs two
jobs concurrently, was the intermittent hang long listed as a known issue; eight copies of it run
side by side hung in most rounds and segfaulted in some. The pool now serves one job at a time,
and a caller that finds it busy runs its loop on its own thread; every loop writes only to
per-index storage, so the answer is the same either way. 160 runs, 0 failures.

### A wall with interfaces no longer leaves the ground beside it held by nothing

The soil nodes on the line of a wall with interfaces were counted as held by the wall, because the
wall runs through them — but such a wall stands on DOFs of its own and reaches those nodes only
through the interface. With the ground removed on both sides, and the interface with it, they
were held by nothing and still not fixed: a row of zeros in the stiffness. PARDISO passed it by
perturbing the pivot and handing the nodes an arbitrary displacement; the portable build refused
it, and the continuous-integration run failed on `test_wall_attachment`. A node now counts as
held by a wall only where the wall moves on its DOFs. The same test gives the same moments
(62.339243 kNm/m) on both backends.

### Hardening Soil is right away from the triaxial compression corner

Measured against published undrained triaxial element tests of four soft-clay parameter sets
(Borgh 2018, Chalmers; HSsmall, m = 1), K0-consolidated and sheared to 3% axial strain:

| | extension q before | extension q now | reference |
|---|---|---|---|
| Clay 4 | 194.1 kPa, p′ 93 → 213 | 73.6 kPa, p′ 94 → 84 | 75.0 kPa, p′ 81 |
| Clay 5 | 237.8 | 100.5 | 104.7 |
| Clay 6 | 323.8 | 132.2 | 138.7 |
| Clay 7 | 441.2 | 212.8 | 223.0 |

An undrained ψ = 0 extension of a normally consolidated clay came out 2 to 2.6 times too strong,
with the mean effective stress doubling on a path where it cannot rise. That is the stress path
under an excavation floor. Every earlier Hardening Soil test walked the triaxial compression
corner, where the model was right and the defects below could not show.

* **The shear flow is the Mohr-Coulomb planes' (Schanz 1998, Ch. 4; Benz 2007, Eqns 7.15-7.22).**
  One direction, (1, R, R), was used for every stress state: the compression corner's shape, with
  a plastic strain along the intermediate stress on every face, and a dilatancy of
  −3 sin ψm/(4 + sin ψm) per unit γp instead of −sin ψm, 27% short at ψm = 5.5°. The plastic
  potentials are now g_ij = (σi − σj)/2 − (σi + σj)/2 sin ψm on the faces, the corners are the
  meeting of two faces (compression: f13 + f12, extension: f13 + f23), decided by the rate as a
  Mohr-Coulomb return decides its edges, and both faces harden the same γp. On the failure
  plateau the yield gradient carries the strength's dependence on σ3.
* **The integrator reads which stress is which at every evaluation.** The major and minor stress
  were taken to be positions 0 and 2 for a whole increment, so once the values swapped inside it
  (the axial stress of an extension passing from major to minor) the flow went to the wrong
  directions.
* **The cap is measured with q̃ = σ1 + (δ − 1)σ2 − δσ3, δ = (3 + sin φ)/(3 − sin φ)** (Schanz 1998,
  Eqns 4.71-4.74) instead of von Mises. On the compression corner the two are identical, so the
  oedometer, the (α, β) calibration and every compression result are unchanged.
* **The cap's drift correction no longer pushes a stress out onto the cap.** It projected on |f|,
  so a state that had ended inside was pulled back out at constant pp along a compressive
  gradient, raising p′. It now corrects only a state outside, along D_e·g with the hardening.
  Integration tolerance 1e-5 and 1e-7 give the same answer (94.5 / 94.5 kPa).
* **Every hardening step is consistent to first order.** The yield function of a shear plane
  depends on σ3 through E_i, q_a and E_ur, and its gradient left that out, so γp ran ahead of the
  stress on every substep and the answer moved with the load-step count: γp by 4% between 4 and
  400 steps of the same oedometer, in this release's predecessor as well. It now moves by 6e-6,
  the returned tangent equals the finite-difference one (42 451 against 42 454 kPa on an
  oedometer increment, non-symmetric terms included), and a Berlin sand triaxial reaches its
  failure plateau at the shipped tolerance (100.000% of q_f, where 99.69% had been recorded as
  integration error). With it: each surface's consistency starts from the value it has, not
  from zero; a substep that starts inside is split where it meets the surface (Sloan, Abbo &
  Sheng 2001); a step that breaks both orderings takes the corner whose result keeps its own
  ordering, instead of pairing the minor stress with one of two equal ones; and the drift
  correction is handed which two directions are equal rather than a corner type read in another
  ordering. Before these, a weightless oedometer split its two equal lateral stresses at
  round-off, the response jumped by 1e-3 when a lateral strain of 1e-8 changed sign, and a
  1e-6 run that takes 26 s took over ten minutes.
* **The committed stress is projected on the trial's principal frame** before the strain
  increment is recovered, instead of both being sorted and subtracted component by component.
* **The initial pp and γp come from the pre-consolidation state** — σ′yy,c = OCR σ′yy (or
  σ′yy + POP), laterals K0nc σ′yy,c — instead of the current state scaled by OCR, and a POP is
  no longer turned into a ratio for every component. An overconsolidated soil is elastic in shear
  up to the deviator it has already carried.
* **The cap is calibrated on the normally consolidated line itself.** The oedometer the (α, β)
  calibration bisects on started isotropic at 2% of p_ref and read K0 at p_ref as though the path
  had forgotten that start. It had not — at the same α the reading was 0.405 from a 2 kPa start
  and 0.453 from a 20 kPa one — so every calibrated cap carried an arbitrary starting point, and
  Berlin sand's K0nc = 0.38 was recorded as "outside what the model can reproduce" with α parked at
  80. Started on the K0nc line, the probe reaches 0.38 exactly (α = 2.16), and a normally
  consolidated oedometer walk holds K0nc at every stress level (0.4264 at 50, 100, 400 and
  1000 kPa for c = 0). The cap now does what it is for in a triaxial test from an isotropic normally
  consolidated state: the secant at q_f/2 is 9.8% softer than E50 there, while the shear hyperbola
  alone reproduces E50 to 0.01%. And Eoed is read AT p_ref, as the response to a small strain
  from a state landed on p_ref, where it was the secant over the step that crossed it — a secant
  spanning ~6% of p_ref, i.e. the tangent a few percent above it, which calibrated every cap about
  1.5% too soft. An oedometer from rest now closes on the closed form as the stress rises (−1.36,
  −0.80, −0.49, −0.30% over 50–800 kPa, the remainder being the start from rest), where it
  levelled off at −0.7%.
* **The cap's hardening modulus reads p_c through the same floor as the stiffness** (10% of
  p_ref). It fell to zero at p_c = 0, so at a stress-free point — the first increment of every
  gravity phase — the first load flowed almost perfectly plastically on a cap of zero size and
  pulled the lateral stresses into tension under a vertical compression (K0 < 0); an equilibrium
  iteration on that stalled at a relative residual of 2e-2 and a Hardening Soil test with a
  tension cut-off ran past 3000 s. Below the floor the hardening law continues as its tangent
  line; p_c itself is not raised, so no shallow point is given an overconsolidation it never had.
* Compression is within −5 to −8% of the reference and systematically low; declared in KV-CST-019.
* New tests `test_hs_multiaxial` (KV-CST-018 flow rule on faces and corners, KV-CST-019 the
  element tests); KV-CST-011's check on the cut-off was a ratio of two recovered stresses that the
  corrected flow moved from 2.0 to 1.17 with the cut-off as wired as before, and now compares the
  two on the refined mesh (2.8).

### A script phase no longer brings back what an earlier phase excavated

Found by rebuilding a published deep-excavation benchmark (Berlin sand, anchored diaphragm wall
with interfaces, dewatering, four excavation stages) from a Python script.

* **Activation was not inherited by the script surface.** `prj.phases` promises that a phase lists
  only what it changes, but a phase's activation vector was written only when that phase toggled
  something of that class — and to the engine an empty vector means *everything active*, not *as
  before*. So the first phase after an excavation that did not mention regions put the excavated
  soil back, and any phase that did not mention structures switched on every structure not yet
  installed, prestressed anchors included. Measured on a linear-elastic block: the excavation
  heaved its floor 13.3 mm and the next phase, which changed nothing, settled it 12.9 mm back. In
  the benchmark every anchor stage undid the excavation before it and the wall moved *away* from
  the pit. A vector is now written whenever the inherited state has anything switched off; a
  file with no deactivation is byte-identical to before. Files written by the GUI were not
  affected. New test `test_python_phase_inheritance` (fails on the old surface: 14.2 mm of
  movement in a phase that changes nothing).
* **The vertex floor is added to every plastic stress point, not substituted for a zero tangent.**
  In the same benchmark the first anchor's prestress phase stopped at 90% as a "collapse" with the
  soil elastic everywhere except at the top of the wall, where cohesionless points sat in the
  corner and tension regions of the yield surface: rank-deficient tangents that are not zero, which
  the floor did not touch. With the floor added to every plastic point once a tangent has been
  refused, all eleven phases complete. Only the tangent changes; the residual is the same.
* The Python `displacement` of a result now says what it is: what that phase did — the increment
  of a plastic or consolidation phase, the peak of a dynamic phase, the mechanism of a Safety
  phase — not the total since the start.

### Cohesionless soil is no longer reported as collapsing, a pushed footing reaches its settlement, and an undrained Safety run is refused

Two plane-strain / axisymmetric phases that stopped while nothing had failed, found by rebuilding
external worked examples:

| | before | now |
|---|---|---|
| smooth rigid wall moved 4 mm from c = 0, φ = 41° sand, 20 steps | "incremental limit (collapse) load" at 1% | full movement; thrust 141.29 kN/m against Rankine's 145.94 (−3.2%, tri6) |
| same, ψ = 20° | "collapse" at 0% | full movement; 142.04 kN/m (−2.7%) |
| rigid footing pushed 50 mm into sand (axisymmetric, 15-noded, c = 1, φ = 30°, ψ = 0), 20 steps | out of iterations at 20% | full settlement in 31 s |
| same, 0.25 m elements | out of iterations at 80% | full settlement, `K2D-A019` |

* **A singular tangent is first asked whether it is a vertex.** With c = 0 the Mohr-Coulomb apex
  is the stress-free state and its tangent is exactly zero, so a node whose stress points had all
  lost their stress carried no stiffness and the factorisation refused. Three halvings later the
  phase printed a collapse load. A refused solve now retries the same increment with the tangent
  of every such point replaced by 10⁻³ of its elastic operator (`kVertexFloor`), and the phase
  keeps that floor. Only the tangent changes — the residual an increment must reach is the same —
  and the answer shows it: 141.29 kN/m at a floor of 10⁻³ and at 10⁻¹. A genuine mechanism still
  ends in a refusal or a stall.
* **A prescribed displacement is predicted into the ground.** It entered each increment through
  the fixed nodes alone, so the first iterate strained one row of elements by the whole increment
  (out-of-balance force 22 times the reference) and Newton spent its budget walking out of
  stresses the converged field never has. Each such increment now starts from the linear
  response to it, `K_ff du = r − Δλ K_fp ū`, with `K_fp ū` measured as the directional derivative
  of the internal force along the prescribed motion; the start is kept only when its residual is
  below the plain start's.
* **A stagnated increment below 10⁻³ is kept and said so (`K2D-A019`).** On refined meshes of
  perfectly plastic soil Newton can stop descending at ~10⁻³ with a fixed set of plastic points,
  far from the 10⁻⁶ the Mohr-Coulomb default asks for. An increment that has spent every retry and
  is about to be cut back is now kept at its best iterate when that iterate is below 10⁻³ of the
  load scale; the phase reports how many and the worst. Measured on the footing: 0.03% on the
  mesh where 10⁻⁶ is reachable. The Safety search is exempt (`K2D-A006`).
* **An undrained soil in a Safety run is refused (`K2D-G017`).** The search re-solves the ground
  from the unstressed state, so every trial loaded an undrained soil with its whole self-weight,
  put that weight into the pore water and reported a factor of safety far too low — with no
  warning. Measured on an embankment on clay rebuilt from a published tutorial (Mohr-Coulomb,
  c′ = 10 kPa, φ′ = 25°): 1.82 drained against 1.8 published, and 0.73 undrained against 1.4.
  The drained factor and the phase's "ignore undrained behaviour" still run; a short-term factor
  needs the search to start from the parent phase, which is not done yet. Undrained (B), whose
  strength does not read the effective stress, is not refused.
* New test `test_python_solver_robustness` (Rankine thrust; footing force the same at 20 and 60
  steps, −0.08%; the undrained Safety refusal and both of its remedies).

### A wall with interfaces carries what is drawn on it, and an excavation unloads it

An anchored sheet pile with interfaces on both sides, rebuilt from a published worked example
(two layers, two excavation stages, two prestressed tiebacks), came out with anchors of 1.7 and
1.3 kN where 35.0 and 88.0 were expected and the wall's largest moment at its free toe. The same
model without interfaces matched the example. Four defects in the wall-with-interfaces path,
each silent:

| | before | now |
|---|---|---|
| anchor forces (upper / lower), kN | 1.7 / 1.3 | 36.1 / 98.7 (example 35.0 / 88.0) |
| largest wall moment, kNm/m | 34.8, **at the free toe** | 46.4 (example 47.0) |
| moment at the free toe, kNm/m | −34.8 | −0.2 |
| largest horizontal movement, mm | −10.2 | −19.1 (example −20.2) |

(tri6, 0.5 m elements; the model without interfaces gives 35.0 / 99.1 kN and 45.2 kNm/m.)

* **The toe was clamped.** The driver fixed the rotation of the wall's toe, so the free end of an
  embedded wall carried a moment. The rotation is free now; the translations stay shared with
  the soil below the toe.
* **What was drawn on the wall acted on the soil beside it.** A wall with interfaces moves on
  DOFs of its own, and an anchor, a point load, a geogrid or a plate found its node by position
  — the soil node on the wall line. An anchor from the wall therefore tensioned the retained
  ground against itself and the wall never felt it. Every node of an active wall's line now maps
  to the wall's own translation (`AnchorElement::end_dof`, `GeogridElement::trans_dof`,
  a plate's `trans_dof`, the point load's equations).
* **The joints of an excavated side went on pushing.** When a phase removes the soil on one side,
  that side's nodes are orphaned and fixed, and the joints there stayed on: they held the wall
  against points fixed in space and kept pressing it with the earth pressure of the removed
  ground, so the excavation never unloaded the wall (the lower anchor ended in compression,
  −2.4 kN). A joint whose ground is gone on either side is now switched off for the phase
  (`InterfaceElement::active`) — no stiffness, no σn0, zero in the report — and keeps its place in
  the arrays so the phase chain's state stays aligned. This holds for soil-soil interfaces too.
* **A one-sided interface was built on both sides.** `iface_pos` / `iface_neg` were only or-ed.
  They are read now: positive is the right of the direction from the first point to the second
  (the side the Studio marks with "+"), and a side without an interface is bonded to the wall.

**An anchor ending on a free embedded beam pulls the beam.** A free connection is the grout body
of a ground anchor, but the anchor's end found the soil node there and the beam took only what
the ground passed to it (10.6% of the anchor's force against 99.6% now, KV-STR-014 (iv)).

**A pile base carried tension.** The embedded beam's foot spring took force both ways, although
the formulation it implements says the base "can only sustain compression, not tension"
(`docs/references/embedded-beam-formulation.md` sec 4): a pulled grout body held the ground by a
single point at its tip, where a Mohr-Coulomb soil can only fail. Rebuilt with embedded-beam
grout bodies (hinged), the anchored sheet pile above stopped at 98% of its last excavation
("mechanism"), while the same wall with geogrid grout bodies solved even with the clay's cohesion
lowered from 10 to 8 kPa. The base now lifts off in tension (reversibly; its plastic state is
kept), and the hinged variant solves: wall 46.4 kNm/m, 19.9 mm, anchors 28.8 / 68.7 kN. On a
linear ground a pile pulled up is now no longer the mirror image of the same pile pushed down:
the tip force pulled / pushed was exactly 1.0000 and is 0.0081 (KV-STR-014 (v)). And **an exactly
horizontal pile had its base at its connection point**: the base is the lower end and the
connection point the upper one, but for a horizontal pile the connection point is the lower-x end
and the base was simply the first point drawn, so a pile drawn left to right had both at its head
and none at its far end. The base is now the end that is not the connection point.

**Still open, measured in the same study.** With FREE embedded-beam grout bodies the anchored
sheet pile does not finish its last excavation (93% at 200 iterations per increment, 96% at
2000, the out-of-balance still falling; 98% "mechanism" with a 0.2 m body), and on a tri15 mesh
the model without interfaces stops at 87% (89% at 2000 iterations, still falling) where tri6
solves. Neither is a collapse load: the equilibrium iteration converges far too slowly near the
limit state, which is a numerical question of its own.

**Interface stiffness can be entered** (`.k2d` **v19**: `structs[i].iface_kn`, `iface_ks`,
kN/m³; 0 = derived from the adjacent soil as before, and written only when set). The derived
value follows the soil and the mesh; a joint whose stiffness is known now gets it. The Python
`plate()` and `interface()` take `iface_material`, `iface_kn` and `iface_ks`, so a wall's
interface strength no longer needs the project's private lists.

**KV-STR-014** (`tests/test_wall_attachment.cpp`) checks each against an exact oracle: the part
of a wall above an excavation on both sides is a cantilever, so its moment is statics at every
station (1.6e-11 under a head load, 3.0e-12 below an anchor); a one-sided interface on an
excavated side is the plain bonded plate (1.4e-11); two interfaces with one side gone are the
other one alone (1.8e-11); an entered k_s is read back as τ/slip exactly. With each fix undone in
turn the case fails: the clamped toe carries 94% of the peak moment (1.6% now), the anchor and the
load reach the wall not at all, the excavated-side joints leave the cantilever's moment 99% off.

### Regions that overlap, nearly meet or pinch out mesh — and a region closes against its neighbours

Several geometries a section drawing produces all the time either meshed wrong without saying so
or never finished meshing. Measured on two-region models, 0.5 m² target element area:

| geometry | before | now |
|---|---|---|
| two regions overlapping in area (union 88 m²) | 76 m² — the 12 m² overlap meshed as a **hole** | 88 m², the later region owns the overlap |
| a 6 m² lens drawn inside a 50 m² layer | 44 m², lens **gone**, a hole in its place | 44 + 6 m² |
| a corner typed 1e-5 m off its neighbour's | **did not finish** (killed after 30 min) | 1 ms, mesh bit-identical to the exact corner |
| a vertex 1e-4 m above a neighbour's edge | **did not finish** | 1 ms |
| a layer pinching out at 5.7° or 1.1° | **did not finish** | 1 ms, 45 + 5 and 49 + 1 m² |
| a single region with a 3° corner | **did not finish** | 1 ms |

Three causes. The inside/outside test counted every polygon edge even-odd, so where two regions
overlapped the count was even and the ground was "outside". Nearly coincident points more than
1e-6 m apart were kept apart, so a corner typed a hair off became a sliver 1e-5 m thick and 10 m
long that the quality refinement tried to resolve. And the refinement kept splitting the triangles
in the corner of an input angle narrower than its 20° bound, where no insertion can help, until its
step cap — with every step scanning all triangles, which is what made the cap take hours.

**The geometry is now noded once, to one tolerance** — 1e-5 of the model's diagonal, 1 mm on a
100 m section (`katai/geometry/planar_graph.hpp`): vertices that close merge, a vertex that close
to an edge splits it, crossings become vertices, and an edge two regions share is one mesh edge
instead of two. **Each piece is kept by who owns its two sides**, with the rule the driver already
used for material and phase activity: the last region containing a point owns it. A piece with
different owners on its two sides is a constraint, and the domain outline when one side has none;
a piece with the same owner on both — a boundary drawn over by a later region — is dropped. So
overlap and lenses are defined, and a gap wider than the tolerance (a 2 mm strip was checked) is
real and stays empty. **Corners narrower than the bound are left at their own angle**
(Miller, Pav & Walkington 2003, the rule Shewchuk's Triangle applies), and the mesh message counts
those elements, because a corner can also be narrow by accident. The refinement now keeps its
candidates incrementally and makes the same choices it made before. Meshed with the old and the
new builder side by side, **23 of the 26 corpus files come out bit for bit identical**. The three
that do not: the two with a shared boundary (`kv-dyn-003`, 80 → 104 elements; `kv-exc-001`,
267 → 266), because the shared edge is no longer a doubled constraint; and `kv-fnd-010`, same 478
elements, whose surface load ends inside the top edge — the edge used to be split at a point
interpolated along it, one unit in the last place off the load's end, and is now split at the
coordinate the load was given.

**Drawing.** A region drawn against existing ones closes along their boundary: start on a
neighbour's edge, click the free vertices, end on a neighbour's edge and right-click. Before, a
region above a layer whose top had several vertices had to be traced through every one of them
again, and drawn with only its end points its straight closing edge left slivers of gap and
overlap. **A new Split tool** cuts every region a drawn line crosses into pieces; a closed line
cuts out the shape it encloses. The pieces keep the material, the coarseness, the per-edge
conditions of the edges they lie on, and in every phase the activation of the region they came
from (`katai/model/region_edit.hpp`).

**Deleting a region left the phase activation flags in place.** A phase activates regions by
index, so deleting region *i* in the Studio shifted every later region down while its flags stayed
put: the excavation of region *i* passed silently to the region after it, in every phase. All four
delete paths now go through `katai::model::erase_polygon`, which moves the flags with the regions.

Verification: KV-GEO-002 (areas per region and boundary length of the mesh against the drawn
union, for all of the geometries above) and KV-GEO-003 (close, split, lens cut, phase flags).

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

**A structural element active in a Safety run was refused** (`K2D-G016`), by the input contract
at `initial.struct` / `phases[i].struct` and again by the engine, with no factor of safety
reported. That refusal did not last the release: the search was then given the structures, and
`K2D-G016` now refuses only a prestressed anchor — see *A Safety run solves its structures* below.

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

### A converged phase could read as a thousandfold miss

The command line printed `converged: force error 1.763e-07 of 1.000e-10 tolerated` for the second
phase of `KV-STR-005`, which had converged, and the Studio's HTML report put the same number under a
**warn** badge beneath a "Converged" row. Nothing was wrong with the run: the ratio its steps
stopped on was 1.763e-13. Two force ratios are measured at every accepted iterate, and the surfaces
printed the one that decides nothing: the stiffness-weighted error, ‖r‖ / (‖f_int‖ + CSP·‖f_const‖),
which is reported because it is the stricter reading near collapse but cannot yet gate — its
denominator has no floor tied to the applied load, and gating on it measured 5 failures in 142 fast
tests. A step stops on the other one, ‖r‖ over the fixed scale max(‖f_ext‖, ‖f_const‖, 1), and a
results file did not store it at all.

**Every surface now leads with the ratio the step stopped on, against the tolerance it was held
to,** and shows the stiffness-weighted error beside it, labelled as reported and not the stopping
test: the command line (a second line under `converged:`), the Studio's text and HTML reports (a
"Weighted force error" row with no badge — a reading that decides nothing cannot warn), and the
Python `Convergence` repr. `Convergence.global_ok()` and `global_recorded()` are new; `force_ok()`
keeps its meaning and is documented as informational. **`.res` v10** stores the stopping ratio.
A version 9 file still loads, and the ratio it never stored reads back as NaN and is shown as "not
recorded" — never as 0, which would be a perfect number nobody measured. `test_results_io` builds a
genuine version 9 file from a version 10 one and checks both halves of that; `test_report` checks
the order, the missing badge and the "not recorded" line; `test_convergence_criteria` checks that the
stopping ratio is recorded and within its tolerance on a converged run.

### A layered model consolidated differently depending on the order of its material list

A consolidation or fully-coupled phase stores water as well as moving it: the pores hold a volume
n/K_w per kPa of pore pressure, and the porosity n is each material's own (`e_init`). Both phases
took one value for the whole mesh — from whichever material came first in the project's list. For a
soft soil that is invisible, because its compressibility dwarfs the water's; for a stiff layer it is
not, and the symptom was an answer that depended on the order of a list. Measured on two 6 m layers
of porosity 0.2 over 0.6 under a 10 kPa surcharge, swapping the two materials in the list moved the
degree of consolidation by

| E_oed | largest change in U |
|---|---|
| 1 MPa | 0.02 percentage points |
| 50 MPa | 0.8 percentage points |
| 500 MPa | 6.8 percentage points, and the undrained pressure from 8.69 to 9.52 kPa |

**The storage term is now per material in both phases.** `KV-CON-004` checks it three ways: after a
first step far shorter than the drainage time each layer holds its own closed-form undrained
pressure, q / (1 + n E_oed / K_w) — 9.5238 and 8.6957 kPa in one column, measured to −0.013% and
+0.000%; the two listing orders give bit-identical answers; and the settlement history follows an
independent finite-volume solution of the two-layer equation written in the test to −0.09%, where
the same solution with one porosity for both layers is 0.97% or 3.76% away. A model whose
materials all share one `e_init` — including every model that never set it — gives exactly the
answer it gave before.

The Non-porous refusal in these phases gave "a single fluid stiffness" as half of its reason. That
half is gone; the refusal stays, because every element still receives a pore-pressure unknown, and
the message now says only that.

### A Safety run solves its structures

The strength-reduction search — the initial Safety procedure or a Safety phase — now solves the
phase's active structural elements together with the soil in every trial, and their self-weight
with them. What the reduction does to each follows from what a strength is. An **interface's**
strength is a soil strength (c_i = R_inter c′, tan φ_i = R_inter tan φ′), so it is reduced exactly
as the soil's: c_i / SRF, tan φ_i / SRF, and its tension cut-off / SRF. A **plate's** M_p and N_p,
an **anchor's** and a **geogrid's** capacity, and an **embedded beam's** skin and base resistance
are the structure's own declared capacities, and stay at their input values. The search starts
from the unstressed state, so an interface starts unstressed with it: the normal stress it is
seeded with in a K0 phase is not carried in.

`KV-STR-010` checks all three rules against closed forms, on the `KV-STR-002` block (4 m × 1 m,
W = 100 kN/m, joint c_w = 2.5 kN/m², φ_w = 26.6°) with the imposed slip replaced by a horizontal
load H on its left face:

| run | closed form | expected | measured |
|---|---|---|---|
| interface, H = 40 kN/m | (B c_w + W tan φ_w) / H | 1.501907 | 1.501636 |
| interface, H = 30 kN/m | same | 2.002542 | 2.002466 |
| interface, c_w = 0 | W tan φ_w / H | 1.251907 | 1.252173 |
| plate on the top, w = 10 kN/m/m | (B c_w + (W + w B) tan φ_w) / H | 2.002669 | 2.002466 |
| anchor pulling back, F = 10 kN, H = 40 | (B c_w + W tan φ_w) / (H − F) | 2.002542 | 2.002466 |
| the same anchor at 2 m spacing | (B c_w + W tan φ_w) / (H − F/2) | 1.716465 | 1.716187 |

Every one is inside a single interval of the search's own bisection, (3.0 − 0.4) / 2¹² = 6.35·10⁻⁴,
which is the band, and the same numbers come out of both linear-solver backends. The answers the
band rejects are far outside it: without the plate's weight 1.5019, with the anchor's capacity
reduced as well 1.7519, without the anchor 1.5019. The plate's weight and the same distributed load
give the same factor bit for bit, and so do the yielded anchor and a constant point force equal to
its capacity; an elastic anchor, which has no capacity, holds the block through the whole search and
the factor is reported as a lower bound. The band also sees geometry: an anchor drawn at mid-height
attached to the nearest node, 3 cm above its line, pulled 0.9° downward and read +0.123% against the
horizontal closed form, and −0.003% against the closed form written on the direction it actually
pulled along.

On the Griffiths & Lane slope, each element against the same mesh with it deactivated
(`test_safety_gui`; direction witnesses, not verification):

| element in the Safety run | factor of safety |
|---|---|
| none (all deactivated) | 1.04080 |
| a geogrid across the slip surface, EA = 1e5 kN/m | 1.20012 |
| an anchor across the slip surface, EA = 1e6 kN | 1.15315 |
| an embedded beam through the slope face | 1.18679 |
| a crest slab, weightless / at w = 150 kN/m/m | 1.04016 / 1.01667 |

An interface as strong as the soil, which before left the two sides of its line unconnected, now
gives 1.036. A model with no active structural element is unchanged bit for bit: `KV-SLP-001`
still reports 1.0103271484374998 and `KV-SLP-002` 1.3835693359375, with the same largest mechanism
displacement.

**One element is still refused: a prestressed anchor** (`K2D-G016`, narrowed to this). Its lock-off
force belongs to a ground that has already moved, and the search starts from the unstressed one — on
the unstressed mesh the force would pull on soil that carries no stress yet. Deactivating the
anchor in the Safety phase runs.

**A Safety run with structures says what starting from the unstressed state means for them**
(`K2D-A018`, a note): each element is present from the first increment of every trial, so it also
carries what the ground's settlement under its own weight does to it, including settlement that, in
the phases built before it, may have happened before it was installed. For the same reason no
structural force is reported for a Safety phase.

**The phase after a Safety phase now continues from the structural state.** A Safety phase has
never committed its stresses forward, but it did not pass the structures on either, so the phase
after it found a parent with no structural state and re-developed every structural force from zero
although nothing had changed. Measured on the slope: an anchor carrying 45.072858952 kN after
gravity carries 45.072858952 kN in a nil phase after a Safety phase.

### A design approach left four ground strengths at their characteristic value

EC7 DA1-C2 and DA3 divide the ground's characteristic shear strength by the M2 partial factors
(EN 1997-1 Table A.4: γ_φ′ = 1.25 on tan φ′, γ_c′ = 1.25, γ_cu = 1.4) and re-solve. Four strengths
were never divided. Each was silent, each erred on the unsafe side, and the report said the approach
had been applied:

| strength | measured under EC7 DA3 | on hand-factored input |
|---|---|---|
| every interface (sliding block, KV-STR-002 geometry) | 59.720 kN/m — the characteristic capacity, bit for bit | 47.743 kN/m |
| the cohesion gradient `c_inc` (footing, c_inc = 20 kPa/m) | 125.467 kN/m — only c_ref divided | 117.718 kN/m |
| Soft Soil (footing) | 61.458 kN/m — the characteristic run, bit for bit | 58.748 kN/m |
| Soft Soil Creep (footing) | 68.357 kN/m — the characteristic run, bit for bit | 64.781 kN/m |

The interfaces are built from characteristic strength before a phase applies its approach, and
nothing factored them afterwards; the gradient lives in a separate profile; the two soft soils keep
their failure line in parameter blocks of their own. **All four are now factored**, each by its own
material's rule, so a gradient or a joint beside an undrained material takes γ_cu.

**A joint beside an undrained clay also kept a friction angle the clay does not have.** Undrained (B)
and (C) solve the soil as a Tresca material (c = s_u, φ = 0), and the input contract says the entered
φ "is ignored for strength". The interface builder read that box anyway: with it left at 26.6° the
sliding block's joint carried 59.72 kN/m, where B s_u = 10. An interface now takes the strength its
material is solved with — c_i = R_inter s_u and φ_i = 0 beside an undrained clay.

`KV-STR-011` holds the joint to closed forms on the checked-in sliding block: under DA3 and DA1-C2,
B c_w / γ_c′ + W tan φ_w / γ_φ′ = 48.061 (measured −0.66%, inside the case's 2%); an Undrained (B) or
(C) joint, B s_u = 10 (−1.04%) whatever its friction box holds; and the same joint under DA3,
B s_u / γ_cu, with the characteristic run's bias to 1e-9. `KV-FND-015` runs each strength twice —
under DA3, and with no design approach on strength factored by hand — and the two agree bit for bit,
while each differs from its characteristic run. For the soft soils the design approach is applied in
both phases of that comparison: applied in the loading phase alone, Soft Soil differs by 1.06e-5,
because the preconsolidation the first phase seeds reads the strength that phase was given.

The same Soft Soil Creep block was also missing from the Safety search's strength reduction. Safety
refuses both soft soils, so no factor of safety was affected; the block is reduced now, so that
lifting the refusal cannot bring the omission back.

### Activating one structure reset every other structure's force

A staged phase carries the structural state of the phase before it — each structure's total
displacement and committed plastic state — so a wall's moment continues from one construction stage
to the next. That carry was accepted only when the structure set was **identical**. Activating or
removing any structure — an anchor row after an excavation stage, a strut, a geogrid with each lift —
dropped it for every structure in the model, and the only statement was a sentence at the end of the
phase message. Measured on a wall behind a 50 kPa surcharge, in a phase whose only change was an
anchor of EA = 1 kN placed far from the wall:

| wall | before the phase | the phase, anchor not activated | the phase, anchor activated |
|---|---|---|---|
| plate with interfaces | max\|M\| 1.914704 kNm/m | 1.914704, nothing moves | **0.921088 (−52%)** |
| plate without interfaces | 17.889413 kNm/m | 17.889413, nothing moves | **10.867339 (−39%)** |

That is the ordinary anchored-excavation sequence. **The state is now matched by drawn structure.** A
structure present in both phases keeps its state; one activated in the phase is installed on the
ground as the parent phase left it, and reads only the displacement since its installation; one
removed releases its force as part of the phase's change. A prestressed anchor activated in a later
phase still applies its lock-off force in that phase. On the wall above, activating the anchor now
leaves the moment at 1.914704 and 17.889413, with nothing moving.

`KV-STR-012` checks the three rules against superposition on a linear model, where every phase is the
solve of its own change: activating an anchor in a phase that changes nothing else moves nothing and
leaves the anchor at zero force; loading after it gives the wall the sum of the two stages' moments
(1.6e-13) and the anchor the second stage's force alone (8.7e-13); removing it again changes the wall
by exactly the moment of its released force applied as a point load (1.5e-11); and a prestressed
anchor activated in the chain equals a slack one with its lock-off force as a point load (8.4e-13).
The nil-phase identity of the existing carry test is unchanged bit for bit.

### A wall with interfaces, or an interface, can be activated in a later phase

Both split the mesh along their line, and the split — with the element — had to exist in every
phase: a phase that left one out was refused. The everyday excavation sequence, a K0 initial phase
and the wall installed in the first construction stage, could not be modelled; the wall had to be in
the geostatic phase itself.

**The mesh is now split in every phase whether the element is active or not**, so every phase has the
same nodes. In a phase where the element is inactive nothing is built on the seam and its two sides
are tied into one unknown each: the ground is continuous and fully permeable there. An element
activated later is installed on the ground as the phase before left it (the rule of the previous
section), and its joint takes over the stress the ground carried across the line: the recovered
committed stress, normal n·σ·n as the joint's σ_n0 and shear t·σ·n as an initial slip, capped by the
joint's strength.

`KV-STR-013` checks it three ways. A wall or an interface inactive in every phase gives the
displacement field of the same line drawn as an inactive geogrid, which splits nothing, bit for bit;
so does a Safety run with the joint deactivated. K0 followed by the wall's activation equals the wall
present from the K0 phase: after a level K0 phase the recovered stress is K0 σ′_v, the seed the wall
has always had, and the surcharge phase after it matches to 2.3·10⁻¹⁵ m and 7.4·10⁻¹³ of the moment,
with a wall of 9.6 kN/m/m settling 1.349268·10⁻³ m under its own weight in the activation phase exactly
as it does in the K0 phase. Activated after a stage that is not geostatic — gravity loading and a
surcharge on Mohr-Coulomb ground — the recovered stress is not the discrete traction exactly, and the
activation moves the ground by 2.15·10⁻⁴ of the stage before it and gives the wall 0.00996 kNm/m;
without the shear seed the same activation moves 3.6·10⁻³ and bends the wall 0.47 kNm/m. On
linear-elastic ground the activation opens the joint where that ground carried tension across the
line, which a joint without tensile strength cannot do.

Two cases are refused: activating a wall or interface in a consolidation, fully-coupled or dynamic
phase (the stress takeover runs in Plastic phases — activate it in one first), and deactivating a
joint whose two sides have moved apart while a structure standing on its line is carried on.

### An undrained phase read the volume changes of earlier stages as pore pressure

An Undrained (A) or (B) material carries the excess pore pressure its undrained phases generate, and
the phases after them are supposed to find that pressure where it was left. They did not: the
pressure was re-derived in every phase as Kw/n times the volume change **since the initial state**,
whatever had produced it — a drained stage, a phase that ignores undrained behaviour, or a
consolidation. And a phase that ignores undrained behaviour deleted the pressure generated before it.
Measured on a laterally confined 1 × 8 m column (E = 5 MPa, ν = 0.3, water at the surface, a 25 kPa
surcharge; displacements of the column top in the phase):

| chain | 0.9.0 and this tree until now | now |
|---|---|---|
| load, ignoring undrained behaviour → an undrained phase that changes nothing | **28.7 mm of heave**, σ′yy at mid-depth −61.8 → −37.6 kPa | 0, unchanged |
| gravity loading ignoring undrained behaviour → 25 kPa undrained | **41.1 mm of heave** under a downward load, σ′yy −36.8 → −2.1 kPa | 1.030 mm down, −0.866 kPa (q H / M_u, −q M′ / M_u) |
| undrained load → a phase ignoring undrained behaviour that changes nothing | 28.7 mm down: the pressure was deleted and the column consolidated in a phase with no time | 0 |
| undrained load → consolidation → an undrained phase that changes nothing, linear-elastic ground | **28.7 mm of heave**: the consolidation undone | 0 |
| the same on Mohr-Coulomb ground | 0.99 mm of heave | 0 |

Gravity loading with undrained behaviour ignored, followed by undrained loading, is the ordinary way
to start an undrained analysis; the heave above is what it produced.

**The excess pore pressure is now a state of each stress point**, carried from phase to phase. An
undrained phase adds to it; a phase that ignores undrained behaviour neither adds to it nor removes
it, and the water adds no stiffness there; a consolidation or fully coupled phase writes the pressure
it ends with back to the stress points (in the fully coupled phase, the share the skeleton feels,
S_eff times the pressure), so the Plastic phase after it starts from that pressure. A Plastic phase
that inherits a pressure starts from its parent's equilibrium with the pressure included, and ramps
only its own change. In a chain in which every phase generates pressure the pressure is computed by
the same expression as before, and every existing undrained and consolidation case passes unchanged.

The same switch had a second effect. Ignoring undrained behaviour cleared the flag that marks a
material undrained, and a design phase reads that flag to give an Undrained (B) or (C) strength the
undrained partial factor γ_cu. **In an EC7 DA1-C2 or DA3 phase that ignored undrained behaviour, an
Undrained (B) soil was factored with γ_c′ instead** (1.25 in place of 1.4). On the checked-in strip
footing on a Tresca clay under DA3, where the water cannot change the capacity, the collapse load
factor was 0.52500 with undrained behaviour ignored and 0.47187 without it; both are now 0.471875.

`KV-CST-017` checks nine chains on the confined column against the closed forms: a phase that
changes nothing moves nothing (0 m, where the defects measured 2.9·10⁻² m and 9.9·10⁻⁴ m); an
undrained load after a drained stage, after gravity loading, after a consolidation
and after a fully coupled phase settles q H / M_u (to 3·10⁻¹⁴); a load while ignoring undrained
behaviour, on top of an undrained one, settles q H / M′ (8·10⁻¹⁶); a consolidation stopped part-way
leaves its remainder to the next consolidation, which reaches the drained settlement; and the footing
above gives the same collapse factor either way. `K2D-A008` now says that the pressure generated
earlier is kept.

### Upgrading from 0.9.0

- Results change for a model with an Undrained (A) or (B) material wherever an undrained phase
  follows a phase that ignores undrained behaviour, a consolidation or a fully coupled phase, and in
  any phase that ignores undrained behaviour after an undrained one. A design phase that ignores
  undrained behaviour now factors an Undrained (B) strength with γ_cu.
- A wall with interfaces, or an interface, may be inactive in a phase; files that kept one active in
  every phase are unchanged.
- A staged phase that activates or removes a structure now carries the state of every other
  structure; results of such chains change, the carried structures' forces most.
- A design phase (EC7 DA1-C2 or DA3) now factors interfaces, the cohesion gradient, Soft Soil and
  Soft Soil Creep. Its results change wherever a model has any of them, in the safe direction.
- An interface beside an Undrained (B) or (C) Mohr-Coulomb material no longer takes a friction
  angle from that material's φ box. Results change wherever that box is not 0.
- A Safety phase (or initial Safety procedure) now solves its active structural elements, so its
  factor of safety changes wherever one is active: 0.9.0 reported the factor of the same model
  without them. A prestressed anchor active in a Safety run is refused (`K2D-G016`); deactivate it
  in the Safety phase. A Safety run with structures raises the note `K2D-A018`.
- Scripts that matched on `K2D-A003` will no longer see it. An anchor whose end stands on a
  non-zero prescribed displacement now reports its force; before, it reported 0.
- Two corpus files are renamed, `tests/corpus/kv-str-002-sliding-block.k2d` and
  `tests/corpus/kv-str-003-beam-bending.k2d`, and so are the examples the installer ships from them.
  The test target that runs the four footing benchmarks is now `test_foundation_benchmarks`, and
  the circular-footing study `study_circular_footing`. Contents and assertions are as described
  above.
- The validation comparison page for the four footing benchmarks is withdrawn; each case is in
  `docs/validation/verification-matrix.md` against its analytical solution.
- Results files are written as `.res` version 10, which 0.9.0 refuses as "from a newer version".
  This build reads the older files; their stopping ratio shows as not recorded.
- The command line prints one more line per converged phase, and scripts that read the number after
  `converged: force error` now get the ratio the step stopped on (before: the stiffness-weighted
  one, still printed on the line below). In Python, `Convergence.force_error` is unchanged; the
  stopping ratio is `global_error`.
- Consolidation and fully-coupled results change for a model whose materials have different
  `e_init`, most for stiff layers (see the table above); models with one `e_init` are unchanged.

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
