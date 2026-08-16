"""An anchored excavation, entirely in Python -- no .k2d file anywhere.

This is the shape of a real job: ground, a diaphragm wall with interfaces on both
faces, a prestressed anchor row, and a staged excavation that installs the anchor
and digs the pit. Every object is created in one call and handed back as a handle,
and the phases activate and deactivate those handles -- which is the whole of
staged construction.

Units are fixed: kN, m, day. Run it with the built package on PYTHONPATH::

    python anchored_excavation.py
"""
import katai

prj = katai.Project("Anchored excavation", mesh_size=2.0, auto_refine=False)

# 1. One soil for the whole site.
sand = prj.materials.mohr_coulomb("Berlin sand", E=3.0e4, nu=0.3,
                                  c=5.0, phi=32.0, gamma=19.0)

# 2. Ground, split so the pit can be dug out of it later. Rectangles take the same
#    per-edge fixity names as polygons; the defaults (base fixed, sides on rollers,
#    surface free) are what a site usually wants.
prj.geometry.rectangle(0, 0, 30, 12, material=sand, name="Below formation")
pit = prj.geometry.rectangle(0, 12, 12, 20, material=sand, name="Pit")
prj.geometry.rectangle(12, 12, 30, 20, material=sand, name="Retained ground")

# 3. The wall: EA and EI per metre out of plane, an interface on each face so the
#    soil can slip against it instead of being glued to it.
wall = prj.structures.plate((12, 20), (12, 6), EA=1.2e7, EI=1.0e5, w=5.0,
                            interfaces="both", name="Diaphragm wall")

# 4. The anchor: a discrete member, so it carries its out-of-plane spacing, and a
#    lock-off force that holds the excavation before anything moves.
anchor = prj.structures.anchor((12, 18), (20, 14), EA=2.0e5, spacing=2.5,
                               prestress=300.0, name="Anchor row 1")

# 5. Staging. The wall is in place from the start -- an embedded wall cannot be
#    switched per phase, and the engine says so honestly if you try. The anchor is
#    installed and the pit removed in the excavation phase.
prj.initial(procedure="k0", exclude=[anchor])
prj.phases.plastic("Excavate to formation", activate=[anchor], deactivate=[pit])

# 6. Run. A refusal is an exception carrying the engine's own words and the stable
#    K2D-* codes, so a script can branch on the reason rather than parse prose.
try:
    job = prj.run(on_phase=lambda i, n, name: print(f"  [{i + 1}/{n}] {name}"))
except katai.Refusal as e:
    raise SystemExit(f"refused: {e}\ncodes: {[d.code for d in e.diagnostics]}")

# 7. Read it. summary() names the phases, tabulates the extremes WITH the place
#    each occurs, lists the structural force envelopes, and repeats any diagnostic
#    the engine raised -- because those change how the numbers should be read.
print()
print(katai.summary(job, prj))

# Anything can also be read directly; the summary is a convenience, not a gate.
last = job.results()[-1]
wall_force = next(s for s in last.struct_forces if s.name == "Diaphragm wall")
print(f"\nwall: |M|max = {wall_force.max_M:.1f} kNm/m, "
      f"|N|max = {wall_force.max_N:.1f} kN/m"
      f"{'  (YIELDED)' if wall_force.yielded else ''}")
print(f"largest movement anywhere: {last.max_disp * 1000:.1f} mm")

# And the same analysis can be handed to the GUI or the command line at any point:
#   prj.save("anchored_excavation.k2d")   ->   katai solve anchored_excavation.k2d

# --- the example checks itself, like every other one in this folder ------------------
# There is no published benchmark for this exact geometry, so what is pinned is what
# must be true of it whatever the numbers come out as. An example that quietly stops
# meaning what it says is worse than no example.
first, last_r = job.results()[0], job.results()[-1]
anchor_force = next(s for s in last_r.struct_forces if s.name == "Anchor row 1")
lock_off = 300.0 / 2.5          # [kN/m] -- a discrete member divided by its spacing
problems = []
if not all(r.ok for r in job.results()):
    problems.append("a phase did not converge")
if not last_r.max_disp > first.max_disp:
    problems.append("excavating moved the ground LESS than the in-situ phase")
if not anchor_force.max_N >= lock_off - 1e-6:
    problems.append(f"anchor force {anchor_force.max_N:.1f} fell below its "
                    f"{lock_off:.0f} kN/m lock-off, which a tensioned member cannot do")
if not wall_force.max_M > 0.0:
    problems.append("the wall carries no moment, so it is not retaining anything")
if problems:
    raise SystemExit("example no longer demonstrates what it claims:\n  - "
                     + "\n  - ".join(problems))
