"""A plate that is given another section in a phase keeps the forces it carries.

A plate lies on a weightless linear-elastic block and is loaded in two steps. In run A it has
the section S1 under the first load, is given S2 (eight times stiffer in bending, four times
axially) in a phase where nothing else happens, and takes the second load with S2. Everything is
linear, so the state at the end is a superposition with a closed statement:

    forces(A, end) = forces(plate S1, first load) + forces(plate S2 from the start, second load)

and in the phase of the change itself nothing may move: the moments and the settlements are
those of the phase before. Multiplying the carried forces by the stiffness ratio -- what a plate
evaluated on its total deformation does when its stiffness is simply replaced -- fails both.
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


S1 = dict(EA=2.0e6, EI=2.0e4)
S2 = dict(EA=8.0e6, EI=1.6e5)
Q1, Q2 = 100.0, 60.0


def run(section, loads, change=None):
    prj = katai.Project("plate section", element="tri6", mesh_size=0.5, auto_refine=False)
    soil = prj.materials.linear_elastic("Soil", E=2.0e4, nu=0.3, gamma=0.0, gamma_sat=0.0)
    prj.geometry.polygon([(0, 0), (12, 0), (12, 5), (0, 5)], material=soil,
                         fix=["full", "horizontal", "free", "horizontal"])
    plate = prj.structures.plate((3, 5), (9, 5), w=0.0, nu=0.0, name="Slab", **section)
    qa = prj.loads.line_load((5, 5), (7, 5), qy=-Q1, name="first")
    qb = prj.loads.line_load((5, 5), (7, 5), qy=-Q2, name="second")
    prj.initial(procedure="k0", exclude=[qa, qb])
    if "a" in loads:
        prj.phases.plastic("first load", activate=[qa])
    if change:
        prj.phases.plastic("new section", plates={plate: change})
    if "b" in loads:
        prj.phases.plastic("second load", activate=[qb])
    job = katai.Job(prj.build())
    ok = job.run()
    return ok, job


def moments(res):
    return [q.M for s in res.struct_forces for q in s.stations]


def axial(res):
    return [q.N for s in res.struct_forces for q in s.stations]


ok_a, a = run(S1, "ab", change=S2)
ok_1, r1 = run(S1, "a")
ok_2, r2 = run(S2, "b")
check(ok_a and ok_1 and ok_2, "the three runs complete"
      + ("" if ok_a and ok_1 and ok_2 else " -- " + (a.message() + r1.message() + r2.message())[:300]))
if ok_a and ok_1 and ok_2:
    ra = a.results()
    m_first, m_change, m_end = moments(ra[1]), moments(ra[2]), moments(ra[3])
    scale = max(abs(m) for m in m_first)
    worst = max(abs(x - y) for x, y in zip(m_first, m_change)) / scale
    check(scale > 1.0 and worst < 1e-6,
          f"in the phase of the change the moments stay what they were: largest {scale:.3f} kNm/m, "
          f"changed by {worst:.2e} of it (replacing the stiffness alone would make them 8x)")
    moved = max(abs(x) for x in ra[2].displacement) / max(abs(x) for x in ra[1].displacement)
    check(moved < 1e-6, f"and nothing moves: the phase displacement is {moved:.2e} of the loading phase's")
    want = [x + y for x, y in zip(moments(r1.results()[1]), moments(r2.results()[1]))]
    worst = max(abs(x - y) for x, y in zip(m_end, want)) / max(abs(m) for m in want)
    check(worst < 1e-6,
          f"under the second load it is the old section's moments plus the new section's: "
          f"largest {max(abs(m) for m in m_end):.3f} against {max(abs(m) for m in want):.3f} kNm/m "
          f"(apart by {worst:.2e})")
    n_want = [x + y for x, y in zip(axial(r1.results()[1]), axial(r2.results()[1]))]
    n_scale = max(max(abs(n) for n in n_want), 1e-9)
    worst_n = max(abs(x - y) for x, y in zip(axial(ra[3]), n_want)) / n_scale
    check(worst_n < 1e-6, f"and the same for the normal force (apart by {worst_n:.2e} of {n_scale:.3f} kN/m)")
    stiff_only = max(abs(m) for m in moments(r2.results()[1]))
    check(abs(max(abs(m) for m in m_end) / stiff_only - 1) > 0.05,
          "which is not what the new section alone would carry under the second load")

sys.exit(1 if failures else 0)
