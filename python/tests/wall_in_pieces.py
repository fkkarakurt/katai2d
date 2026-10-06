"""A wall drawn in pieces is the wall drawn in one.

A wall with interfaces moves on its own degrees of freedom and hangs from its toe, the one node
it shares with the ground. A wall drawn in two pieces -- at a layer boundary, at the excavation
level -- was built as two walls, and the upper one took the point where they meet as ITS toe:
the wall was tied to the soil on one side at every point it had been drawn through. On an 8 m
cantilever wall in front of a 4 m excavation the largest moment was 20.5 kNm/m in one piece and
18.1 in two, and the head moved 11.8 against 9.7 mm, with nothing to say so.

The cut is put on a node the mesh has anyway (the excavation level), so both runs are solved on
the same mesh and the two answers must agree to the solver's tolerance, not to a mesh's.
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def wall(cuts):
    prj = katai.Project("wall in pieces", element="tri6", mesh_size=1.0, auto_refine=False)
    soil = prj.materials.mohr_coulomb("Soil", E=20000, nu=0.3, c=5.0, phi=30.0, psi=0.0,
                                      gamma=18.0, rinter_rigid=False, Rinter=0.6)
    g = prj.geometry
    g.polygon([(0, 0), (20, 0), (20, 6), (10, 6), (10, 10), (0, 10)], material=soil,
              name="Retained", fix=["full", "horizontal", "free", "free", "free", "horizontal"])
    pit = g.polygon([(10, 6), (20, 6), (20, 10), (10, 10)], material=soil, name="Pit",
                    fix=["free", "horizontal", "free", "free"])
    ys = [10.0] + list(cuts) + [2.0]
    pieces = [prj.structures.plate((10, ys[i]), (10, ys[i + 1]), EA=5e6, EI=5e4, w=0.0, nu=0.0,
                                   interfaces="both", name=f"W{i}") for i in range(len(ys) - 1)]
    prj.initial(procedure="k0", exclude=pieces)
    prj.phases.plastic("Wall", activate=pieces)
    prj.phases.plastic("Excavate", deactivate=[pit])
    job = katai.Job(prj.build())
    ok = job.run()
    if not ok:
        return False, job.message(), 0, 0.0, 0.0
    r = job.results()[-1]
    st = [q for s in r.struct_forces for q in s.stations]
    return (True, "", len(r.node_x), max(abs(q.M) for q in st),
            1000.0 * max(abs(q.ux) for q in st))


ok1, msg1, n1, m1, u1 = wall([])
ok2, msg2, n2, m2, u2 = wall([6.0])
check(ok1 and ok2, "the wall runs drawn in one piece and in two" +
      ("" if ok1 and ok2 else " -- " + (msg1 or msg2)[:200]))
check(ok1 and ok2 and n1 == n2, f"on the same mesh ({n1} and {n2} nodes)")
check(ok1 and ok2 and abs(m2 / m1 - 1) < 5e-3,
      f"with the same largest moment: {m1:.2f} and {m2:.2f} kNm/m ({100 * (m2 / m1 - 1):+.3f}%; "
      f"was -12%)")
check(ok1 and ok2 and abs(u2 / u1 - 1) < 5e-3,
      f"and the same head movement: {u1:.3f} and {u2:.3f} mm ({100 * (u2 / u1 - 1):+.3f}%; "
      f"was -18%)")


# A wall through two soils. Its interfaces used to be those of the soil at the wall's midpoint,
# over its whole length and on both faces; drawn in two pieces cut at the layer boundary, each
# piece took its own soil. Each interface element now takes the soil it lies against, so the two
# drawings are the same model. The layers differ in strength and in Rinter, so they differ in
# what the interface is given.
def layered(cuts, e_box=None):
    prj = katai.Project("wall through two soils", element="tri6", mesh_size=1.0, auto_refine=False)
    kw = {} if e_box is None else {"E": e_box}
    clay = prj.materials.hardening_soil("Clay", E50=6000, Eoed=6000, Eur=18000, m=0.8, c=8.0,
                                        phi=24.0, psi=0.0, gamma=17.0, nu_ur=0.2, p_ref=100.0,
                                        Rf=0.9, rinter_rigid=False, Rinter=0.5, **kw)
    sand = prj.materials.mohr_coulomb("Sand", E=40000, nu=0.3, c=1.0, phi=34.0, psi=4.0,
                                      gamma=19.0, rinter_rigid=False, Rinter=0.7)
    g = prj.geometry
    g.polygon([(0, 0), (20, 0), (20, 5), (10, 5), (0, 5)], material=sand, name="Sand",
              fix=["full", "horizontal", "free", "free", "horizontal"])
    g.polygon([(0, 5), (10, 5), (10, 7), (10, 10), (0, 10)], material=clay, name="Clay",
              fix=["free", "free", "free", "free", "horizontal"])
    low = g.polygon([(10, 5), (20, 5), (20, 7), (10, 7)], material=clay, name="Clay in front",
                    fix=["free", "horizontal", "free", "free"])
    pit = g.polygon([(10, 7), (20, 7), (20, 10), (10, 10)], material=clay, name="Pit",
                    fix=["free", "horizontal", "free", "free"])
    ys = [10.0] + list(cuts) + [2.0]
    pieces = [prj.structures.plate((10, ys[i]), (10, ys[i + 1]), EA=5e6, EI=5e4, w=0.0, nu=0.0,
                                   interfaces="both", name=f"W{i}") for i in range(len(ys) - 1)]
    prj.initial(procedure="k0", exclude=pieces)
    prj.phases.plastic("Wall", activate=pieces)
    prj.phases.plastic("Excavate", deactivate=[pit])
    job = katai.Job(prj.build())
    ok = job.run()
    if not ok:
        return False, job.message(), 0.0, 0.0
    st = [q for s in job.results()[-1].struct_forces for q in s.stations]
    return True, "", max(abs(q.M) for q in st), 1000.0 * max(abs(q.ux) for q in st)


oka, msga, ma, ua = layered([])
okb, msgb, mb, ub = layered([5.0])
check(oka and okb, "a wall through clay into sand runs drawn in one piece and cut at the boundary"
      + ("" if oka and okb else " -- " + (msga or msgb)[:200]))
# The clay is a hardening soil, solved to its tolerated error of 3e-3: the two drawings are the
# same model with their unknowns in a different order, so they take different iteration paths
# and agree to about that, not to the last digit (0.19% measured; it was 2.8% while the
# tolerated error was 1e-2, which is what moved the default).
check(oka and okb and abs(mb / ma - 1) < 5e-3 and abs(ub / ua - 1) < 5e-3,
      f"and is the same wall: {ma:.2f} and {mb:.2f} kNm/m, {ua:.3f} and {ub:.3f} mm "
      f"({100 * (mb / ma - 1):+.3f}%, {100 * (ub / ua - 1):+.3f}%)")

# The stiffness of an interface beside a hardening soil is that soil's own (E50). It was read
# from the Young's modulus box, which a hardening material never uses: the answer moved with a
# number the model does not have.
okc, msgc, mc, uc = layered([], e_box=130000.0)
check(okc and oka and abs(mc / ma - 1) < 1e-9 and abs(uc / ua - 1) < 1e-9,
      f"the unused Young's modulus box of a hardening soil changes nothing: {ma:.4f} and "
      f"{mc:.4f} kNm/m, {ua:.4f} and {uc:.4f} mm" + ("" if okc else " -- " + msgc[:200]))

sys.exit(1 if failures else 0)
