"""Resetting the displacements resets the small-strain history.

A strip load is placed on a sand with small-strain stiffness and then a second, larger one. The
second phase is run three ways:

  * plainly -- the sand meets it with the stiffness the first load degraded it to;
  * with ``reset_small_strain=True`` -- it meets it at G0;
  * with ``reset_displacements=True`` -- a new event read on its own, which has to meet the sand
    at G0 as well: the strain the earlier phases left is not the history that event meets.

So the third settles exactly as the second, both settle less than the first, and the phase says
what it did (K2D-M005). A soil without small-strain stiffness is indifferent to all three.
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def run(small, **second):
    prj = katai.Project("small strain reset", element="tri6", mesh_size=1.0, auto_refine=False)
    if small:
        sand = prj.materials.hs_small("Sand", E50=25e3, Eoed=25e3, Eur=75e3, G0=80e3, gamma07=0.2e-3,
                                      c=1.0, phi=31, psi=1, m=0.5, nu_ur=0.2, p_ref=100, Rf=0.9,
                                      gamma=17, gamma_sat=20)
    else:
        sand = prj.materials.hardening_soil("Sand", E50=25e3, Eoed=25e3, Eur=75e3, c=1.0, phi=31,
                                            psi=1, m=0.5, nu_ur=0.2, p_ref=100, Rf=0.9, gamma=17,
                                            gamma_sat=20)
    prj.geometry.polygon([(0, 0), (8, 0), (8, 6), (0, 6)], material=sand)
    q1 = prj.loads.line_load((0, 6), (2, 6), qy=-30.0, name="q1")
    q2 = prj.loads.line_load((0, 6), (2, 6), qy=-40.0, name="q2")
    prj.initial(procedure="k0", exclude=[q1, q2])
    prj.phases.plastic("first", activate=[q1])
    prj.phases.plastic("second", activate=[q2], **second)
    job = katai.Job(prj.build())
    ok = job.run()
    return ok, job


def settlement(job):
    r = job.results()[-1]
    i = min(range(len(r.node_x)), key=lambda k: r.node_x[k] ** 2 + (r.node_y[k] - 6.0) ** 2)
    return -r.displacement[2 * i + 1]


def codes(job):
    return [d.code for d in job.results()[-1].diagnostics]


ok_a, a = run(True)
ok_b, b = run(True, reset_small_strain=True)
ok_c, c = run(True, reset_displacements=True)
check(ok_a and ok_b and ok_c, "the three second phases run"
      + ("" if ok_a and ok_b and ok_c else " -- " + (a.message() + b.message() + c.message())[:300]))
if ok_a and ok_b and ok_c:
    sa, sb, sc = settlement(a), settlement(b), settlement(c)
    check(sa > 0 and sb < 0.97 * sa,
          f"reset small strain: the second load settles {1000 * sb:.3f} mm against {1000 * sa:.3f} mm")
    check(abs(sc / sb - 1.0) < 1e-9,
          f"reset displacements settles exactly as much: {1000 * sc:.6f} and {1000 * sb:.6f} mm")
    check("K2D-M005" in codes(b) and "K2D-M005" in codes(c) and "K2D-M005" not in codes(a),
          "both resetting phases say that the history was reset (K2D-M005), the plain one does not")

ok_d, d = run(False)
ok_e, e = run(False, reset_displacements=True)
check(ok_d and ok_e, "the same two phases run on a soil without small-strain stiffness")
if ok_d and ok_e:
    sd, se = settlement(d), settlement(e)
    check(abs(se / sd - 1.0) < 1e-12 and "K2D-M005" not in codes(e),
          f"and there resetting the displacements changes nothing and says nothing: "
          f"{1000 * sd:.6f} and {1000 * se:.6f} mm")

sys.exit(1 if failures else 0)
