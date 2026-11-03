"""Planar sliding of a block on a jointed rock's plane: the factor does not depend on the dilatancy of the slip.

The block slides along a persistent set dipping out of the face, so its factor is the limit-equilibrium closed form of
plane failure (Hoek, E. & Bray, J. (1981). Rock Slope Engineering, 3rd ed., ch. 7: F = (c A + W cos(p) tan(phi)) / (W sin(p))),
and slip along the plane needs no volume change: psi does not enter it. With psi < phi the kinks of the planes used
to stop this model well short of that. The self-weight stopped at 83% with psi 0 and 86% with psi 5 (15-noded, 2 m), and a band of
the mesh softens long before any mechanism, so the strength reduction path dips: with psi 10 it reported 1.03 with nothing moving by
more than 18 micrometres, on an exact cycle of Newton iterates. The dip is passed now and the path rises to the factor.
"""
import math
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


H, FACE, DIP, GAMMA = 10.0, 60.0, 35.0, 25.0
E, NU = 2.0e6, 0.25
C, PHI = 20.0, 30.0
t, r = math.tan, math.radians
AREA = 0.5 * H * (H / t(r(DIP)) - H / t(r(FACE)))
W = GAMMA * AREA
A = H / math.sin(r(DIP))
F_CLOSED = (C * A + W * math.cos(r(DIP)) * t(r(PHI))) / (W * math.sin(r(DIP)))


def slope(psi, mesh=2.0):
    prj = katai.Project(f"block psi {psi}", element="tri15", mesh_size=mesh, corner_refine=False, auto_refine=False)
    rock = prj.materials.jointed_rock("Jointed rock", E=E, nu=NU, gamma=GAMMA,
                                      joints=[dict(dip=DIP, c=C, phi=PHI, psi=psi, sigma_t=0.0)], tensile_cutoff=True)
    tan = t(r(DIP))
    xt, yt = 15.0, 15.0
    xc = xt + H / t(r(FACE))
    xp = xt + H / tan
    x0, x1, y0 = 10.0, 33.0, 10.0
    yl = yt - (xt - x0) * tan
    base = prj.materials.linear_elastic("Intact rock", E=E, nu=NU, gamma=GAMMA)
    prj.geometry.polygon([(x0, y0), (x1, y0), (x1, yt + H), (xp, yt + H), (xt, yt), (x0, yl)], material=base, name="Base",
                         fix=["full", "horizontal", "free", "free", "free", "horizontal"])
    prj.geometry.polygon([(xt, yt), (xp, yt + H), (xc, yt + H)], material=rock, name="Block", fix=["free", "free", "free"])
    prj.initial(procedure="gravity")
    prj.phases.safety("FoS")
    job = katai.Job(prj.build())
    job.run()
    return job.results()


print(f"closed form {F_CLOSED:.4f}")
for psi in (30.0, 10.0, 5.0, 0.0):
    res = slope(psi)
    check(len(res) == 2 and res[0].ok, f"psi {psi:.0f}: the self-weight is carried")
    if len(res) == 2:
        fos = res[1].fos
        print(f"   psi {psi:4.0f}: factor {fos:.4f} ({100.0 * (fos / F_CLOSED - 1.0):+.2f}%)")
        check(res[1].ok and 0.0 <= fos / F_CLOSED - 1.0 < 0.025, f"psi {psi:.0f}: the factor of plane failure on this mesh (0 to +2.5%)")

print("\n" + ("OK" if not failures else f"{failures} check(s) failed"))
sys.exit(1 if failures else 0)
