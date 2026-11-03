"""A consolidation step is closed by its second correction, not its first.

A layer of Modified Cam-Clay consolidates under a load (the case of Chen, Yu, Abousleiman & Kees 2024,
arXiv:2410.04623, Figs 2-3: 20 m doubly drained, half modelled, OCR 2, 200 kPa on 50 kPa). In one-dimensional
consolidation the effective stress and the pore pressure add up to the total stress at every depth and every time, so
that sum is an answer known without a solver. A first correction that happened to land just under the force
tolerance used to close its step: the excess pore pressure dropped 2.74 kPa in that step against 0.44 in its
neighbours, the sum stayed at 247.73 kPa from there on instead of 250, and the excess pore pressure at mid-depth
after 100 days came out 3.1% low (67.78 kPa against 69.96 from an independent finite-volume solve), the
settlement hardly moving.
"""
import math
import sys

import numpy as np

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


LAM, KAP, M, NU = 0.15, 0.03, 1.2, 0.278
V0 = 2.258
PHI = math.degrees(math.asin(3 * M / (6 + M)))
H = 10.0
SV, Q0 = 50.0, 200.0
K = 1.96e-8 / 9.8 * 86400.0 * 9.81
SCHED = [(1.0, 20), (5.0, 40), (10.0, 25), (25.0, 30), (50.0, 25), (75.0, 25), (100.0, 25)]

prj = katai.Project("consolidation balance", element="tri15", mesh_size=0.5, auto_refine=False)
clay = prj.materials.modified_cam_clay("MCC", lam_star=LAM / V0, kap_star=KAP / V0, phi_cs=PHI, nu_ur=NU,
                                       gamma=0.0, gamma_sat=0.0, OCR=2.0, k0=1.0, k0nc=1.0, k=K,
                                       drainage="undrained_a")
dead = prj.materials.linear_elastic("Dead", E=1e5, nu=0.3, gamma=SV)
prj.geometry.polygon([(0, 0), (1, 0), (1, H), (0, H)], material=clay, name="Layer",
                     fix=["vertical", "horizontal", "free", "horizontal"],
                     flow=["closed", "closed", ("head", H), "closed"])
block = prj.geometry.polygon([(0, H), (1, H), (1, H + 1), (0, H + 1)], material=dead, name="Block",
                             fix=["free", "horizontal", "free", "horizontal"])
top = prj.loads.line_load((0, H), (1, H), qy=-SV, name="top")
q = prj.loads.line_load((0, H), (1, H), qy=-Q0, name="q0")
drain = prj.dewatering.drain((0, H), (1, H), head=H, name="drain")
prj.phases.plastic("swap", activate=[top], deactivate=[block], ignore_undrained=True)
t = 0.0
for i, (tend, steps) in enumerate(SCHED):
    kw = dict(activate=[q, drain], load_at_start=True) if i == 0 else {}
    prj.phases.consolidation(f"c{i}", duration=tend - t, steps=steps, **kw)
    t = tend
prj.initial(procedure="k0", exclude=[top, q, drain])
job = katai.Job(prj.build())
check(job.run(), "the column consolidates to 100 days")
res = job.results()
last = res[-1]
y = np.array(last.node_y)
inside = y <= H + 1e-9
s = np.array(last.stress)
ep = np.array(last.excess_pore)
total = -s[inside, 1] - ep[inside]
u_mid = -ep[inside][np.argmin(y[inside])]
print(f"      u at mid-depth {u_mid:.3f} kPa (independent 69.96), sum {total.min():.3f} .. {total.max():.3f} kPa")
check(abs(u_mid / 69.96 - 1.0) < 0.015, "the excess pore pressure at mid-depth, within 1.5% of the independent solve")
check(total.min() > 249.0 and total.max() < 251.0, "effective stress and pore pressure add up to the 250 kPa, within 1 kPa")
sys.exit(1 if failures else 0)
