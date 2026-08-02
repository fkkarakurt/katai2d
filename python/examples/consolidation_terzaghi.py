"""One-dimensional consolidation with KATAI 2D -- reading a time series.

The problem: Terzaghi's classic column (Terzaghi 1943). A 12 m clay layer,
drained at the ground surface only, takes a 10 kPa surcharge at t = 0+ and
consolidates: the excess pore pressure the load created bleeds off through
the top and the surface settles with time. The closed form is the series

    U(Tv) = 1 - sum_j (2/M^2) exp(-M^2 Tv),   M = (2j+1) pi/2,

with Tv = cv t / H^2 and cv = k E_oed / gamma_w.

Units everywhere: kN, m, day. Run with the built package on PYTHONPATH::

    python consolidation_terzaghi.py

The same case is pinned in the verification corpus (KV-CON-002): the FEM
degree of consolidation stays within 3% of the series on this mesh.
"""
import math

import katai

E_OED = 1000.0    # oedometer modulus [kN/m2] (nu = 0 makes E the constrained modulus)
K = 0.1           # permeability [m/day]
GAMMA_W = 9.81    # water unit weight the engine uses [kN/m3]
H = 12.0          # drainage path [m]: drained at the top ONLY, so the full height
Q = 10.0          # surcharge [kN/m2]
CV = K * E_OED / GAMMA_W          # consolidation coefficient [m2/day]
S_INF = Q * H / E_OED             # final (fully drained) settlement [m]

# 1. The column: 1 m wide, 12 m tall, fine mesh so the early pore front resolves.
prj = katai.Project("Terzaghi column", mesh_size=0.4, auto_refine=False)

# 2. One clay; nu = 0 so E is exactly the oedometer modulus of the closed form.
clay = prj.materials.linear_elastic("Oedometer clay", E=E_OED, nu=0.0,
                                    gamma=16.0, gamma_sat=18.0, k=K)

# 3. Geometry with FLOW boundary conditions: water may leave through the top
#    (its head stays at the surface, 12 m); every other face is sealed.
prj.geometry.rectangle(0.0, 0.0, 1.0, H, material=clay, name="Column",
                       flow={"top": ("head", H), "bottom": "closed",
                             "right": "closed", "left": "closed"})

# 4. The surcharge exists as a load object but is EXCLUDED from the initial
#    state -- it must arrive at t = 0+, undrained, inside the consolidation
#    phase, or there would be no excess pore pressure to dissipate.
q = prj.loads.line_load((0.0, H), (1.0, H), qy=-Q, name="Surcharge")
prj.initial(procedure="k0", exclude=[q])

# 5. Consolidate to Tv ~ 2 (practically complete) in 120 equal time steps.
prj.phases.consolidation("Consolidation", duration=2.0 * H * H / CV, steps=120,
                         activate=[q])

job = prj.run()
res = job.results()[-1]

# 6. The time series: settlement of the drained surface at every time step.
times = list(res.consol_time)              # [day]
settle = list(res.consol_settlement)       # [m]


def terzaghi_U(tv, terms=200):
    return 1.0 - sum((2.0 / m_ ** 2) * math.exp(-m_ ** 2 * tv)
                     for m_ in ((2 * j + 1) * math.pi / 2.0 for j in range(terms)))


print("   t [day]     Tv     U (FEM)   U (Terzaghi)")
worst = 0.0
for tv_target in (0.2, 0.4, 0.6, 0.9):
    t_target = tv_target * H * H / CV
    i = min(range(len(times)), key=lambda n: abs(times[n] - t_target))
    tv = CV * times[i] / (H * H)
    u_fem, u_ref = settle[i] / S_INF, terzaghi_U(CV * times[i] / (H * H))
    worst = max(worst, abs(u_fem - u_ref))
    print(f"   {times[i]:7.2f}   {tv:5.3f}   {u_fem:.4f}    {u_ref:.4f}")

print(f"final settlement: {settle[-1]:.4f} m   (closed form s_inf = {S_INF:.4f} m)")
print(f"largest |U_FEM - U_Terzaghi| at the sampled times: {worst:.4f}")

# The corpus band for this exact case, on this exact mesh.
if worst > 0.03:
    raise SystemExit("consolidation drifted outside the corpus band -- investigate")
