"""Staged construction with KATAI 2D -- deactivating soil, phase by phase.

The problem: a 10 m two-layer profile; the upper 4 m stratum is excavated in
one staged phase across the full width. Because the column stays laterally
confined, the pit floor's rebound has a closed form -- one-dimensional
elastic unloading with the constrained (oedometer) modulus:

    heave u = + gamma_f h_exc H_rem / E_oed,
    E_oed   = E (1 - nu) / ((1 + nu)(1 - 2 nu)).

Units everywhere: kN, m, day. Run with the built package on PYTHONPATH::

    python staged_excavation.py

The same case is pinned in the verification corpus (KV-EXC-001): the FEM
heave matches the closed form within 2% on this mesh (measured +0.0%).
"""
import katai

E, NU = 1.0e4, 0.3
E_OED = E * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU))   # constrained modulus
GAMMA_F, H_EXC, H_REM = 17.0, 4.0, 6.0                     # fill weight, cut, remaining
HEAVE = GAMMA_F * H_EXC * H_REM / E_OED                    # closed-form rebound [m]

# 1. Two stacked strata. Only the region handle of the upper one matters later:
#    staged construction works by (de)activating regions by handle.
prj = katai.Project("Staged excavation", mesh_size=1.5, auto_refine=False)
lower_m = prj.materials.linear_elastic("Lower stratum", E=E, nu=NU,
                                       gamma=18.0, gamma_sat=20.0)
upper_m = prj.materials.linear_elastic("Upper stratum (excavated)", E=E, nu=NU,
                                       gamma=GAMMA_F, gamma_sat=20.0)
prj.geometry.rectangle(0.0, 0.0, 20.0, 6.0, material=lower_m, name="Lower")
upper = prj.geometry.rectangle(0.0, 6.0, 20.0, 10.0, material=upper_m,
                               name="Upper", bottom="free")

# 2. The initial state is the FULL geometry under K0; the staged phase then
#    removes the upper stratum. Everything a phase does not mention is
#    inherited -- deactivate lists only what changes, like the GUI's files.
prj.phases.plastic("Excavate", deactivate=[upper])

# 3. Run and read the excavation phase's displacement field.
job = prj.run(on_phase=lambda i, n, name: print(f"  phase {i + 1}/{n}: {name}"))
res = job.results()[-1]

# 4. The pit floor is y = 6; its rebound is the largest displacement in the
#    model (pure unloading, everything moves up).
print(f"pit-floor heave (FEM):        {res.max_disp * 1000:7.2f} mm")
print(f"1D elastic unloading (hand):  {HEAVE * 1000:7.2f} mm")
print(f"difference:                   {(res.max_disp / HEAVE - 1) * 100:+7.2f} %")

# The corpus band for this exact case, on this exact mesh.
if abs(res.max_disp / HEAVE - 1.0) > 0.02:
    raise SystemExit("excavation heave drifted outside the corpus band -- investigate")
