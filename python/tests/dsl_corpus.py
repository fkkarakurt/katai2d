"""The engineer-facing DSL's acceptance test: every corpus case, BYTE-identical.

Each checked-in .k2d was built by a C++ programmatic builder and is pinned
byte-identical to it (test_input_corpus). Here the same seven cases are built
through `katai.Project` -- the easy surface -- and their canonical JSON must
equal the checked-in file BYTE FOR BYTE. One contract, three authors: the C++
builder, the file, and now the script. A mismatch prints where the two texts
diverge, which names the field the DSL got wrong.

End to end, one case also RUNS: the DSL-built Griffiths & Lane slope must give
exactly the same factor of safety as the file-loaded run in the same process.
"""
import sys

import katai

CORPUS = sys.argv[1]
GAMMA_W = 9.81  # the engine's water unit weight [kN/m3] (analysis constants)

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def expect_bytes(prj, filename):
    got = katai.project_to_json(prj.build())
    with open(f"{CORPUS}/{filename}", encoding="utf-8") as fh:
        want = fh.read()
    if got == want:
        check(True, f"{filename}: DSL == checked-in file, byte for byte")
        return
    n = next((i for i, (a, b) in enumerate(zip(got, want)) if a != b),
             min(len(got), len(want)))
    check(False, f"{filename}: first divergence at byte {n}: "
                 f"DSL ...{got[max(0, n - 30):n + 30]}... vs "
                 f"file ...{want[max(0, n - 30):n + 30]}...")


# ---------------------------------------------------------------- KV-CON-002 --
prj = katai.Project("KV-CON-002 Terzaghi column", mesh_size=0.4, auto_refine=False)
clay = prj.materials.linear_elastic("Oedometer clay", E=1000.0, nu=0.0,
                                    gamma=16.0, gamma_sat=18.0, k=0.1)
prj.geometry.rectangle(0.0, 0.0, 1.0, 12.0, material=clay, name="Column",
                       flow={"top": ("head", 12.0), "bottom": "closed",
                             "right": "closed", "left": "closed"})
q = prj.loads.line_load((0.0, 12.0), (1.0, 12.0), qy=-10.0, name="Surcharge")
prj.initial(procedure="k0", exclude=[q])
cv = 0.1 * 1000.0 / GAMMA_W
prj.phases.consolidation("Consolidation", duration=2.0 * 12.0 * 12.0 / cv,
                         steps=120, activate=[q])
expect_bytes(prj, "kv-con-002-terzaghi-column.k2d")

# ---------------------------------------------------------------- KV-FND-008 --
prj = katai.Project("KV-FND-008 strip load", mesh_size=1.0)
soil = prj.materials.linear_elastic("Weightless elastic", E=30000.0, nu=0.3,
                                    gamma=0.0)
prj.geometry.rectangle(0.0, 0.0, 40.0, 20.0, material=soil, name="Half-plane")
prj.loads.line_load((18.0, 20.0), (22.0, 20.0), qy=-100.0, name="Strip q")
prj.initial(procedure="gravity")
expect_bytes(prj, "kv-fnd-008-strip-load.k2d")

# ---------------------------------------------------------------- KV-FND-009 --
prj = katai.Project("KV-FND-009 Flamant line load", element="tri15", mesh_size=1.0)
soil = prj.materials.linear_elastic("Weightless elastic", E=30000.0, nu=0.3,
                                    gamma=0.0)
prj.geometry.rectangle(0.0, 0.0, 40.0, 20.0, material=soil, name="Half-plane")
prj.loads.point_load((20.0, 20.0), qy=-100.0, name="Line load P")
prj.initial(procedure="gravity")
expect_bytes(prj, "kv-fnd-009-flamant-line-load.k2d")

# ---------------------------------------------------------------- KV-NUM-003 --
prj = katai.Project("KV-NUM-003 K0 geostatic block", mesh_size=1.0,
                    auto_refine=False)
sand = prj.materials.linear_elastic("Sand", E=1.0e4, nu=0.3, gamma=17.0,
                                    gamma_sat=20.0, phi=30.0, c=1.0)
prj.geometry.rectangle(0.0, 0.0, 20.0, 10.0, material=sand, name="Block")
prj.water.table(10.0)
expect_bytes(prj, "kv-num-003-k0-geostatic-block.k2d")

# ---------------------------------------------------------------- KV-CST-001 --
prj = katai.Project("KV-CST-001 undrained column", mesh_size=1.0,
                    auto_refine=False)
clay = prj.materials.linear_elastic("Weightless undrained clay", E=1.0e4, nu=0.3,
                                    gamma=0.0, drainage="undrained_a")
prj.geometry.rectangle(0.0, 0.0, 2.0, 10.0, material=clay, name="Column",
                       right="normal", left="normal")
prj.loads.line_load((0.0, 10.0), (2.0, 10.0), qy=-50.0, name="Surcharge")
prj.initial(procedure="gravity")
expect_bytes(prj, "kv-cst-001-undrained-column.k2d")

# ---------------------------------------------------------------- KV-SLP-001 --
def build_slope():
    prj = katai.Project("KV-SLP-001 Griffiths-Lane slope", mesh_size=3.0,
                        auto_refine=False)
    soil = prj.materials.mohr_coulomb("Griffiths-Lane soil", E=1.0e5, nu=0.3,
                                      c=3.0, phi=19.6, gamma=20.2, gamma_sat=20.0)
    prj.geometry.polygon(
        [(20.0, 20.0), (70.0, 20.0), (70.0, 35.0), (50.0, 35.0),
         (30.0, 25.0), (20.0, 25.0)],
        material=soil, name="Slope",
        fix=["full", "horizontal", "free", "free", "free", "horizontal"])
    prj.initial(procedure="safety")
    return prj


expect_bytes(build_slope(), "kv-slp-001-griffiths-lane-slope.k2d")

# ---------------------------------------------------------------- KV-EXC-001 --
prj = katai.Project("KV-EXC-001 staged excavation", mesh_size=1.5,
                    auto_refine=False)
lower_m = prj.materials.linear_elastic("Lower stratum", E=1.0e4, nu=0.3,
                                       gamma=18.0, gamma_sat=20.0)
upper_m = prj.materials.linear_elastic("Upper stratum (excavated)", E=1.0e4,
                                       nu=0.3, gamma=17.0, gamma_sat=20.0)
prj.geometry.rectangle(0.0, 0.0, 20.0, 6.0, material=lower_m, name="Lower")
upper = prj.geometry.rectangle(0.0, 6.0, 20.0, 10.0, material=upper_m,
                               name="Upper", bottom="free")
prj.phases.plastic("Excavate", deactivate=[upper])
expect_bytes(prj, "kv-exc-001-staged-excavation.k2d")

# ---------------------------------------------------------------- KV-DYN-002 --
# f_1 is derived EXPRESSION FOR EXPRESSION as in the C++ builder (dyn_f1 in
# test_input_corpus.cpp): the file stores the derived numbers, and byte-identity
# across authors needs the same IEEE operations in the same order.
import math  # noqa: E402

DY_E, DY_NU, DY_GAMMA, DY_H = 208000.0, 0.3, 19.62, 20.0
_g = DY_E / (2.0 * (1.0 + DY_NU))
_rho = DY_GAMMA / 9.81
DY_F1 = math.sqrt(_g / _rho) / (4.0 * DY_H)

prj = katai.Project("KV-DYN-002 resonant shear column", mesh_size=0.4,
                    auto_refine=False)
soil = prj.materials.linear_elastic("Shear column soil", E=DY_E, nu=DY_NU,
                                    gamma=DY_GAMMA)
prj.geometry.rectangle(0.0, 0.0, 2.0, DY_H, material=soil, name="Column",
                       right="vertical", left="vertical")
prj.phases.dynamic("Resonance", duration=8.0, steps=800, frequency=DY_F1,
                   damping=0.05, rayleigh=(DY_F1, 3.0 * DY_F1))
prj.phases.dynamic("Off resonance", duration=9.0, steps=720,
                   frequency=DY_F1 / 3.0, damping=0.05,
                   rayleigh=(DY_F1, 3.0 * DY_F1))
expect_bytes(prj, "kv-dyn-002-resonant-column.k2d")

# ---------------------------------------------------------------- KV-FLW-001 --
FW_L, FW_D, FW_H1, FW_H2 = 10.0, 6.0, 5.0, 1.0
prj = katai.Project("KV-FLW-001 Charny unconfined dam", mesh_size=0.35,
                    auto_refine=False)
fill = prj.materials.linear_elastic("Dam fill", E=1.0e4, nu=0.3, gamma=18.0,
                                    gamma_sat=20.0, k=0.5)
prj.geometry.polygon(
    [(0.0, 0.0), (FW_L, 0.0), (FW_L, FW_H2), (FW_L, FW_D), (0.0, FW_D),
     (0.0, FW_H1)],
    material=fill, name="Dam",
    fix=["full", "free", "free", "free", "free", "free"],
    flow=["closed", ("head", FW_H2), "seepage", "closed", "closed",
          ("head", FW_H1)])
expect_bytes(prj, "kv-flw-001-charny-unconfined-dam.k2d")

# ---------------------------------------------------------------- KV-DYN-003 --
# Both authors read tests/data/elcentro-1940-ns.dat with the same parse and the
# same g -> m/s2 conversion (a * 9.81), so the file's seventeen-digit record
# bytes agree across authors.
import os  # noqa: E402

_times, _rec = [], []
with open(os.path.join(CORPUS, "..", "data", "elcentro-1940-ns.dat"),
          encoding="utf-8") as fh:
    for line in fh:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            tv, av = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        _times.append(tv)
        _rec.append(av * 9.81)
_dt = _times[1] - _times[0]

prj = katai.Project("KV-DYN-003 El Centro two-layer column", mesh_size=1.0,
                    auto_refine=False)
soft = prj.materials.linear_elastic("Soft upper layer", E=2.0 * 1.3 * 25920.0,
                                    nu=0.3, gamma=1.8 * 9.81)
stiff = prj.materials.linear_elastic("Stiff lower layer", E=2.0 * 1.3 * 189000.0,
                                     nu=0.3, gamma=2.1 * 9.81)
prj.geometry.rectangle(0.0, 12.0, 2.0, 20.0, material=soft, name="Upper",
                       bottom="free", right="vertical", left="vertical")
prj.geometry.rectangle(0.0, 0.0, 2.0, 12.0, material=stiff, name="Lower",
                       right="vertical", left="vertical")
prj.phases.dynamic("ElCentro", duration=(len(_rec) - 1) * _dt,
                   steps=min(len(_rec) - 1, 20000), wave="record", record=_rec,
                   record_dt=_dt, damping=0.05, rayleigh=(3.0, 9.0),
                   compliant_base=True)
expect_bytes(prj, "kv-dyn-003-el-centro-two-layer.k2d")

# ---------------------------------------------------------------- KV-FND-010 --
prj = katai.Project("KV-FND-010 Prandtl strip footing", element="tri15",
                    mesh_size=0.4, auto_refine=False)
clay = prj.materials.mohr_coulomb("Weightless Tresca clay", E=10000.0, nu=0.3,
                                  c=10.0, phi=0.0, psi=0.0, gamma=0.0,
                                  tension_cutoff=False)
prj.geometry.rectangle(0.0, 0.0, 6.0, 4.0, material=clay, name="Half-space")
q = prj.loads.line_load((2.4, 4.0), (3.6, 4.0), qy=-60.0, name="Footing pressure")
prj.initial(procedure="gravity", exclude=[q])
prj.phases.plastic("Load to collapse", activate=[q])
expect_bytes(prj, "kv-fnd-010-prandtl-strip-footing.k2d")

# ---------------------------------------------------------------- KV-FND-011 --
prj = katai.Project("KV-FND-011 Gibson strip load", element="tri15",
                    mesh_size=0.15, auto_refine=False)
gibson = prj.materials.linear_elastic("Gibson soil", E=0.01, nu=0.495,
                                      gamma=0.0, E_inc=299.0, y_ref=4.0)
prj.geometry.rectangle(0.0, 0.0, 7.0, 4.0, material=gibson, name="Layer")
q = prj.loads.line_load((0.0, 4.0), (1.0, 4.0), qy=-10.0, name="Strip q")
prj.initial(procedure="gravity", exclude=[q])
prj.phases.plastic("Strip load", activate=[q])
expect_bytes(prj, "kv-fnd-011-gibson-strip-load.k2d")

# ---------------------------------------------------------------- KV-FND-012 --
GR_G, GR_NU = 500.0, 1.0 / 3.0
prj = katai.Project("KV-FND-012 Giroud rigid footing", element="tri15",
                    mesh_size=0.5, auto_refine=False)
soil = prj.materials.linear_elastic("Elastic soil", E=2.0 * GR_G * (1.0 + GR_NU),
                                    nu=GR_NU, gamma=0.0)
prj.geometry.rectangle(0.0, 0.0, 7.0, 4.0, material=soil, name="Half-space")
footing = prj.displacements.line((0.0, 4.0), (1.0, 4.0), uy=-0.010, name="Footing")
prj.initial(procedure="gravity", exclude=[footing])
prj.phases.plastic("Indent", activate=[footing])
expect_bytes(prj, "kv-fnd-012-giroud-rigid-footing.k2d")

# ------------------------------------- the prescribed flux, from the easy surface --
# KV-FLW-002 verifies the physics in C++; what is checked here is that an engineer
# can ASK for it in one line, and that a model with no flux edge still writes the
# empty array the corpus files above are pinned to.
prj = katai.Project("KV-FLW-002 prescribed flux", mesh_size=1.0, auto_refine=False)
sand = prj.materials.linear_elastic("Sand", E=2.0e4, nu=0.3, gamma=18.0,
                                    gamma_sat=20.0, k=2.0)
prj.geometry.rectangle(0.0, 0.0, 4.0, 10.0, material=sand, name="Column",
                       flow={"bottom": ("head", 12.0), "top": ("flux", 0.05)})
built = prj.build()
text = katai.project_to_json(built)
check('"edge_flow":[1,0,3,0]' in text, "the DSL writes the flux edge kind (3)")
check(list(built.polygons[0].edge_flux) == [0.0, 0.0, 0.05, 0.0],
      "and the rate beside it, on that edge only")
check(katai.validate_project(built).ok, "the flux model validates")
check('"edge_flux":[]' in katai.project_to_json(build_slope().build()),
      "a model with no flux edge still writes the empty array the corpus pins")
try:
    prj.geometry.rectangle(0.0, 0.0, 1.0, 1.0, material=sand, flow={"top": "flux"})
    check(False, '"flux" without a value is refused')
except ValueError as exc:
    check("m/day" in str(exc), f'"flux" without a value is refused: {exc}')

# ------------------------------------ the staged options, from the easy surface --
# KV-EXC-002 and KV-CST-005 verify what these DO; this is the one line an engineer
# has to write to ask for them.
prj = katai.Project("staged options", mesh_size=2.0, auto_refine=False)
soil = prj.materials.mohr_coulomb("Clay", E=5000.0, nu=0.3, c=20.0, phi=25.0,
                                  gamma=17.0, gamma_sat=19.0, drainage="undrained_a")
prj.geometry.rectangle(0.0, 0.0, 10.0, 5.0, material=soil, name="Block")
prj.initial()
prj.phases.plastic("Half a lift", apply_fraction=0.5, ignore_undrained=True)
built = prj.build()
check(built.phases[0].sum_mstage == 0.5, "the DSL sets the staged-construction fraction")
check(built.phases[0].ignore_undrained, "and the ignore-undrained switch")
check('"mstage":0.5' in katai.project_to_json(built)
      and '"ignoreund":true' in katai.project_to_json(built),
      "both reach the file, where a reviewer can see them")

# ------------------------------- the pore fluid's stiffness, from the easy surface --
# KV-CST-006 verifies what these DO; this is how an engineer with a measured Skempton B
# says so instead of accepting PLAXIS's 0.495 for a soil that was never asked about.
prj = katai.Project("undrained stiffness", mesh_size=2.0, auto_refine=False)
soft = prj.materials.mohr_coulomb("Soft clay", E=5000.0, nu=0.3, c=20.0, phi=25.0,
                                  gamma=17.0, drainage="undrained_a",
                                  und_mode=1, skempton_B=0.9)
stiff = prj.materials.mohr_coulomb("Stiff clay", E=5.0e4, nu=0.25, c=60.0, phi=25.0,
                                   gamma=19.0, drainage="undrained_a", nu_u=0.498)
prj.geometry.rectangle(0.0, 0.0, 10.0, 3.0, material=soft, name="Soft")
prj.geometry.rectangle(0.0, 3.0, 10.0, 6.0, material=stiff, name="Stiff")
prj.initial()
built = prj.build()
check(built.materials[0].und_mode == 1 and built.materials[0].skempton_B == 0.9,
      "the DSL sets Skempton's B on one material")
check(built.materials[1].nu_u == 0.498,
      "and the equivalent undrained Poisson ratio on the other")
check('"und_mode":1' in katai.project_to_json(built)
      and '"skempton_B":0.9' in katai.project_to_json(built)
      and '"nu_u":0.498' in katai.project_to_json(built),
      "all three reach the file, per material")
check(katai.validate_project(built).ok, "two differently-watered clays validate")

# ---------------------------------- Undrained (C): the total-stress drainage type --
# KV-CST-007 verifies what it DOES; this is the one word that asks for it.
prj = katai.Project("undrained C", mesh_size=2.0, auto_refine=False)
total = prj.materials.mohr_coulomb("Clay (total stress)", E=1.5e4, nu=0.495, c=60.0,
                                   phi=0.0, gamma=18.0, gamma_sat=20.0,
                                   drainage="undrained_c")
prj.geometry.rectangle(0.0, 0.0, 10.0, 5.0, material=total, name="Clay")
prj.initial()
built = prj.build()
check(built.materials[0].drainage == katai._core.Drainage.UndrainedC,
      "the DSL asks for a total-stress analysis by name")
check('"drainage":4' in katai.project_to_json(built), "and it reaches the file as drainage 4")
check(katai.validate_project(built).ok, "an undrained total-stress clay validates")

# --------------------------------- wells and drains, from the easy surface --
# KV-FLW-003 verifies what they DO; this is the sentence an engineer writes to dewater
# a pit: pump this much out of that line, and hold the head down along this one.
prj = katai.Project("dewatering", mesh_size=2.0, auto_refine=False)
sand = prj.materials.mohr_coulomb("Sand", E=3.0e4, nu=0.3, c=1.0, phi=32.0,
                                  gamma=18.0, k=2.0)
prj.geometry.rectangle(0.0, 0.0, 40.0, 8.0, material=sand, name="Aquifer",
                       flow={"left": ("head", 10.0), "right": ("head", 10.0)})
pump = prj.dewatering.well((20.0, 0.0), (20.0, 8.0), q=1.6, h_min=2.0, name="Pump")
trench = prj.dewatering.drain((5.0, 6.0), (10.0, 6.0), head=6.0, name="Trench")
prj.water.table(10.0)
prj.initial(exclude=[pump])
prj.phases.plastic("Dewater", activate=[pump])
built = prj.build()
check(len(built.hydros) == 2, "the DSL adds a well and a drain")
check(built.hydros[0].q == 1.6 and built.hydros[0].h_min == 2.0,
      "with the well's discharge and its floor")
check(built.hydros[1].head == 6.0, "and the drain's head")
check(list(built.initial.hydro_active) == [0, 1] and
      list(built.phases[0].hydro_active) == [1, 1],
      "and switches the well on in the phase that pumps, like any other object")
check('"hydros":[' in katai.project_to_json(built) and '"hydro":[0,1]' in katai.project_to_json(built),
      "both reach the file, with their per-phase activity")
check(katai.validate_project(built).ok, "the dewatered model validates")

# ------------------------------------------------- end to end: the slope RUNS --
job_dsl = build_slope().run()
file_prj, _ = katai.load_project(f"{CORPUS}/kv-slp-001-griffiths-lane-slope.k2d")
job_file = katai.run(file_prj)
fos_dsl = job_dsl.results()[-1].fos
fos_file = job_file.results()[-1].fos
check(fos_dsl == fos_file,
      f"DSL run == file run, bit for bit (FoS = {fos_dsl:.3f})")
check(abs(fos_dsl - 0.99) < 0.08 * 0.99,
      "FoS within the declared band of the published benchmark")

if failures:
    sys.exit(f"{failures} check(s) failed")
print("\nOK: the easy surface writes exactly the contract the corpus pins")
