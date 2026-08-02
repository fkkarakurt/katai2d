"""Python smoke: the corpus runs through the bindings, and refusals stay honest.

Not the oracle record -- the physics numbers are declared and asserted by
test_input_corpus (C++) and the corpus files themselves. This test pins the
PLUMBING: every corpus case loads and validates through katai._core, a
single-phase and a staged case RUN through katai.run() with progress, results
come back as arrays of the right shape, and an invalid project raises Refusal
carrying the validator's field path.
"""
import sys

import katai

CORPUS = sys.argv[1]
CASES = [
    "kv-con-002-terzaghi-column.k2d",
    "kv-fnd-008-strip-load.k2d",
    "kv-fnd-009-flamant-line-load.k2d",
    "kv-num-003-k0-geostatic-block.k2d",
    "kv-cst-001-undrained-column.k2d",
    "kv-slp-001-griffiths-lane-slope.k2d",
    "kv-exc-001-staged-excavation.k2d",
]

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


# (1) Every corpus case loads clean and validates clean through the bindings.
for name in CASES:
    project, notes = katai.load_project(f"{CORPUS}/{name}")
    report = katai.validate_project(project)
    check(not notes and report.ok(), f"{name} loads + validates clean")

# (2) A single-phase case runs; results have honest shapes.
project, _ = katai.load_project(f"{CORPUS}/kv-fnd-008-strip-load.k2d")
announced = []
job = katai.run(project, on_phase=lambda i, n, name: announced.append((i, n, name)))
check(job.state() == katai.JobState.Done, "strip-load job Done")
res = job.results()[-1]
check(announced == [(0, 1, "Initial phase")], "progress announced the single phase")
n = len(res.node_x)
check(n > 0 and len(res.displacement) == 2 * n, "displacement is a 2n vector")
check(res.stress.shape == (n, 3), "stress is (n, 3)")
check(res.max_disp > 0.0, "the load produced displacement")
check(katai.project_to_json(project) == katai._core.project_to_json(project),
      "package and _core agree on the canonical JSON")

# (3) A staged case runs: two phases, consolidation series present.
project, _ = katai.load_project(f"{CORPUS}/kv-con-002-terzaghi-column.k2d")
job = katai.run(project)
check(len(job.results()) == 2, "Terzaghi ran initial + consolidation")
check(len(job.results()[-1].consol_time) > 0, "consolidation time series present")

# (4) Refusal: an invalid project raises with the field path, honestly.
project, _ = katai.load_project(f"{CORPUS}/kv-fnd-008-strip-load.k2d")
project.materials[0].E = -1.0
try:
    katai.run(project)
    check(False, "invalid project raised Refusal")
except katai.Refusal as e:
    check("materials[0].E" in str(e), "Refusal names the field path")
    check(e.report is not None and not e.report.ok(), "Refusal carries the report")

# (5) Provenance is visible.
check(isinstance(katai.backend_name(), str) and katai.backend_name() != "",
      f"backend: {katai.backend_name()}")

if failures:
    sys.exit(f"{failures} check(s) failed")
print("\nOK: the corpus runs through Python, and the bindings cannot change the answer")
