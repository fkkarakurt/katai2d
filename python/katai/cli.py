"""``katai`` from pip -- the same command line as the native executable.

Installed as a console script by the wheel, this module mirrors ``cli/main.cpp``
line for line: the same three commands, the same output text, and the SAME exit
codes (0 ok, 2 usage, 3 file error, 4 input contract refused, 5 solve failed) --
scripts written against one front end must not break on the other. Anything that
holds here is pinned by the wheel build's verification step.
"""
import sys

from . import _core


def _print_issues(issues):
    for i in issues:
        sev = "error" if i.severity == _core.Severity.Error else "warning"
        print(f"{sev}: {i.path}: {i.message}")


def _any_error_note(notes):
    return any(i.severity == _core.Severity.Error for i in notes)


def _print_diagnostics(result):
    """What the run did that the file does not literally say: a line clipped to
    the soil, a fallback taken, an object refused. The input contract is checked
    before a mesh exists, so this is the only report on everything decided
    against the mesh -- printing it is what keeps a discarded load from passing
    as a clean run."""
    for d in result.diagnostics:
        sev = ("refusal" if d.severity == _core.DiagnosticSeverity.Refusal else
               "warning" if d.severity == _core.DiagnosticSeverity.Warning else "note   ")
        subject = f"{d.subject}: " if d.subject else ""
        print(f"  {sev} {d.code}  {subject}{d.message}")


def _usage():
    sys.stderr.write(
        "katai -- KATAI 2D command line\n"
        "usage:\n"
        "  katai solve <file.k2d> [--out <file.res>]\n"
        "  katai validate <file.k2d>\n"
        "  katai info\n"
        "exit codes: 0 ok, 2 usage, 3 file error, 4 input contract refused, "
        "5 solve failed\n")
    return 2


def _load(path):
    """Load + surface the reader's forward-version notes. Returns (project,
    notes, 0) or (None, None, exit_code)."""
    try:
        pr, notes = _core.load_project(path)
    except ValueError as e:
        sys.stderr.write(f"error: {e}\n")
        return None, None, 3
    _print_issues(notes)
    return pr, notes, 0


def _cmd_validate(path):
    pr, notes, rc = _load(path)
    if rc:
        return rc
    rep = _core.validate_project(pr)
    sys.stdout.write(str(rep))
    if not rep.ok() or _any_error_note(notes):
        print(f"REFUSED: the project fails the input contract "
              f"({rep.errors} error(s), {rep.warnings} warning(s))")
        return 4
    print(f"OK: {path} satisfies the input contract ({rep.warnings} warning(s))")
    return 0


def _cmd_solve(path, out):
    pr, notes, rc = _load(path)
    if rc:
        return rc
    if _any_error_note(notes):
        print("REFUSED: the file carries reader ERROR notes (see above); "
              "a solve must not proceed on clamped input")
        return 4

    job = _core.Job(pr)
    job.set_on_phase(lambda cur, total, name:
                     sys.stderr.write(f"phase {cur + 1}/{total}: {name}\n"))
    if not job.run():
        if not job.report().ok():
            sys.stdout.write(str(job.report()))
            print(f"REFUSED: {job.message()}")
            return 4
        # Every phase that ran, including the one that stopped the job: its
        # diagnostics name the object at fault, so they survive the failure path.
        for i, r in enumerate(job.results()):
            if r.diagnostics:
                print(f"phase {i + 1}:")
                _print_diagnostics(r)
        print(f"FAILED: {job.message()}")
        return 5

    res = job.results()
    total_s = sum(t.seconds for t in job.timings())
    for i, r in enumerate(res):
        line = f"phase {i + 1}/{len(res)}: ok  max|u| = {r.max_disp:.6e} m"
        if r.fos >= 0.0:
            line += f"  FoS {'>' if r.fos_is_lower_bound else '='} {r.fos:.3f}"
        print(line)
        _print_diagnostics(r)
    print(f"solved {len(res)} phase(s) in {total_s:.2f} s")

    if out:
        model_hash = _core.fnv1a64(_core.project_to_json(pr))
        try:
            _core.save_results(out, model_hash, res)
        except ValueError as e:
            sys.stderr.write(f"error: {e}\n")
            return 3
        print(f"results written to {out}")
    return 0


def _cmd_info():
    print(f"KATAI 2D {_core.__version__} ({_core.__version_date__})")
    print(f"  .k2d project file version: {_core.PROJECT_FILE_VERSION}")
    print(f"  .res results file version: {_core.RESULTS_FILE_VERSION}")
    print(f"  linear solver backend: {_core.backend_name()}")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if not argv:
        return _usage()
    cmd = argv[0]
    if cmd == "info":
        return _cmd_info() if len(argv) == 1 else _usage()
    if cmd == "validate":
        return _cmd_validate(argv[1]) if len(argv) == 2 else _usage()
    if cmd == "solve":
        out, file = "", ""
        i = 1
        while i < len(argv):
            a = argv[i]
            if a == "--out" and i + 1 < len(argv):
                out = argv[i + 1]
                i += 2
            elif a.startswith("-") or file:
                return _usage()
            else:
                file = a
                i += 1
        return _cmd_solve(file, out) if file else _usage()
    return _usage()


if __name__ == "__main__":
    raise SystemExit(main())
