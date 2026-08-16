"""Readable output for a scripted run.

The raw result carries everything -- displacements, stresses, reactions, structural
force envelopes, the engine's diagnostics -- but as arrays, and a script that wants
to SEE what happened otherwise has to hand-roll a numpy scan and guess the shape of
each object. That is the difference between data being available and output being
usable, and it is the whole of this module.

Nothing here computes anything the engine did not already produce: it selects,
labels, and formats. Two functions do the work::

    print(katai.summary(job, prj))     # every phase, named, with its extremes
    ux = katai.extremes(res)["u_x"]    # (min, at, max, at) for one field

This is the SCRIPTING output -- plain text, no layout, meant for a terminal, a log
or a notebook. The formatted client deliverable (corporate header, typeset notation,
figures, print pagination) is the Studio report and stays there.

Units are the schema's own and fixed: kN, m, day.
"""

import math

from . import _core

# Nodal fields, in the order a reader wants them: what moved, then what it cost.
# (name, unit, extractor). Kept here rather than imported so this module has no
# dependency on the GUI's field enum.
_FIELDS = (
    ("|u|",     "m",   lambda d, s, p, n: math.hypot(d[2 * n], d[2 * n + 1])),
    ("u_x",     "m",   lambda d, s, p, n: d[2 * n]),
    ("u_y",     "m",   lambda d, s, p, n: d[2 * n + 1]),
    ("sig'_xx", "kPa", lambda d, s, p, n: s[n][0]),
    ("sig'_yy", "kPa", lambda d, s, p, n: s[n][1]),
    ("sig_xy",  "kPa", lambda d, s, p, n: s[n][2]),
    ("tau_max", "kPa", lambda d, s, p, n: math.hypot(0.5 * (s[n][0] - s[n][1]), s[n][2])),
    ("p_w",     "kPa", lambda d, s, p, n: p[n] if n < len(p) else 0.0),
)

_KIND = {0: "plate", 1: "anchor", 2: "geogrid", 3: "pile", 4: "interface"}
_SEVERITY = {0: "NOTE", 1: "WARN", 2: "STOP"}


def extremes(result):
    """``{field: (vmin, (x, y), vmax, (x, y))}`` over every node of one phase.

    The location is part of the answer: "the largest settlement is 43 mm" and
    "...and it is under the wall toe, not the excavation floor" are two different
    findings, and a bare maximum only gives the first. A field that is the same
    everywhere reports ``None`` for both places rather than inventing one.
    """
    d = result.displacement
    s = result.stress
    p = list(result.pore)
    x, y = list(result.node_x), list(result.node_y)
    n_nodes = len(x)
    out = {}
    for name, unit, get in _FIELDS:
        if n_nodes == 0:
            continue
        lo = hi = get(d, s, p, 0)
        nlo = nhi = 0
        for n in range(1, n_nodes):
            v = get(d, s, p, n)
            if v < lo: lo, nlo = v, n
            if v > hi: hi, nhi = v, n
        flat = lo == hi
        out[name] = (lo, None if flat else (x[nlo], y[nlo]),
                     hi, None if flat else (x[nhi], y[nhi]))
    return out


def _fmt_point(pt):
    return "  -  " if pt is None else f"({pt[0]:.3g}, {pt[1]:.3g})"


def summary(what, project=None, *, fields=True, structures=True, diagnostics=True):
    """A readable account of a run: ``what`` is a Job or a single result.

    Pass the ``project`` (the :class:`katai.Project` or the built schema project)
    and each phase is named the way the script named it -- a result does not carry
    its own phase name, so without it the phases can only be numbered.

    Set ``fields=False`` / ``structures=False`` / ``diagnostics=False`` to leave a
    block out; by default nothing that could change how a number is read is hidden.
    """
    results = what.results() if hasattr(what, "results") else [what]
    names = _phase_names(project, len(results))
    return "\n".join(
        _one(r, names[i], fields, structures, diagnostics)
        for i, r in enumerate(results))


def _phase_names(project, n):
    """Phase names in result order: the initial phase, then the staged ones."""
    if project is None:
        return [f"Phase {i}" for i in range(n)]
    pr = project.build() if hasattr(project, "build") else project
    names = ["Initial phase"] + [ph.name for ph in pr.phases]
    names += [f"Phase {i}" for i in range(len(names), n)]
    return names[:n]


def _one(r, name, fields, structures, diagnostics):
    L = [f"=== {name} " + "=" * max(0, 60 - len(name))]
    status = "converged" if r.ok else "NOT CONVERGED"
    L.append(f"  status        {status}")
    if not r.ok and r.message:
        L.append(f"  reason        {r.message}")
    L.append(f"  max |u|       {r.max_disp:.6e} m")
    if r.load_factor != 1.0:
        L.append(f"  load factor   {r.load_factor:.4f}"
                 "   <- the fraction that reached equilibrium")
    if r.fos >= 0.0:
        L.append(f"  factor of safety {'>' if r.fos_is_lower_bound else '='} {r.fos:.4f}")

    # The engine's qualifications first: they change how everything below is read.
    if diagnostics and r.diagnostics:
        L.append("  diagnostics")
        for d in r.diagnostics:
            tag = _SEVERITY.get(int(d.severity), "NOTE")
            subj = f"  [{d.subject}]" if d.subject else ""
            L.append(f"    {tag}  {d.code}{subj}")
            L.append(f"          {d.message}")

    if fields:
        ex = extremes(r)
        if ex:
            L.append("  extremes                    min          at"
                     "                 max          at")
            for fname, unit, _ in _FIELDS:
                if fname not in ex:
                    continue
                lo, plo, hi, phi = ex[fname]
                L.append(f"    {fname:<9s} [{unit:>3s}]  {lo:>12.5g}  {_fmt_point(plo):<18s}"
                         f" {hi:>12.5g}  {_fmt_point(phi)}")

    if structures and r.struct_forces:
        L.append("  structural forces (envelope |max|)")
        for sf in r.struct_forces:
            kind = _KIND.get(int(sf.kind), "?")
            flag = "   YIELDED" if sf.yielded else ""
            L.append(f"    {sf.name:<22s} {kind:<9s} N {sf.max_N:>11.4g} kN/m"
                     f"   Q {sf.max_Q:>11.4g} kN/m   M {sf.max_M:>11.4g} kNm/m{flag}")

    if r.consol_time:
        L.append(f"  consolidation  {len(r.consol_time)} steps, "
                 f"t = {r.consol_time[-1]:.4g} day, settlement "
                 f"{r.consol_settlement[-1]:.6e} m, excess pore "
                 f"{r.consol_excess_pore[-1]:.4g} kPa")
    return "\n".join(L)
