#!/usr/bin/env python3
"""Input-safety gate: nothing a user drew may be skipped without saying so.

The driver turns a .k2d project into engine objects one build loop at a time --
loads, plates, geogrids, anchors, walls, interfaces, piles, prescribed
displacements. Each loop used to skip whatever it could not attach to the mesh
with a bare ``continue``, and a skipped object produces no error: the phase runs,
converges, and answers a model the engineer never drew. A load drawn five
centimetres above the ground surface simply vanished; on a weightless benchmark
that shows as a zero displacement, and in a real model the self-weight settlement
hides it completely.

That class of defect is closed one call at a time in the driver. This gate keeps
it closed: inside a loop over user-authored objects, every ``continue`` must be
one of three things, and the third one costs a sentence.

  1. SELECTION -- the guard filters on identity or activity ("not a plate",
     "not installed in this phase"). Skipping an object the phase never asked
     for is not a drop, so these are recognised by their predicates
     (SELECTION_PATTERNS below).
  2. DIAGNOSED -- a warn() / refuse() / diag() call within LOOKBACK lines above
     it. The user is told, with a stable code.
  3. ANNOTATED -- the line, or a comment just above it, carries
     ``silent-drop-ok: <reason>``. Use this only when nothing is lost, and say
     WHY: usually because the input contract already refuses the case at its
     field path, so the object cannot reach here in a run that got this far.

Anything else fails, with the file and line. The point is not the count; it is
that a new bare ``continue`` in a build loop cannot be added without a decision.

A whole loop can also be declared out of scope, on or just above its header:

    // silent-drop-scope: none -- this scans region pairs to decide whether the
    // K0 field is level; it builds nothing, so skipping a pair discards nothing.

Prefer that to annotating four continues that all say the same thing, and put
it where the reader meets the loop.

Usage:
    python scripts/check_silent_drop.py           # gate
    python scripts/check_silent_drop.py --list    # print every classified continue
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# The composition points where a user's drawing becomes an engine object. Adding a
# second driver to this list is how the gate follows the code.
SOURCES = ("kernel/jobs/src/driver.cpp",)

# A loop is in scope when it iterates over user-authored objects. `walls` and
# `soil_ifaces` are the driver's own staging vectors, but each entry is one drawn
# structure, so dropping one drops a wall the engineer put in the file.
USER_COLLECTIONS = (
    "pr.structs", "pr.loads", "pr.disps", "pr.polygons", "pr.embedded",
    "pr.materials", "pr.phases", "walls", "soil_ifaces",
)

# Guards that select rather than discard: the object is not of this kind, or not
# installed in this phase. These describe an object the phase never asked for.
SELECTION_PATTERNS = (
    r"\.kind\s*!=",
    r"!\s*\w*_on\s*\(",           # !struct_on(si), !load_on(li), !poly_on(pi), !disp_on(di)
    r"!\s*\w*active\w*\s*\(",
    r"plate_is_wall\[",
    r"\bskip\w*\b",
)

DIAGNOSTIC_CALL = re.compile(r"\b(warn|refuse|diag)\s*\(\s*R\s*,")
ANNOTATION = re.compile(r"silent-drop-ok:\s*(\S.*\S)")
SCOPE_OUT = re.compile(r"silent-drop-scope:\s*none\s*(?:--|—)\s*(\S.*\S)")
SCOPE_LOOKBACK = 4     # lines above a loop header that may carry its scope note
LOOP_HEAD = re.compile(r"^\s*(for|while)\s*\(")
CONTINUE = re.compile(r"\bcontinue\s*;")
LOOKBACK = 25          # lines above a continue that may carry its diagnostic
ANNOTATION_LOOKBACK = 6


def strip_noise(line: str) -> str:
    """Braces inside a string literal or a comment must not move the depth."""
    line = re.sub(r'"(\\.|[^"\\])*"', '""', line)
    line = re.sub(r"'(\\.|[^'\\])*'", "''", line)
    return re.sub(r"//.*$", "", line)


def loop_stack_at(lines: list[str]) -> list[list[int]]:
    """For each line, the stack of enclosing loop-header line numbers (0-based).

    Brace counting, not parsing: it is enough because the driver is formatted in
    the house style (a loop body's braces are balanced on their own lines) and
    because a miscount would show up immediately as a misclassified continue.
    """
    stacks: list[list[int]] = []
    stack: list[tuple[int, int]] = []      # (header line, depth OUTSIDE the body)
    depth = 0
    pending: int | None = None             # a loop header whose '{' has not been seen
    for i, raw in enumerate(lines):
        text = strip_noise(raw)
        stacks.append([h for h, _ in stack])
        if pending is not None and "{" in text:
            stack.append((pending, depth))
            pending = None
        elif LOOP_HEAD.match(text) and "{" in text:
            stack.append((i, depth))
        elif LOOP_HEAD.match(text):
            # A single-statement loop body (no brace) cannot contain a continue
            # that this gate cares about, so only brace bodies are tracked.
            pending = i if not text.rstrip().endswith(";") else None
        depth += text.count("{") - text.count("}")
        while stack and depth <= stack[-1][1]:
            stack.pop()
    return stacks


def classify(lines: list[str], idx: int, loop_head: int) -> tuple[str, str]:
    """(verdict, detail) for the continue on line idx inside loop_head."""
    line = lines[idx]
    if any(re.search(p, line) for p in SELECTION_PATTERNS):
        return "selection", line.strip()
    for j in range(idx, max(idx - ANNOTATION_LOOKBACK, loop_head) - 1, -1):
        hit = ANNOTATION.search(lines[j])
        if hit:
            return "annotated", hit.group(1)
    for j in range(idx - 1, max(idx - LOOKBACK, loop_head) - 1, -1):
        if DIAGNOSTIC_CALL.search(lines[j]):
            return "diagnosed", lines[j].strip()
    return "SILENT", line.strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list", action="store_true", help="print every classified continue")
    args = ap.parse_args()

    failures: list[str] = []
    counts = {"selection": 0, "diagnosed": 0, "annotated": 0, "scoped-out": 0, "SILENT": 0}
    for rel in SOURCES:
        path = REPO_ROOT / rel
        if not path.is_file():
            print(f"{rel}: not found -- the gate names a file that moved")
            return 2
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        stacks = loop_stack_at(lines)
        for i, raw in enumerate(lines):
            if not CONTINUE.search(strip_noise(raw)):
                continue
            stack = stacks[i]
            if not stack:
                continue
            head = stack[-1]                       # continue binds to the INNERMOST loop
            if not any(c in lines[head] for c in USER_COLLECTIONS):
                continue
            if any(SCOPE_OUT.search(lines[j])
                   for j in range(max(head - SCOPE_LOOKBACK, 0), head + 1)):
                counts["scoped-out"] += 1
                if args.list:
                    print(f"{rel}:{i + 1}: [scoped-out] {raw.strip()}")
                continue
            verdict, detail = classify(lines, i, head)
            counts[verdict] += 1
            if args.list:
                print(f"{rel}:{i + 1}: [{verdict}] {detail}")
            if verdict == "SILENT":
                failures.append(
                    f"{rel}:{i + 1}: a user-authored object is skipped in silence:\n"
                    f"    {raw.strip()}\n"
                    f"    (loop at line {head + 1}: {lines[head].strip()})\n"
                    "    Emit a diagnostic -- refuse(R, ...) when the object would carry "
                    "nothing, warn(R, ...) when it is used differently than drawn -- or, if "
                    "nothing is lost, say why with a 'silent-drop-ok: <reason>' comment.")

    if failures:
        print("check_silent_drop: FAILED\n")
        for f in failures:
            print(f + "\n")
        return 1
    print(f"check_silent_drop: OK -- {counts['selection']} selection, "
          f"{counts['diagnosed']} diagnosed, {counts['annotated']} annotated, "
          f"{counts['scoped-out']} in loops declared out of scope; "
          "0 silent drops in the driver's build loops.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
