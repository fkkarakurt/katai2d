"""No other FE program is named in anything a user of this program reads.

The rule is the maintainer's and it is about the PRODUCT surface: the Studio's
labels, tooltips and reports, the command line's output, and the Python package's
docstrings -- everything `help()`, a tooltip or a printed page can show. It is not
about the validation record, where naming the program a number was compared
against is the whole point of the comparison, nor about a source comment that
records where a default came from.

So this gate reads STRING LITERALS AND DOCSTRINGS ONLY, and it needs a real
scanner to do that: a line-based one would be fooled by a `//` inside a string
and by a name sitting after one. What it walks is given on the command line --
the open tree's `cli/` and `python/katai/`, and the Studio's `app/` when the
Studio tree is present.

The point of a gate rather than a one-time edit: a sentence written next year
under a tooltip is exactly where this would come back, and nobody greps for a
rule they were not there for.
"""

import re
import sys
from pathlib import Path

NAMES = re.compile(
    r"\b(PLAXIS|Plaxis|OpenSees|FLAC|ADONIS|HYRCAN|GEO ?5|GTS ?NX|Abaqus|ANSYS|"
    r"Midas|MIDAS|Rocscience|GeoStudio|SEEP/W|Settle3|Slide2|RS2|Optum|ZSoil)\b",
    re.I,
)


def cpp_literals(text):
    """Yield (line_no, literal) for every C++ string literal outside a comment."""
    out = []
    line_no = 1
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == "\n":
            line_no += 1
            i += 1
        elif c == "/" and i + 1 < n and text[i + 1] == "/":
            while i < n and text[i] != "\n":
                i += 1
        elif c == "/" and i + 1 < n and text[i + 1] == "*":
            i += 2
            while i + 1 < n and not (text[i] == "*" and text[i + 1] == "/"):
                if text[i] == "\n":
                    line_no += 1
                i += 1
            i += 2
        elif c == '"':
            start, start_line = i, line_no
            i += 1
            while i < n and text[i] != '"':
                if text[i] == "\\":
                    i += 1
                elif text[i] == "\n":
                    line_no += 1
                i += 1
            out.append((start_line, text[start:i + 1]))
            i += 1
        else:
            i += 1
    return out


def py_strings(text):
    """Yield (line_no, string) for every Python string literal, docstrings included."""
    out = []
    pat = re.compile(r'(?s)("""(?:[^"\\]|\\.|"(?!""))*"""|\'\'\'.*?\'\'\'|'
                     r'"(?:[^"\\\n]|\\.)*"|\'(?:[^\'\\\n]|\\.)*\')')
    for m in pat.finditer(text):
        # A '#' comment cannot open a string, but a string can contain '#'. Only skip a
        # match whose line begins with a comment marker before the quote.
        line_start = text.rfind("\n", 0, m.start()) + 1
        before = text[line_start:m.start()]
        if "#" in before and before.split("#", 1)[0].count('"') % 2 == 0 \
                and before.split("#", 1)[0].count("'") % 2 == 0:
            continue
        out.append((text.count("\n", 0, m.start()) + 1, m.group(1)))
    return out


def scan(root):
    findings = []
    root = Path(root)
    files = []
    if root.is_file():
        files = [root]
    else:
        for ext in ("*.cpp", "*.hpp", "*.h", "*.py"):
            files += [p for p in root.rglob(ext) if "third_party" not in p.parts
                      and "build" not in p.parts and "__pycache__" not in p.parts]
    for path in sorted(files):
        text = path.read_text(encoding="utf-8", errors="replace")
        pieces = py_strings(text) if path.suffix == ".py" else cpp_literals(text)
        for line_no, lit in pieces:
            hit = NAMES.search(lit)
            if hit:
                findings.append((path, line_no, hit.group(0), lit.strip()[:120]))
    return findings


def main(argv):
    roots = argv[1:]
    if not roots:
        print("usage: product_text.py <dir-or-file> [...]")
        return 2
    findings = []
    for r in roots:
        if Path(r).exists():
            findings += scan(r)
        else:
            print(f"  (skipped, absent: {r})")
    for path, line_no, name, lit in findings:
        print(f"  {path}:{line_no}: names {name} in user-visible text")
        print(f"      {lit}")
    if findings:
        print(f"\nproduct_text: {len(findings)} competitor name(s) in what a user reads")
        return 1
    print("ok: no other FE program is named in the product's own text")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
