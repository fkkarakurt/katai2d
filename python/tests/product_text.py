"""No other program is named anywhere in this project, and no other vendor's manual is cited.

The rule is the maintainer's. It covers everything the project contains and publishes: source
comments as well as messages, the Python package and its bindings, the command line, the
validation record, the changelog, the README, the Studio's labels, reports and installer. An
engineering statement stands on its own -- a formula written out, a definition, a closed form,
an academic primary source -- and needs no other program to lend it weight.

HOW THE GATE READS. Every tracked text file under each root given on the command line is split
into word tokens: letters and digits apart, camelCase and snake_case apart, lower-cased. Each
token, and each run of two or three tokens separated only by white space, a comment marker, a
slash or a hyphen, is compared with a set of SHA-256 digests. The names are kept as digests so
that this file does not itself become a list of other programs. To add one, lower-case it, drop
everything that is not a letter or a digit, and add the digest of what is left; for a spelling
that writes a number apart from a word, add the digest with an underscore at that gap as well.

Only TRACKED files are read (``git ls-files`` in each root), so a private working note that git
does not carry is not in scope. Without git the roots are walked instead, skipping third-party
code, build and distribution output, and the private note directories.

The point of a gate rather than a one-time edit: a sentence written next year is exactly where
this would come back, and nobody greps for a rule they were not there for.
"""

import hashlib
import re
import subprocess
import sys
from pathlib import Path

DIGESTS = frozenset({
    "4dc291ad2a05eb712070be663cfc47b0f2ad94cb19f665c647e5fcb7c715877c",
    "3fd725a89717d2175289be8011e2b983ab6428bc0064f788f1a9d6eb014c2253",
    "fee6fd0cc1707a180d7981bec65872f618dc323a577e3ac95a24a3b24e8bf52d",
    "559054694ada46f586fa9ad4cc7874eed1077162eae10e0d66d70671915d219f",
    "0e364c206974b5d9564601b5a8c51a27131ccc514af2af4b8904cc0b4e3bf06b",
    "920ca0bde4744c77125c3118a1107af43930979aa2cbc70d1dd6b1ca442c7703",
    "c8cff8e3d42a449a60e4bf1e813b3253a8c9fd428b9384e03a7bdbce1d389d3f",
    "c9fe0a3d6a791aab584bacca52ca327a835c9e10ccd8198448b86522e9da5d89",
    "88f2c31fd6c910978f00f7ad9696cd24907a713b6a19c070ae95c15466faf7ea",
    "232b00645540ae6d21ce88a90e6791f01495adaf959c7d645419470d0d938267",
    "6be8f2da65ac4631618b022a9fab609377eaad5fa14e850cde24f63fe9eb0029",
    "9cc78f0ad9614b073638536561207e9257289f643e72d4e8173a1466fa540f5f",
    "01de1183502ff18d17906c32b2df6bde9b4ce67f92eb0be2ffd6813da3a92c2c",
    "0d9d62c60663ff617373f39814d2ec83ea8c50fe4aa3075513c3dab272c3a274",
    "3a9023d5efaf19c2ef7d8f82da962a6a1a91c89911298c94fbad7d8525db1b8f",
    "8b7293f6dd93853c83bf8c87ec5349a41789fa166a9431f108edb538e6fd82fd",
    "6bd496a3198b9685e5cf867723cad4e5fedbd7e91d68c7bc0c9ecf5da90ffef7",
    "a7d324a4b53ea3878089c268660a04874f40e805f6f74b28bb448073625ea6f0",
    "d7e7f5b711747635399d8c443e56c3df614f5d7af0ed005b2a2d7a6d231e89c9",
    "b48ab14c941506c9e1ba8a0bdedf7292ef390fc9d8a2c3744621ec33afda0289",
    "d40ad4b65864a753f19cab171e92ff0932cdee3250e021b80984360a3c777986",
    "ef82db7ccafc38195df0d0d5f1d270a8cbf29dfeb6122fe6ed4742392c3fde3f",
    "b2b452f81f336c6d4e552dbb139c4d28a62a32c8a058a59e1f9b7bf08a64f871",
    "b604257892f6c012bf913cef1edd2d3a489850457173b549795f4551eb6b45e9",
    "70bf119c8264eed16123dc0cb09f70cae1d6a6a23360a694013ae56e1e8885ea",
    "25c9f1d0e767f78a909e05e841e79bff4b14d4006d0fe6d9e62b4d433f973fff",
    "052224a445e410b071ebdf2668bd24bf21b7912654908dcedd1260eb80046889",
    "61c869b2a67e8b41363dee62342973a022a1dbb24a1a75b54e9a994e3c89b232",
    "6df2ddcd93f7719baec31d3e230f4767b7224ac9fcddb727dbb8fabac9c05017",
    "659353384faa06ebba3f13b09e22c985871f81c77b36f11702c63973178e3525",
    "f75d60b570154c433f1455db73e4f611ba376f5f2e98cd2c4dd3af14e17b204c",
    "af5fa6ec9bf26b46ece0620c0cbc3a16f602353ba91e5a3a385e58c0d86b40dc",
    "b9fb9042c23c00fec65e70305b2d1271d2baabb48b98593d0e0d8f38402c9cde",
    "3073daf535a6fa12954d2d189538a1a71c2b1944c5ac60e15620684a3360c35c",
    "4fd6275b5f39fcd986b980aa470ebdd1ff8df7310183e95a2953d330e69d2a2e",
    "d678cc3f5994039e15fb1aee17e5eebda4850da7e5aa21141b52161c47e0053b",
})
MAX_LEN = 20   # no name is longer, so longer candidates are not hashed

TOKEN = re.compile(r"[A-Z]?[a-z]+|[A-Z]+(?![a-z])|[0-9]+")
# What may stand between two tokens of one name: white space, comment markers, a slash, a hyphen.
JOIN = re.compile(r"[\s/#*\-]*(?://)?[\s/#*\-]*")

SKIP_PARTS = {"third_party", ".git", "__pycache__", "node_modules"}
SKIP_TOP = ("build", "dist", "studio")
PRIVATE = ("docs/internal/", "docs/references/")


def _digest(s):
    return hashlib.sha256(s.encode("ascii")).hexdigest()


def tracked_files(root):
    try:
        out = subprocess.run(["git", "-C", str(root), "ls-files", "-z"], capture_output=True,
                             check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    return [root / p for p in out.decode("utf-8", "replace").split("\0") if p]


def walked_files(root):
    files = []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        rel = p.relative_to(root).as_posix()
        if SKIP_PARTS.intersection(p.parts) or rel.startswith(PRIVATE):
            continue
        if rel.split("/", 1)[0].startswith(SKIP_TOP):
            continue
        files.append(p)
    return files


def scan_text(text):
    """Return (offset, line) for every name found in `text`."""
    toks = [(m.start(), m.end(), m.group(0).lower()) for m in TOKEN.finditer(text)]
    hits = []
    for i, (start, end, tok) in enumerate(toks):
        cand, last_end = tok, end
        for n in range(3):
            if n > 0:
                j = i + n
                sep = text[last_end:toks[j][0]] if j < len(toks) else None
                if sep is None or not JOIN.fullmatch(sep):
                    break
                # A number written apart from a word is marked, so "phase 2" and a name that is
                # a word run into a digit hash differently.
                mark = "_" if sep and (toks[j][2].isdigit() or cand[-1].isdigit()) else ""
                cand += mark + toks[j][2]
                last_end = toks[j][1]
            if len(cand) > MAX_LEN:
                break
            if _digest(cand) in DIGESTS:
                hits.append(start)
                break
    return hits


def scan(root):
    root = Path(root).resolve()
    files = tracked_files(root)
    if files is None:
        files = walked_files(root)
    findings = []
    for path in sorted(files):
        rel = path.relative_to(root).as_posix()
        if SKIP_PARTS.intersection(Path(rel).parts) or rel.startswith(PRIVATE):
            continue
        try:
            data = path.read_bytes()
        except OSError:
            continue
        if b"\0" in data[:8192]:
            continue   # binary
        text = data.decode("utf-8", "replace")
        seen = set()
        for off in scan_text(text):
            line_no = text.count("\n", 0, off) + 1
            if line_no in seen:
                continue
            seen.add(line_no)
            line = text.splitlines()[line_no - 1].strip()
            findings.append((path, line_no, line[:140]))
    return findings


def main(argv):
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(errors="replace")   # a console code page must not end the report
    roots = argv[1:]
    if not roots:
        print("usage: product_text.py <root> [...]")
        return 2
    findings = []
    for r in roots:
        if Path(r).exists():
            findings += scan(r)
        else:
            print(f"  (skipped, absent: {r})")
    for path, line_no, line in findings:
        print(f"  {path}:{line_no}: {line}")
    if findings:
        print(f"\nproduct_text: {len(findings)} line(s) name another program or cite its manual")
        return 1
    print("ok: no other program is named and no other vendor's manual is cited")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
