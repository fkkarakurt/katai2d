#!/usr/bin/env python3
"""Source-language gate: fail the build on non-English prose in source files.

KATAI 2D is written by a Turkish-speaking author but is intended to be read,
audited and cited internationally. Every comment, identifier and diagnostic
message in the source tree must therefore be English. Localized user-facing
text is the one legitimate exemption, and it must be declared (see the escapes
at the end of this docstring), never silent.

Detection is deliberately three-tier so that the gate never cries wolf:

  1. Double-encoded UTF-8 (mojibake) anywhere on a line is a definite hit, and
     it is checked FIRST because it hides the other two: a comment written as
     UTF-8 and re-encoded from a legacy-codepage reading of its own bytes has
     no Turkish letters left in it and no recognisable words either. Two such
     comments sat in ``kernel/`` -- a subtree this gate reports as clean --
     until the check existed. The line is reported RECOVERED, since its stored
     form says nothing to a reader.
  2. Turkish-specific letters (g-breve, dotless i, s-cedilla, ...) anywhere on a
     line are a definite hit. C++ identifiers cannot contain them, so any
     occurrence outside an exempt string literal is non-English prose.
  3. ASCII-folded Turkish function words (``icin``, ``degil``, ``ancak``, ...)
     inside a *comment* are a definite hit. The word list contains only tokens
     that are not English words, so a single occurrence is conclusive. Ambiguous
     tokens that collide with English ("her", "var", "son", "once", "on", "an")
     are intentionally absent.

Three comment syntaxes are recognised: C-style, hash-style, and *prose* files
(``.md``, ``.tex``), which have no comment syntax at all -- the whole line is
prose, so both tiers apply to all of it. Documentation and the manuals ship with
the product and are read by the same audience as the code, so they are inside
the gate. Only genuinely Turkish content is exempted: a citation to a
Turkish-language standard, or the localisation table itself.

Usage:
    python scripts/check_language.py                  # gate: exit 1 on any hit
    python scripts/check_language.py --report         # census, always exit 0
    python scripts/check_language.py --report --by-file
    python scripts/check_language.py kernel/materials # restrict to a subtree

A line may be exempted with a trailing ``katai-lang: allow`` marker, and whole
paths may be exempted in ``scripts/language_policy.json``. Both escapes are for
localized user-facing text only, never for comments.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import unicodedata
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
POLICY_PATH = Path(__file__).resolve().parent / "language_policy.json"
BASELINE_PATH = Path(__file__).resolve().parent / "language_baseline.json"

# Turkish letters that cannot appear in a C++ identifier or in English prose.
TURKISH_LETTERS = re.compile(r"[ğĞıİşŞüÜöÖçÇ]")  # katai-lang: allow (the detector's own alphabet)

# ASCII-folded Turkish function words with no English homograph. One occurrence
# inside a comment is conclusive. Keep this list free of English collisions.
TURKISH_WORDS = frozenset(
    """
    icin degil ancak yani bkz olarak gibi kadar sonra uzere boylece cunku
    ayrica ornegin gerekir gereken gerekli kullanilir kullanan kullanim
    hesaplanir hesaplanan yerine hepsi hicbir sadece yalnizca dolayi
    nedeniyle sayesinde tarafindan olsun olmali olmasi olmadigi yapilir
    yapilan yapmak almak vermek bunun bunlar sunun hangi nedenle icinde
    disinda uzerinde altinda arasinda birlikte ayni farkli buyuk kucuk
    dogru yanlis eksik fazla artik simdilik henuz zaten mutlaka kesinlikle
    dikkat uyari hata cozum deger degeri degerler dosya dosyasi klasor
    satir sutun boyut duzey seviye asama adim adimlar surec islem islemi
    sinir sinirlar kosul kosullar durum durumu tanim tanimli tanimlanan
    ozellik ozellikler yontem yontemi kaynak kaynagi cikti girdi baslangic
    bitis sonuc sonucu ornek katsayi carpan bolen toplam ortalama gerilme
    birim zemin kazik duvar tabaka basinc yuk yukleme mesnet kiris
    """.split()
)

# One extension family per comment syntax.
CPP_SUFFIXES = {".cpp", ".hpp", ".h", ".cc", ".cxx", ".inl", ".ipp"}
HASH_SUFFIXES = {".py", ".cmake", ".ps1", ".sh", ".yml", ".yaml", ".toml"}
HASH_NAMES = {"CMakeLists.txt"}
PROSE_SUFFIXES = {".md", ".tex"}

CPP, HASH, PROSE = "cpp", "hash", "prose"

ALLOW_MARKER = re.compile(r"katai-lang:\s*allow")


def comment_style(path: Path) -> str | None:
    """Which comment syntax applies to this file, or None if it is not scanned."""
    if path.suffix in PROSE_SUFFIXES:
        return PROSE
    if path.suffix in HASH_SUFFIXES or path.name in HASH_NAMES:
        return HASH
    if path.suffix in CPP_SUFFIXES:
        return CPP
    return None


@dataclass(frozen=True)
class Hit:
    path: Path
    line_no: int
    reason: str
    text: str

    def render(self, root: Path) -> str:
        return f"{display_path(self.path, root)}:{self.line_no}: [{self.reason}] {self.text.strip()[:110]}"


def display_path(path: Path, root: Path) -> str:
    """Repo-relative when possible, absolute otherwise (ad-hoc paths, probes)."""
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def demojibake(line: str) -> str | None:
    """Recover a double-encoded line, or None if the line is not one.

    Text written once as UTF-8 and then re-encoded from a legacy-codepage
    reading of those bytes -- s-cedilla stored as the pair A-ring, Y-diaeresis
    -- defeats both detectors above: the two-byte form of every Turkish letter
    becomes a pair of Latin-1 ones, and the folded words never form. The damage
    is recognisable, though -- encoding the line back through the codepage that
    misread it and decoding it as UTF-8 returns the original text, which no
    line that was written correctly can do.

    Found this way: two comments inside kernel/, a subtree the gate reports as
    clean, that had been sitting in its scan set unread. cp1252 is tried first
    because that is the Windows console default that produces these here; the
    two codecs differ only in 0x80-0x9F, and that is exactly where the second
    byte of s-cedilla and g-breve lands.
    """
    for codec in ("cp1252", "latin-1"):
        try:
            recovered = line.encode(codec).decode("utf-8")
        except (UnicodeEncodeError, UnicodeDecodeError):
            continue
        return recovered if recovered != line else None
    return None


def fold(text: str) -> str:
    """ASCII-fold so that 'icin' matches both 'icin' and its accented spelling."""
    swapped = text.replace("ı", "i").replace("İ", "I")  # katai-lang: allow (dotless-i folding)
    decomposed = unicodedata.normalize("NFKD", swapped)
    return "".join(c for c in decomposed if not unicodedata.combining(c)).lower()


def comment_spans(line: str, in_block: bool, style: str) -> tuple[str, bool]:
    """Return the comment text on this line and the block-comment state after it.

    String literals are skipped so that a ``//`` inside a URL literal or a ``#``
    inside a format string is not mistaken for a comment. Prose files have no
    code to separate from commentary, so the whole line is returned.
    """
    if style == PROSE:
        return line, False
    if style == HASH:
        idx = line.find("#")
        return (line[idx + 1 :] if idx >= 0 else ""), False

    out: list[str] = []
    i, n = 0, len(line)
    quote: str | None = None
    while i < n:
        ch = line[i]
        if in_block:
            end = line.find("*/", i)
            if end < 0:
                out.append(line[i:])
                return "".join(out), True
            out.append(line[i:end])
            i, in_block = end + 2, False
            continue
        if quote is not None:
            if ch == "\\":
                i += 2
                continue
            if ch == quote:
                quote = None
            i += 1
            continue
        if ch in "\"'":
            quote = ch
            i += 1
            continue
        if line.startswith("//", i):
            out.append(line[i + 2 :])
            return "".join(out), False
        if line.startswith("/*", i):
            i, in_block = i + 2, True
            continue
        i += 1
    return "".join(out), in_block


def load_policy() -> dict:
    if POLICY_PATH.exists():
        return json.loads(POLICY_PATH.read_text(encoding="utf-8"))
    return {"exempt_paths": [], "exempt_line_patterns": []}


def iter_sources(roots: list[Path], skip: list[re.Pattern]) -> list[Path]:
    found: list[Path] = []
    for root in roots:
        candidates = [root] if root.is_file() else sorted(root.rglob("*"))
        for path in candidates:
            if not path.is_file():
                continue
            if comment_style(path) is None:
                continue
            rel = display_path(path, REPO_ROOT)
            if any(p.search(rel) for p in skip):
                continue
            found.append(path)
    return found


def scan(path: Path, exempt_line: list[re.Pattern]) -> list[Hit]:
    style = comment_style(path)
    if style is None:
        return []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (UnicodeDecodeError, OSError):
        return []

    hits: list[Hit] = []
    in_block = False
    for line_no, line in enumerate(lines, start=1):
        comment, in_block = comment_spans(line, in_block, style)
        if ALLOW_MARKER.search(line) or any(p.search(line) for p in exempt_line):
            continue
        recovered = demojibake(line)
        if recovered is not None:
            hits.append(Hit(path, line_no, "MOJIBAKE", recovered))
            continue
        if TURKISH_LETTERS.search(line):
            hits.append(Hit(path, line_no, "TR-CHAR", line))
            continue
        if comment:
            words = set(re.findall(r"[a-z]+", fold(comment)))
            found = words & TURKISH_WORDS
            if found:
                hits.append(Hit(path, line_no, f"TR-WORD:{sorted(found)[0]}", line))
    return hits


def census(hits: list[Hit]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for hit in hits:
        key = display_path(hit.path, REPO_ROOT)
        counts[key] = counts.get(key, 0) + 1
    return counts


def run_ratchet(hits: list[Hit]) -> int:
    """Allow the existing debt, forbid adding to it.

    The full translation sweep spans many commits; until it lands, the gate must
    still be able to run on every build. The ratchet therefore fails only when a
    file gains non-English lines or a previously clean file loses that status.
    """
    if not BASELINE_PATH.exists():
        print(f"check_language: no baseline at {display_path(BASELINE_PATH, REPO_ROOT)} - run --write-baseline first.")
        return 1

    baseline: dict[str, int] = json.loads(BASELINE_PATH.read_text(encoding="utf-8"))
    current = census(hits)

    regressions = [(path, current[path], baseline.get(path, 0))
                   for path in sorted(current)
                   if current[path] > baseline.get(path, 0)]
    if regressions:
        for path, now, was in regressions:
            print(f"{path}: {was} -> {now} non-English lines (regression)")
        print(f"\ncheck_language: FAILED - {len(regressions)} file(s) gained non-English prose. "
              f"Write new comments in English.")
        return 1

    improved = sum(baseline.get(p, 0) - current.get(p, 0) for p in baseline)
    total = sum(current.values())
    cleared = [p for p in baseline if p not in current]
    note = f" ({improved} lines translated, {len(cleared)} file(s) now clean)" if improved else ""
    print(f"check_language: ratchet OK - {total} non-English lines remain in {len(current)} files{note}.")
    return 0


def main() -> int:
    # Reported lines carry mathematical notation and non-Latin letters. The
    # Windows console defaults to a legacy codepage that cannot encode them, so
    # force UTF-8 and degrade gracefully rather than dying inside a build gate.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")

    ap = argparse.ArgumentParser(description="Fail on non-English prose in KATAI sources.")
    ap.add_argument("paths", nargs="*", default=None, help="files or directories (default: repo source roots)")
    ap.add_argument("--report", action="store_true", help="print a census and always exit 0")
    ap.add_argument("--by-file", action="store_true", help="with --report, one line per file instead of per hit")
    ap.add_argument("--limit", type=int, default=40, help="max hits printed in gate mode")
    ap.add_argument("--ratchet", action="store_true",
                    help="fail only on regressions against scripts/language_baseline.json")
    ap.add_argument("--write-baseline", action="store_true", help="record the current census as the ratchet baseline")
    args = ap.parse_args()

    policy = load_policy()
    skip = [re.compile(p) for p in policy.get("exempt_paths", [])]
    exempt_line = [re.compile(p) for p in policy.get("exempt_line_patterns", [])]

    roots = [Path(p).resolve() for p in args.paths] if args.paths else [
        REPO_ROOT / "kernel",
        REPO_ROOT / "tests",
        REPO_ROOT / "scripts",
        REPO_ROOT / "docs",
        REPO_ROOT / "LaTeX",
        REPO_ROOT / "CMakeLists.txt",
    ]
    roots = [r for r in roots if r.exists()]

    hits: list[Hit] = []
    for path in iter_sources(roots, skip):
        hits.extend(scan(path, exempt_line))

    if args.report:
        if args.by_file:
            per_file: dict[Path, int] = {}
            for hit in hits:
                per_file[hit.path] = per_file.get(hit.path, 0) + 1
            for path, count in sorted(per_file.items(), key=lambda kv: -kv[1]):
                print(f"{count:5d}  {path.relative_to(REPO_ROOT).as_posix()}")
        else:
            for hit in hits:
                print(hit.render(REPO_ROOT))
        files = len({h.path for h in hits})
        print(f"\n{len(hits)} non-English lines in {files} files.")
        return 0

    if args.write_baseline:
        BASELINE_PATH.write_text(json.dumps(census(hits), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"check_language: baseline written - {len(hits)} lines in {len(census(hits))} files.")
        return 0

    if args.ratchet:
        return run_ratchet(hits)

    if not hits:
        print("check_language: clean - all source prose is English.")
        return 0

    for hit in hits[: args.limit]:
        print(hit.render(REPO_ROOT))
    if len(hits) > args.limit:
        print(f"... and {len(hits) - args.limit} more (run with --report for the full list)")
    files = len({h.path for h in hits})
    print(f"\ncheck_language: FAILED - {len(hits)} non-English lines in {files} files.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
