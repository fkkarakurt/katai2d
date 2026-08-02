# Contributing to KATAI 2D

Thank you for considering it. This project has one non-negotiable rule and a
few practical ones; all of them exist to protect the same thing — **numbers
someone may build on must be right**.

## The non-negotiable: verification

A capability that is not verified is not in the product. Every change that
touches physics must come with a test checked against one of:

1. a **closed-form solution**, stated in full in the test;
2. an **independent computation path** that shares no code with the solver;
3. a **published benchmark**, with the primary source cited.

Self-consistency ("the code agrees with itself") does not count. New
verified capabilities declare their reference in the test's declaration
block — the verification matrix and its bibliography are generated from
those blocks (`python scripts/check_references.py --write`) and a suite gate
fails when they drift.

## Building and testing

```powershell
.\scripts\build.ps1 -Configure -Test     # configure + build + the full suite
```

Every test must pass before and after your change — the suite, not a
subset, is the specification. The source-tree gates run inside it:
`check_language` (English-only prose), `check_architecture` (the layer
contract), `check_k2d_spec` (writer ↔ spec ↔ schema), and the reference
registry. Run `.\scripts\check_composition.ps1` when your change could
affect what the build is made of.

## Practical rules

- **English only** in code, comments, messages and documentation.
- **Honest refusals**: when an input or a state cannot be handled
  correctly, refuse with a clear message. Never approximate silently — a
  plausible wrong number is worse than a stop.
- **Units are fixed** (kN, m, day, degrees) and pass through verbatim.
- Match the style of the code around you; keep one logical change per pull
  request.
- The `.k2d` format is a contract: a field you add must appear in the
  writer, the reader, `docs/k2d-format.md` and `docs/k2d.schema.json` — the
  `check_k2d_spec` gate holds you to it.

## Sign-off (DCO)

Contributions are accepted under the
[Developer Certificate of Origin](https://developercertificate.org/). Add a
`Signed-off-by: Your Name <email>` line to your commits (`git commit -s`) to
certify you have the right to submit the work under the project's
Apache-2.0 license.

## Reporting problems

Open an issue with the smallest input that reproduces the problem — a
`.k2d` file or a short Python script is ideal. If a number looks wrong,
say what you expected and why (a source is gold). "It refused my input" may
be working as intended; the message text should make that judgeable.
