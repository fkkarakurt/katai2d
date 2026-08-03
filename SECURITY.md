# Security policy

## Supported versions

KATAI 2D is pre-1.0: security fixes land on `main` and ship with the next
release. The latest release is the only supported version.

## What counts as a security bug

`katai` runs locally, makes no network connections and collects nothing. The
main attack surface is file parsing: a crafted `.k2d` project file or `.res`
results file that makes the reader crash, hang, exhaust memory or execute
anything — in short, anything worse than the documented refusal (exit codes
3 and 4) — is a security bug. So is anything that lets `install.ps1` install
something other than the published, hash-verified release.

A wrong *number* is not handled here: numerical defects are verification
issues and are worked in the open, against primary sources — use the
"Numerical result question" issue template for those.

## How to report

Report vulnerabilities privately through GitHub's
[report a vulnerability](https://github.com/fkkarakurt/katai2d/security/advisories/new)
form. If you cannot use it, mail <fatihkucukkarakurt@gmail.com> with
`katai2d security` in the subject line.

Include what a good bug report would include: the smallest input that
demonstrates the problem, and what it does. You will get an acknowledgement
within 7 days. Please keep the report private until a fix is released; you
will be credited in the advisory unless you prefer otherwise.
