"""katai.liquefaction against the AFAD worked example (TBDY 2018 Annex 16B; KV-LIQ-001).

The kernel test checks the chain step by step; this one checks that the scripting surface gives
an engineer the same answer from a boring log: the factor at every evaluated depth and the
screening of the others.
"""
import sys

import katai.liquefaction as liq

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


# Slide 44: water table 1.80 m, gamma 17 / 20, gamma_w 9.8, S_DS 0.75, Mw 7.5; C_R as printed.
CR = {3.00: 0.85, 4.50: 0.85, 6.00: 0.95, 9.00: 1.00, 10.50: 1.00, 12.00: 1.00}
rows = liq.boring(
    tests=[(3.00, 16, 5), (4.50, 39, 16), (6.00, 20, 16), (9.00, 18, 22), (10.50, 16, 32),
           (12.00, 20, 89)],
    water_depth=1.80, gamma=17.0, gamma_sat=20.0, gamma_w=9.8, SDS=0.75, Mw=7.5,
    rod_factor=lambda z: CR[z], susceptible=lambda z: z != 12.00)
by = {r["depth"]: r for r in rows}
for z, fs in ((3.00, 0.90), (6.00, 1.05), (9.00, 0.83), (10.50, 0.79)):
    r = by[z]
    print("   %.2f m: N1,60f %.2f, FS %.4f, %s" % (z, r["N160f"], r["FS"], r["verdict"]))
    check(abs(r["FS"] - fs) <= 0.005 and r["verdict"] == "triggered",
          "%.2f m: factor %.2f, liquefies" % (z, fs))
check(by[4.50]["verdict"] == "not required" and abs(by[4.50]["N160"] - 42.5) <= 0.05,
      "4.50 m: N1,60 42.5 >= 30, not required (16.6.5)")
check(by[12.00]["verdict"] == "not susceptible", "12.00 m: not a liquefiable soil (16.6.2)")
check(abs(liq.rod_factor_tbdy(7.0) - 0.95) < 1e-12 and abs(liq.rod_factor_tbdy(10.0) - 1.0) < 1e-12,
      "Table 16B.1: C_R 0.95 for a 7 m rod, 1.00 from 10 m")

if failures:
    print("\n%d check(s) failed" % failures)
    sys.exit(1)
print("\nOK: katai.liquefaction reproduces the AFAD worked example")
