# elcentro-1940-ns.dat — provenance

**Record:** Imperial Valley earthquake, 1940-05-18 (UTC May 19, 04:36), El Centro Terminal
Substation station, **North-South (S00E) component** — the literature's classic reference
accelerogram.

**File:** two columns `time [s]  acceleration [g]`, dt = 0.02 s, **1560 samples**
(t = 0 → 31.18 s). Downloaded from: vibrationdata.com/elcentro.dat (Tom Irvine,
vibrationdata — a widely used digitization), 2026-07-20. The file is kept AS-IS
(integrity); parsing happens in the test.

**Identity verification (measured on the day of download):** PGA = **0.31882 g @
t = 2.02 s** — matches the published classic characteristic value (≈0.319 g, peak ~2 s)
exactly. A ±2-3% spread between digitizations is normal (USGS corrected series, NGA RSN6
etc. carry slightly different samples); this file's role in V&V is not "the single
official truth" but **a real record with verified identity**. The tests check the PGA
against the published value and the spectrum shape against published wide bands
(test_real_record.cpp).
