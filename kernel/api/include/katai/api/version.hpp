#pragma once
// KATAI 2D version identity -- the ONE source of truth. Every user-facing surface (window title,
// splash, About, text/HTML report headers, `katai info`, katai.__version__, the wheel) reads THESE
// constants; a hardcoded version string anywhere else is a bug (two of them already went stale at
// "0.3-beta" once). Traceability is a V&V requirement, not cosmetics: a printed report must say
// which build produced it, or a number on someone's desk cannot be matched to the validation
// record that covers it. The identity lives in the published facade because every front end needs
// it -- a version only the GUI could read was the bug in the making this move removes.
//
// Scheme: MAJOR.MINOR.PATCH per the release gates in docs/internal/ROADMAP.md sec 4 ("-dev" while the
// gate's checklist is open; drop it in the gate-closing commit, which also updates kVersionDate).

namespace katai::api {

inline constexpr const char* kVersion = "0.6.1";
inline constexpr const char* kVersionDate = "2026-08-07";

// "KATAI 2D 0.6.0-dev" -- the canonical short identity for titles and report headers.
inline constexpr const char* kAppName = "KATAI 2D";

}  // namespace katai::api
