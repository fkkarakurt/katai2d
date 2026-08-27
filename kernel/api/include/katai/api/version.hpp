#pragma once
// KATAI 2D version identity -- the ONE source of truth. Every user-facing surface (window title,
// splash, About, text/HTML report headers, `katai info`, katai.__version__, the wheel) reads THESE
// constants; a hardcoded version string anywhere else is a bug (two of them already went stale at
// "0.3-beta" once). Traceability is a V&V requirement, not cosmetics: a printed report must say
// which build produced it, or a number on someone's desk cannot be matched to the validation
// record that covers it. The identity lives in the published facade because every front end needs
// it -- a version only the GUI could read was the bug in the making this move removes.
//
// Scheme: MAJOR.MINOR.PATCH per the project's release gates ("-dev" while the
// gate's checklist is open; drop it in the gate-closing commit, which also updates kVersionDate).
// kVersionDate moves whenever kVersion does, opening a line as well as closing one: a report
// printed from a development build must not carry the previous release's date, which is the same
// traceability argument that put the identity here in the first place.
//
// 0.9.0's checklist is open. What is already in it is in CHANGELOG.md under [Unreleased]; what is
// not yet done is the numerical-controls line N-3 onward (automatic step size from the iteration
// band, then arc-length). Drop the "-dev" in the commit that closes the gate.

namespace katai::api {

inline constexpr const char* kVersion = "0.9.0-dev";
inline constexpr const char* kVersionDate = "2026-08-27";

// "KATAI 2D <kVersion>" -- the canonical short identity for titles and report headers.
inline constexpr const char* kAppName = "KATAI 2D";

}  // namespace katai::api
