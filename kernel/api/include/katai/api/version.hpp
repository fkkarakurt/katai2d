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
// 0.9.0's checklist is CLOSED, and closed one item short of what it opened with on purpose. The
// numerical-controls line N-3 onward -- automatic step size from the iteration band, then
// arc-length -- is NOT in it: both change how every existing answer is reached, so they belong at
// the head of a line rather than at the end of one, where the re-measurement they force can be
// made deliberately (V0.9.0-FOUNDATIONS-AUDIT.md 3.2, 3.3). What the gate does close is in
// CHANGELOG.md under [0.9.0]. The next line reopens as "-dev" with N-3 as its first item.

namespace katai::api {

inline constexpr const char* kVersion = "0.9.0";
inline constexpr const char* kVersionDate = "2026-08-31";

// "KATAI 2D <kVersion>" -- the canonical short identity for titles and report headers.
inline constexpr const char* kAppName = "KATAI 2D";

}  // namespace katai::api
