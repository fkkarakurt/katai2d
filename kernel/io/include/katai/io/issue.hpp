#pragma once
// The input contract's issue vocabulary (Stage C). One Issue = one violation at
// one field path, reported the way an engineer reads a file: `path` names the
// exact JSON location (`materials[2].Rf`), `message` states the rule in plain
// engineering language and echoes the offending value with its unit. Producers:
// the reader's forward-version notes (project_io.hpp) and the validator
// (validate.hpp). Consumers: the CLI, the GUI and any script front end -- they
// print a report, they do not interpret it.
//
// This header is new code, so it lives in the convention namespace katai::io
// (the module's two legacy headers keep their historical spellings until the
// dedicated namespace-rename step; ARCHITECTURE.md records that exception).

#include <string>
#include <utility>
#include <vector>

namespace katai::io {

enum class Severity {
    Warning,   // suspicious or silently-ignored input; the solve is still defensible
    Error      // the file cannot produce a defensible solve as written
};

struct Issue {
    Severity severity = Severity::Error;
    std::string path;      // JSON field path in the .k2d file, e.g. "materials[2].Rf"
    std::string message;   // engineer-readable; echoes the offending value and unit
};

// A validation outcome: every issue in file order, with the error/warning split
// precomputed so front ends can gate on ok() without re-counting.
struct ValidationReport {
    std::vector<Issue> issues;
    int errors = 0;
    int warnings = 0;

    bool ok() const { return errors == 0; }   // warnings do not block a solve

    void add(Severity s, std::string path, std::string message) {
        (s == Severity::Error ? errors : warnings) += 1;
        issues.push_back(Issue{s, std::move(path), std::move(message)});
    }

    // One line per issue, `severity: path: message`, ready for a console or log.
    std::string to_string() const {
        std::string out;
        for (const Issue& i : issues) {
            out += (i.severity == Severity::Error) ? "error: " : "warning: ";
            out += i.path;
            out += ": ";
            out += i.message;
            out += '\n';
        }
        return out;
    }
};

}  // namespace katai::io
