#pragma once
// FILE PATHS ARE UTF-8. Every path that reaches this tree as a std::string -- from the Studio's
// file dialog, the command line, a Python str -- is UTF-8, and it is opened as UTF-8 here, on every
// platform. A narrow std::ifstream(std::string) on Windows reads the bytes in the process's ANSI
// code page instead (cp1254 on a Turkish system), so until 2026-10-22 a project saved in a folder
// or under a name with letters outside ASCII (a user name with a c-cedilla, a file name with a
// dotless i) could not be opened by the program that saved it, or was written under a mangled
// name. std::filesystem::path built from char8_t is UTF-8 by definition, and the standard
// streams accept it.

#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>

namespace katai::io {

inline std::filesystem::path utf8_path(const std::string& s) {
    return std::filesystem::path(std::u8string(s.begin(), s.end()));
}

// Write `bytes` to `path` (UTF-8) so that the file at `path` is either the old one or the whole
// new one, never a part: the bytes go to a temporary file beside the target, which replaces the
// target only after it was written and closed without an error.
inline bool write_file_atomically(const std::string& path, const std::string& bytes, std::string* err) {
    const std::filesystem::path target = utf8_path(path);
    std::filesystem::path temp = target;
    temp += ".saving";
    {
        std::ofstream f(temp, std::ios::binary | std::ios::trunc);
        if (!f) { if (err) *err = "cannot open file for writing: " + path; return false; }
        f.write(bytes.data(), (std::streamsize)bytes.size());
        f.close();
        if (!f) {
            std::error_code ec;
            std::filesystem::remove(temp, ec);
            if (err) *err = "write failed: " + path;
            return false;
        }
    }
    std::error_code ec;
    std::filesystem::rename(temp, target, ec);   // replaces an existing target
    if (ec) {
        std::filesystem::remove(temp, ec);
        if (err) *err = "cannot replace " + path + " (is it open in another program?)";
        return false;
    }
    return true;
}

}  // namespace katai::io
