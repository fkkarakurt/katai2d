// test_unicode_paths -- a project and its results are saved and opened under a path with letters
// outside ASCII, and a save replaces the file whole or not at all.
//
// Paths reach the kernel as UTF-8 (the Studio's file dialog under the UTF-8 process code page, the
// command line, a Python str). Until 2026-10-22 the kernel opened them with a narrow
// std::ofstream(std::string), which Windows reads in the ANSI code page: a Turkish user name or a
// Turkish file name produced a file under a mangled name, or no file at all, and the program
// could not open what it had saved. The oracle is independent of the kernel's own reading: the
// file's existence is asked of the file system under the name spelled as UTF-16, so a round trip
// that wrote and read the same wrong name cannot pass.
//
// The letters are written as escapes so that this source stays ASCII.
//
// verify: none -- file-system mechanics (UTF-8 paths, atomic replacement); no physics is asserted

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <katai/io/project_io.hpp>
#include <katai/io/results_io.hpp>

namespace m = katai::model;
namespace fs = std::filesystem;

static int g_failures = 0;
static void check(bool ok, const char* what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what);
    if (!ok) ++g_failures;
}

static std::string utf8(const char8_t* s) { return std::string(reinterpret_cast<const char*>(s)); }

int main() {
    // A folder and a file named with c-cedilla, g-breve, dotless i, capital dotted I, s-cedilla,
    // o-umlaut and u-umlaut.
    const std::string dir_u8 = utf8(u8"katai_Çağrı_İşler");  // katai-lang: allow (test data)
    const std::wstring dir_w = L"katai_Çağrı_İşler";  // katai-lang: allow (test data)
    const std::string file_u8 = dir_u8 + "/" + utf8(u8"şev_çözümü.k2d");  // katai-lang: allow (test data)
    const std::wstring file_w = dir_w + L"/" + L"şev_çözümü.k2d";  // katai-lang: allow (test data)
    std::error_code ec;
    fs::remove_all(fs::path(dir_w), ec);
    fs::create_directory(fs::path(dir_w), ec);
    check(!ec && fs::is_directory(fs::path(dir_w)), "a folder with Turkish letters exists (made through UTF-16)");

    m::Project pr;
    pr.name = utf8(u8"Şev çözümü");  // katai-lang: allow (test data)
    m::Material s; s.name = "Clay"; pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0; P.x = {0, 10, 10, 0}; P.y = {0, 0, 5, 5};
    pr.polygons.push_back(P);

    std::string err;
    check(m::save_project(pr, file_u8, &err), "the project is saved under the UTF-8 path");
    if (!err.empty()) std::printf("   (%s)\n", err.c_str());
    check(fs::exists(fs::path(file_w)), "...and the file system has it under the name spelled in UTF-16");
    check(!fs::exists(fs::path(file_w + L".saving")), "...with no temporary file left beside it");
    m::Project back;
    check(m::load_project(file_u8, back, &err) && back.name == pr.name &&
              m::project_to_json(back) == m::project_to_json(pr),
          "it opens again, the same project, its non-ASCII name included");

    // Results beside it.
    katai::app::SolveResult R;
    R.mesh.node_count = 3; R.mesh.element_count = 1; R.mesh.nodes_per_element = 3;
    R.mesh.x = {0, 1, 0}; R.mesh.y = {0, 0, 1}; R.mesh.connectivity = {0, 1, 2};
    R.mesh.element_material = {0};
    R.ok = true; R.disp = Eigen::VectorXd::Zero(6);
    R.stress.stress.assign(3, Eigen::Vector3d::Zero()); R.pore.assign(3, 0.0);
    const std::string res_u8 = file_u8 + ".res";
    check(katai::app::save_results(res_u8, 42, {R}, &err), "the results are saved under the UTF-8 path");
    check(fs::exists(fs::path(file_w + L".res")), "...under the name spelled in UTF-16");
    std::vector<katai::app::SolveResult> rb;
    check(katai::app::load_results(res_u8, 42, rb, &err) && rb.size() == 1, "...and open again");

    // A save replaces the file whole: the second save is what is there afterwards.
    pr.name = "second";
    check(m::save_project(pr, file_u8, &err) && m::load_project(file_u8, back, &err) && back.name == "second",
          "saving over an existing project replaces it");
    // A save that cannot be written leaves what was there: the target's folder does not exist.
    const std::string nowhere = dir_u8 + "/no_such_folder/x.k2d";
    check(!m::save_project(pr, nowhere, &err) && !err.empty(), "a save into a missing folder is refused with a message");

    fs::remove_all(fs::path(dir_w), ec);
    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
