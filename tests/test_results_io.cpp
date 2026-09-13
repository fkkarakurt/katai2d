// Binary results file round-trip (results_io.hpp). FEM results are big numeric arrays -- they
// live in a compact binary sibling of the JSON model file (a double costs 8 bytes here vs ~20 as
// JSON text, with zero parse cost). Checks: exact round-trip of every field across two phases
// (raw IEEE bytes -> bitwise equality), the model-hash STALE guard, and corrupt/truncated
// rejection without a crash.
#include <katai/io/results_io.hpp>

#include <cstdio>
#include <cstdlib>

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

katai::mesh::Mesh tiny_mesh() {
    katai::mesh::Mesh m;
    m.node_count = 6; m.element_count = 1; m.nodes_per_element = 6;
    m.x = {0, 1, 0, 0.5, 0.5, 0};
    m.y = {0, 0, 1, 0, 0.5, 0.5};
    m.connectivity = {0, 1, 2, 3, 4, 5};
    m.element_material = {0};
    m.bottom_nodes = {0, 3, 1};
    m.boundary_nodes = {0, 1, 2, 3, 4, 5};
    return m;
}

std::vector<katai::app::SolveResult> build_phases() {
    std::vector<katai::app::SolveResult> ph(2);
    const katai::mesh::Mesh m = tiny_mesh();
    for (int k = 0; k < 2; ++k) {
        auto& R = ph[k];
        R.mesh = m;
        R.ok = true;
        R.nil_step = k == 1;
        R.max_disp = 0.0123 + k;
        R.load_factor = 1.0;
        R.fos = k == 1 ? 1.43 : -1.0;
        R.message = k == 0 ? "initial \"phase\"" : "excavate";
        // v6: the reason a solve stopped short. A load factor without it cannot be read -- the
        // same 62% is a capacity when a mechanism formed and a setting when the iteration budget
        // ran out -- so it has to survive the round trip like every other qualification.
        R.stopped_by = k == 1 ? katai::core::NewtonResult::Abandonment::IterationBudget
                              : katai::core::NewtonResult::Abandonment::None;
        // v7: WHICH convergence criteria the accepted iterate satisfied. Same argument as
        // stopped_by one line up: "converged" is a claim about the accuracy of the numbers in
        // this file, and a reopened result that drops the evidence keeps the claim.
        R.convergence.measured = true;
        R.convergence.csp = k == 1 ? 0.0976543 : 1.0;
        R.convergence.force_error = 3.705e-7 * (k + 1);
        R.convergence.moment_error = k == 1 ? 4.25e-4 : 0.0;
        R.convergence.has_moment = k == 1;
        R.convergence.tolerated = 1e-6;
        R.convergence.plastic_points = 347 * k;
        R.convergence.plastic_inaccurate = 337 * k;
        R.convergence.elastic_points = 1471;
        R.convergence.nl_elastic_points = 96 * k;
        R.convergence.nl_elastic_inaccurate = 5 * k;
        R.convergence.worst_plastic_error = 2.334e-2 * k;
        R.convergence.worst_nl_elastic_error = 7.69e-5 * k;
        // v8: the structural half of the same family.
        R.convergence.iface_points = 47 * k;
        R.convergence.iface_inaccurate = 14 * k;
        R.convergence.worst_iface_error = 2.05e-6 * k;
        R.convergence.feet = k;
        R.convergence.foot_force_error = 5.58e-9 * k;
        // v10: the ratio the step stopped on. A reopened result that kept only the
        // stiffness-weighted error reported a number nothing had been compared against.
        R.convergence.global_error = 6.0917253e-8 * (k + 1);
        // v9: WHAT a consolidation phase was asked to end on, and what it reached. A result that
        // says "90% consolidated" without saying which 90% (a pressure ratio, not the settlement
        // ratio) is the number without its definition -- the same reason stopped_by is stored.
        R.consol_stop = k == 1 ? katai::core::ConsolidationStop::DegreeOfConsolidation
                               : katai::core::ConsolidationStop::TimeInterval;
        R.consol_stop_met = k == 1;
        R.consol_pore_reference = 9.99833 * k;
        R.consol_degree_reached = 0.90025 * k;
        R.disp = Eigen::VectorXd::LinSpaced(12, -0.5 + k, 0.75);
        R.stress.stress.resize(6);
        for (int n = 0; n < 6; ++n)
            R.stress.stress[n] = Eigen::Vector3d(-10.5 * n - k, -20.25 * n, 3.125 * n);
        R.pore = {0, 1.5, 3, 4.5, 6, 7.5};
        if (k == 1) R.active = {0};   // excavated element
    }
    katai::app::StructForce sf;
    sf.name = "Wall"; sf.kind = 0; sf.yielded = true;
    sf.max_N = 12.5; sf.max_Q = 3.25; sf.max_M = 44.0;
    sf.stations.push_back({0.0, 10.0, 4.0, -12.5, 3.25, 0.0});
    sf.stations.push_back({2.5, 10.0, 6.5, -6.0, 1.0, 44.0});
    ph[1].struct_forces.push_back(sf);
    katai::app::InterfaceResult ir;
    ir.name = "Wall interface"; ir.any_slip = true;
    ir.slip_checked = true;   // v5: the Coulomb branch really ran (static / nonlinear dynamic)
    ir.max_abs_tau = 29.5; ir.max_abs_sigma_n = 161.0; ir.max_abs_slip = 4.1e-3;
    katai::core::InterfaceStation is0; is0.s = 0.0; is0.x = 10.0; is0.y = 2.0;
    is0.tau = 12.0; is0.sigma_n = -80.0; is0.slip = 1.0e-3; is0.gap = -2.0e-4; is0.slipping = false;
    katai::core::InterfaceStation is1; is1.s = 3.0; is1.x = 10.0; is1.y = 5.0;
    is1.tau = 29.5; is1.sigma_n = -161.0; is1.slip = 4.1e-3; is1.gap = 1.0e-4; is1.slipping = true;
    ir.stations.push_back(is0); ir.stations.push_back(is1);
    ph[1].interface_forces.push_back(ir);
    return ph;
}

bool same(const katai::app::SolveResult& a, const katai::app::SolveResult& b) {
    if (a.ok != b.ok || a.nil_step != b.nil_step || a.message != b.message) return false;
    if (a.max_disp != b.max_disp || a.load_factor != b.load_factor || a.fos != b.fos) return false;
    if (a.stopped_by != b.stopped_by) return false;
    const auto& ca = a.convergence; const auto& cb = b.convergence;
    if (ca.measured != cb.measured || ca.csp != cb.csp || ca.force_error != cb.force_error ||
        ca.moment_error != cb.moment_error || ca.has_moment != cb.has_moment ||
        ca.tolerated != cb.tolerated || ca.plastic_points != cb.plastic_points ||
        ca.plastic_inaccurate != cb.plastic_inaccurate ||
        ca.elastic_points != cb.elastic_points ||
        ca.nl_elastic_points != cb.nl_elastic_points ||
        ca.nl_elastic_inaccurate != cb.nl_elastic_inaccurate ||
        ca.worst_plastic_error != cb.worst_plastic_error ||
        ca.worst_nl_elastic_error != cb.worst_nl_elastic_error ||
        ca.iface_points != cb.iface_points || ca.iface_inaccurate != cb.iface_inaccurate ||
        ca.worst_iface_error != cb.worst_iface_error || ca.feet != cb.feet ||
        ca.foot_force_error != cb.foot_force_error) return false;
    // NaN-safe: a field that was not recorded on both sides is the same absence.
    if (!(ca.global_error == cb.global_error ||
          (ca.global_error != ca.global_error && cb.global_error != cb.global_error)))
        return false;
    if (a.consol_stop != b.consol_stop || a.consol_stop_met != b.consol_stop_met ||
        a.consol_pore_reference != b.consol_pore_reference ||
        a.consol_degree_reached != b.consol_degree_reached) return false;
    if (a.disp.size() != b.disp.size() || (a.disp - b.disp).cwiseAbs().maxCoeff() != 0.0) return false;
    if (a.stress.stress.size() != b.stress.stress.size()) return false;
    for (size_t i = 0; i < a.stress.stress.size(); ++i)
        if (a.stress.stress[i] != b.stress.stress[i]) return false;
    if (a.pore != b.pore || a.active != b.active) return false;
    if (a.mesh.x != b.mesh.x || a.mesh.y != b.mesh.y ||
        a.mesh.connectivity != b.mesh.connectivity ||
        a.mesh.element_material != b.mesh.element_material ||
        a.mesh.bottom_nodes != b.mesh.bottom_nodes ||
        a.mesh.boundary_nodes != b.mesh.boundary_nodes) return false;
    if (a.struct_forces.size() != b.struct_forces.size()) return false;
    for (size_t i = 0; i < a.struct_forces.size(); ++i) {
        const auto& x = a.struct_forces[i]; const auto& y = b.struct_forces[i];
        if (x.name != y.name || x.kind != y.kind || x.yielded != y.yielded) return false;
        if (x.max_N != y.max_N || x.max_Q != y.max_Q || x.max_M != y.max_M) return false;
        if (x.stations.size() != y.stations.size()) return false;
        for (size_t q = 0; q < x.stations.size(); ++q) {
            const auto& s = x.stations[q]; const auto& t = y.stations[q];
            if (s.s != t.s || s.x != t.x || s.y != t.y || s.N != t.N || s.Q != t.Q || s.M != t.M)
                return false;
        }
    }
    if (a.interface_forces.size() != b.interface_forces.size()) return false;
    for (size_t i = 0; i < a.interface_forces.size(); ++i) {
        const auto& x = a.interface_forces[i]; const auto& y = b.interface_forces[i];
        if (x.name != y.name || x.any_slip != y.any_slip || x.slip_checked != y.slip_checked) return false;
        if (x.max_abs_tau != y.max_abs_tau || x.max_abs_sigma_n != y.max_abs_sigma_n ||
            x.max_abs_slip != y.max_abs_slip) return false;
        if (x.stations.size() != y.stations.size()) return false;
        for (size_t q = 0; q < x.stations.size(); ++q) {
            const auto& s = x.stations[q]; const auto& t = y.stations[q];
            if (s.s != t.s || s.x != t.x || s.y != t.y || s.tau != t.tau || s.sigma_n != t.sigma_n ||
                s.slip != t.slip || s.gap != t.gap || s.slipping != t.slipping) return false;
        }
    }
    return true;
}

}  // namespace

int main() {
    std::printf("Binary results file round-trip\n\n");
    const auto phases = build_phases();
    const std::uint64_t hash = katai::app::fnv1a64("{\"katai2d\":1,...canonical model...}");
    const char* path = "test_results_tmp.res";

    std::string err;
    check(katai::app::save_results(path, hash, phases, &err), "results saved");
    if (!err.empty()) std::printf("  (%s)\n", err.c_str());

    std::vector<katai::app::SolveResult> back;
    check(katai::app::load_results(path, hash, back, &err), "results loaded (matching model hash)");
    check(back.size() == 2, "both phases restored");
    if (back.size() == 2) {
        check(same(phases[0], back[0]), "phase 0 bitwise identical (mesh, fields, message)");
        check(same(phases[1], back[1]), "phase 1 bitwise identical (active mask, FoS, diagrams)");
    }

    // STALE guard: a different model hash must reject the results with an honest message.
    err.clear();
    std::vector<katai::app::SolveResult> stale;
    check(!katai::app::load_results(path, hash + 1, stale, &err) && !err.empty(),
          "stale results rejected (model changed)");
    std::printf("   (%s)\n", err.c_str());

    // Corruption: truncate the file -> rejected, no crash, no partial garbage.
    {
        std::ifstream in(path, std::ios::binary);
        std::ostringstream ss; ss << in.rdbuf();
        const std::string data = ss.str();
        std::printf("   file size: %d bytes (binary)\n", (int)data.size());
        std::ofstream out(path, std::ios::binary);
        out.write(data.data(), (std::streamsize)(data.size() / 2));
    }
    err.clear();
    std::vector<katai::app::SolveResult> cut;
    check(!katai::app::load_results(path, hash, cut, &err) && !err.empty(),
          "truncated file rejected");

    // A version 9 file is the same bytes without the v10 field, the last double of each phase
    // record. Reading one must not invent the ratio the step stopped on -- it comes back NaN,
    // "not recorded" -- and every other field must come back exactly as before.
    {
        err.clear();
        check(katai::app::save_results(path, hash, phases, &err), "results saved again");
        std::string data;
        {
            std::ifstream in(path, std::ios::binary);
            std::ostringstream ss; ss << in.rdbuf();
            data = ss.str();
        }
        bool located = true;
        for (const auto& P : phases) {
            const double v = P.convergence.global_error;
            const std::string bytes(reinterpret_cast<const char*>(&v), sizeof v);
            const std::size_t at = data.find(bytes);
            located = located && at != std::string::npos &&
                      data.find(bytes, at + 1) == std::string::npos;
            if (at != std::string::npos) data.erase(at, sizeof v);
        }
        check(located, "the v10 field is found exactly once per phase");
        const std::uint32_t v9 = 9;
        std::memcpy(&data[4], &v9, sizeof v9);   // "K2DR", then the version
        {
            std::ofstream out(path, std::ios::binary);
            out.write(data.data(), (std::streamsize)data.size());
        }
        std::vector<katai::app::SolveResult> old;
        err.clear();
        check(katai::app::load_results(path, hash, old, &err) && old.size() == 2,
              "a version 9 file still loads");
        if (old.size() == 2) {
            check(!old[0].convergence.global_recorded() && !old[1].convergence.global_recorded() &&
                      !old[1].convergence.global_ok(),
                  "...and the ratio it never stored reads as NOT RECORDED, not as a met zero");
            bool rest_same = true;
            for (int k = 0; k < 2; ++k) {
                auto expect = phases[k];
                expect.convergence.global_error = old[k].convergence.global_error;
                rest_same = rest_same && same(expect, old[k]);
            }
            check(rest_same, "...and every other field is unchanged");
        }
    }
    std::remove(path);

    if (g_failures == 0) {
        std::printf("\nOK: binary results round-trip exact + stale/corrupt rejection\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
