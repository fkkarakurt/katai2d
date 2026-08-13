// The SSI structural-dynamics bodies (elastic stiffness, consistent mass, self-weight,
// seismic influence vector), compiled ONCE (section 5.2 batch 3d). The embedded-beam
// coupling is templated on the soil element type and instantiates here for tri6/tri15;
// consumers see declarations only.
#include <katai/analysis/structural_dynamics.hpp>

namespace katai::core {

namespace detail {

// The EQUATION index of a plate translation DOF. An INDEPENDENT wall (embedded wall) gives
// trans_dof[.]>=0 -> that extra DOF; otherwise the plate translation DOF is SHARED with the soil
// node (bonded plate-in-soil). The SAME rule as the plate assembly in nonlinear_solver.cpp (if
// the two ever diverge, test_ssi_dynamics' cross-verification catches it).
inline int plate_trans_eq(const DofMap& dofs, int trans_gdof, int node, int comp) {
    return dofs.equation(trans_gdof >= 0 ? trans_gdof : dofs.global_dof(node, comp));
}

// Scatter a symmetric element matrix onto the free DOFs (a fixed DOF, eq<0, is dropped).
template <int N, class Mat>
void scatter(const std::array<int, N>& eq, const Mat& Ke, math::SparseMatrixBuilder& builder) {
    for (int a = 0; a < N; ++a) {
        if (eq[a] < 0) continue;
        for (int b = 0; b < N; ++b)
            if (eq[b] >= 0) builder.add_entry(eq[a], eq[b], Ke(a, b));
    }
}

// Embedded beam (pile row) -- ELASTIC contribution. TEMPLATED on the soil element type: the
// skin/foot coupling is mesh-INCOMPATIBLE (a beam point sits INSIDE a soil element at some
// (xi,eta)) -> the soil shape functions N_s are needed. Hence the dispatch runs off
// mesh.nodes_per_element (the assembler.cpp pattern). The u=0 elastic branch of
// nonlinear_solver.cpp's embedded-beam assembly: the beam via plate::stiffness (already
// elastic), skin D = k_a(t x t) + k_n(n x n) (NO axial cap T_max), foot D = D_foot(t x t)
// (k_n=0).
template <class E>
void assemble_embedded_beams(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::vector<ebeam::EmbeddedBeam>& beams,
                             math::SparseMatrixBuilder& builder) {
    // One coupling point: an (eq, cx, cy) list + local D -> scatter wJ.c'.D.c (the elastic form
    // of axial_couple).
    auto couple = [&](const int* eq, const double* cx, const double* cy, int nc,
                      const Eigen::Vector2d& tang, double k_a, double k_n, double wJ) {
        const Eigen::Vector2d nrm(-tang(1), tang(0));
        const Eigen::Matrix2d D = k_a * (tang * tang.transpose()) + k_n * (nrm * nrm.transpose());
        for (int d = 0; d < nc; ++d) {
            if (eq[d] < 0) continue;
            const Eigen::RowVector2d cd(cx[d], cy[d]);
            for (int e2 = 0; e2 < nc; ++e2)
                if (eq[e2] >= 0)
                    builder.add_entry(eq[d], eq[e2],
                                      wJ * (cd * D * Eigen::Vector2d(cx[e2], cy[e2]))(0, 0));
        }
    };
    for (const auto& eb : beams) {
        for (const auto& el : eb.elements) {          // (1) the beam (Timoshenko), on its own extra DOFs
            plate::NodeCoords Xe;
            std::array<int, 9> eq;
            for (int k = 0; k < 3; ++k) {
                Xe(k, 0) = eb.node_x[el[k]]; Xe(k, 1) = eb.node_y[el[k]];   // the beam's OWN coordinates
                eq[3 * k + 0] = ebeam::trans_eq(eb, el[k], 0, dofs);
                eq[3 * k + 1] = ebeam::trans_eq(eb, el[k], 1, dofs);
                eq[3 * k + 2] = dofs.equation(eb.dof_phi[el[k]]);
            }
            detail::scatter<9>(eq, plate::stiffness(Xe, eb.props), builder);
        }
        constexpr int NC = 6 + 2 * E::kNodeCount;
        for (const auto& sp : eb.skin) {              // (2) skin: beam <-> soil (mesh-incompatible)
            if (!sp.ok) continue;
            const auto Ns = E::shape_functions(sp.xi_s, sp.eta_s);
            std::array<int, NC> eq; std::array<double, NC> cx, cy; int nc = 0;
            for (int i = 0; i < 3; ++i) {
                eq[nc] = ebeam::trans_eq(eb, sp.beam_node[i], 0, dofs); cx[nc] = sp.Nb(i); cy[nc] = 0.0; ++nc;
                eq[nc] = ebeam::trans_eq(eb, sp.beam_node[i], 1, dofs); cx[nc] = 0.0; cy[nc] = sp.Nb(i); ++nc;
            }
            for (int j = 0; j < E::kNodeCount; ++j) {
                const int sn = mesh.node_of(sp.soil_elem, j);
                eq[nc] = dofs.equation(dofs.global_dof(sn, 0)); cx[nc] = -Ns(j); cy[nc] = 0.0; ++nc;
                eq[nc] = dofs.equation(dofs.global_dof(sn, 1)); cx[nc] = 0.0; cy[nc] = -Ns(j); ++nc;
            }
            couple(eq.data(), cx.data(), cy.data(), nc, sp.tang, sp.k_a, sp.k_n, sp.wJ);
        }
        if (eb.foot.D_foot > 0.0 && eb.foot.ok) {     // (3) foot (axial spring; wJ=1, k_n=0)
            const auto Ns = E::shape_functions(eb.foot.xi_s, eb.foot.eta_s);
            constexpr int NF = 2 + 2 * E::kNodeCount;
            std::array<int, NF> eq; std::array<double, NF> cx, cy; int nc = 0;
            eq[nc] = ebeam::trans_eq(eb, eb.foot.beam_node, 0, dofs); cx[nc] = 1.0; cy[nc] = 0.0; ++nc;
            eq[nc] = ebeam::trans_eq(eb, eb.foot.beam_node, 1, dofs); cx[nc] = 0.0; cy[nc] = 1.0; ++nc;
            for (int j = 0; j < E::kNodeCount; ++j) {
                const int sn = mesh.node_of(eb.foot.soil_elem, j);
                eq[nc] = dofs.equation(dofs.global_dof(sn, 0)); cx[nc] = -Ns(j); cy[nc] = 0.0; ++nc;
                eq[nc] = dofs.equation(dofs.global_dof(sn, 1)); cx[nc] = 0.0; cy[nc] = -Ns(j); ++nc;
            }
            couple(eq.data(), cx.data(), cy.data(), nc, eb.foot.tang, eb.foot.D_foot, 0.0, 1.0);
        }
    }
}

}  // namespace detail

void assemble_structural_stiffness(const mesh::Mesh& mesh, const DofMap& dofs,
                                          const Structures& structures,
                                          math::SparseMatrixBuilder& builder) {
    // --- Plate (3-node, tri6 edge): 9 DOF [u_x,u_y,phi]x3. Translation shared with the soil or
    // an independent wall DOF; rotation is the plate's own extra DOF. The plate is ELASTIC ->
    // the tangent = K_p always.
    for (const auto& pe : structures.plates) {
        plate::NodeCoords Xe;
        std::array<int, 9> eq;
        for (int k = 0; k < 3; ++k) {
            Xe(k, 0) = mesh.x[pe.nodes[k]];
            Xe(k, 1) = mesh.y[pe.nodes[k]];
            eq[3 * k + 0] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 0], pe.nodes[k], 0);
            eq[3 * k + 1] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 1], pe.nodes[k], 1);
            eq[3 * k + 2] = dofs.equation(pe.rot_dof[k]);
        }
        detail::scatter<9>(eq, plate::stiffness(Xe, pe.props), builder);
    }

    // --- 5-node (quartic) plate (tri15 edge): 15 DOF; the same logic as the 3-node one.
    for (const auto& pe : structures.plates5) {
        plate::NodeCoords5 Xe;
        std::array<int, 15> eq;
        for (int k = 0; k < 5; ++k) {
            Xe(k, 0) = mesh.x[pe.nodes[k]];
            Xe(k, 1) = mesh.y[pe.nodes[k]];
            eq[3 * k + 0] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 0], pe.nodes[k], 0);
            eq[3 * k + 1] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 1], pe.nodes[k], 1);
            eq[3 * k + 2] = dofs.equation(pe.rot_dof[k]);
        }
        detail::scatter<15>(eq, plate::stiffness5(Xe, pe.props), builder);
    }

    // --- Anchor (axial spring, ELASTIC branch): k = EA/L, U = (u_b-u_a).dir -> K = k.(g x g).
    // fixed-end (node_b<0): the far end is fixed in the RELATIVE frame = a support moving with
    // the base (the rock/soil-anchored root idealisation; recorded in the formulation document
    // sec 10.3).
    for (const auto& an : structures.anchors) {
        const Eigen::Vector2d Xa(mesh.x[an.node_a], mesh.y[an.node_a]);
        const Eigen::Vector2d Xb = an.node_b >= 0
            ? Eigen::Vector2d(mesh.x[an.node_b], mesh.y[an.node_b]) : an.fixed_point;
        const Eigen::Vector2d dvec = Xb - Xa;
        const double Lgeom = dvec.norm();
        if (Lgeom < 1e-30) continue;
        const Eigen::Vector2d dir = dvec / Lgeom;
        const double kk = an.EA / (an.L > 0.0 ? an.L : Lgeom);
        const std::array<int, 4> eq = {dofs.equation(dofs.global_dof(an.node_a, 0)),
                                       dofs.equation(dofs.global_dof(an.node_a, 1)),
                                       an.node_b >= 0 ? dofs.equation(dofs.global_dof(an.node_b, 0)) : -1,
                                       an.node_b >= 0 ? dofs.equation(dofs.global_dof(an.node_b, 1)) : -1};
        const double g[4] = {-dir(0), -dir(1), dir(0), dir(1)};   // ∂U/∂u
        Eigen::Matrix<double, 4, 4> Ke;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) Ke(i, j) = kk * g[i] * g[j];
        detail::scatter<4>(eq, Ke, builder);
    }

    // --- Geogrid (axial membrane, ELASTIC branch): K = int EA B'B ds, 2-point Gauss (weight 1).
    {
        const auto gxi = geogrid::gauss_xi();
        for (const auto& ge : structures.geogrids) {
            geogrid::NodeCoords Xe;
            std::array<int, 6> eq;
            for (int k = 0; k < 3; ++k) {
                Xe(k, 0) = mesh.x[ge.nodes[k]];
                Xe(k, 1) = mesh.y[ge.nodes[k]];
                eq[2 * k + 0] = dofs.equation(dofs.global_dof(ge.nodes[k], 0));
                eq[2 * k + 1] = dofs.equation(dofs.global_dof(ge.nodes[k], 1));
            }
            Eigen::Matrix<double, 6, 6> Ke = Eigen::Matrix<double, 6, 6>::Zero();
            for (int q = 0; q < geogrid::kGaussCount; ++q) {
                const auto kin = geogrid::axial_kin(Xe, gxi[q]);
                Ke += (kin.J * ge.props.EA) * kin.Be.transpose() * kin.Be;   // 2-point weight = 1
            }
            detail::scatter<6>(eq, Ke, builder);
        }
    }

    // --- Interface (zero thickness, ELASTIC branch): Newton-Cotes node pairs (Day & Potts) ->
    // 4 DOF per point (soil x,y + structure x,y). B rows: s=[-c,-s,+c,+s], n=[+s,-c,-s,+c].
    // Elastic tangent: D = k_s.(a x a) + k_n.(b x b) (the Coulomb return stays on the elastic
    // branch at zero relative displacement -> the same tangent; no slip in v1 dynamics, see the
    // constitutive limit in the header intro).
    // 3-node interface.
    {
        const auto ncpts = iface::nc_points();
        for (const auto& ie : structures.interfaces) {
            iface::NodeCoords Xe;
            for (int k = 0; k < 3; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
            for (int q = 0; q < iface::kPointCount; ++q) {
                const int nd = ncpts[q].node;
                const auto fr = iface::edge_frame(Xe, ncpts[q].xi);
                const double c = fr.c, s = fr.s, wJ = ncpts[q].w * fr.J;
                const std::array<int, 4> eq = {dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 0)),
                                               dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 1)),
                                               dofs.equation(ie.struct_dof[2 * nd + 0]),
                                               dofs.equation(ie.struct_dof[2 * nd + 1])};
                const double a[4] = {-c, -s, c, s};
                const double b[4] = {s, -c, -s, c};
                Eigen::Matrix<double, 4, 4> Ke;
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        Ke(i, j) = wJ * (ie.props.ks * a[i] * a[j] + ie.props.kn * b[i] * b[j]);
                detail::scatter<4>(eq, Ke, builder);
            }
        }
    }
    // 5-node interface (tri15 edge).
    {
        const auto ncpts5 = iface::nc_points5();
        for (const auto& ie : structures.interfaces5) {
            iface::NodeCoords5 Xe;
            for (int k = 0; k < 5; ++k) { Xe(k, 0) = mesh.x[ie.soil_nodes[k]]; Xe(k, 1) = mesh.y[ie.soil_nodes[k]]; }
            for (int q = 0; q < iface::kPointCount5; ++q) {
                const int nd = ncpts5[q].node;
                const auto fr = iface::edge_frame5(Xe, ncpts5[q].xi);
                const double c = fr.c, s = fr.s, wJ = ncpts5[q].w * fr.J;
                const std::array<int, 4> eq = {dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 0)),
                                               dofs.equation(dofs.global_dof(ie.soil_nodes[nd], 1)),
                                               dofs.equation(ie.struct_dof[2 * nd + 0]),
                                               dofs.equation(ie.struct_dof[2 * nd + 1])};
                const double a[4] = {-c, -s, c, s};
                const double b[4] = {s, -c, -s, c};
                Eigen::Matrix<double, 4, 4> Ke;
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        Ke(i, j) = wJ * (ie.props.ks * a[i] * a[j] + ie.props.kn * b[i] * b[j]);
                detail::scatter<4>(eq, Ke, builder);
            }
        }
    }
    // --- Embedded beam (pile row). Templated on the soil element type (mesh-incompatible skin)
    // -> dispatch.
    if (!structures.embedded_beams.empty()) {
        if (mesh.nodes_per_element == Tri15Element::kNodeCount)
            detail::assemble_embedded_beams<Tri15Element>(mesh, dofs, structures.embedded_beams, builder);
        else
            detail::assemble_embedded_beams<Tri6Element>(mesh, dofs, structures.embedded_beams, builder);
    }
}

void assemble_structural_mass(const mesh::Mesh& mesh, const DofMap& dofs,
                                     const Structures& structures,
                                     math::SparseMatrixBuilder& builder) {
    // --- Plate (3-node): m_ij = int rho N_i N_j ds. N quadratic -> N_iN_j degree 4; 3-point
    // Gauss (exact to degree 5) full integration.
    {
        constexpr double g3 = 0.7745966692414834;   // sqrt(3/5)
        const std::array<double, 3> xi3{-g3, 0.0, g3};
        const std::array<double, 3> w3{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
        for (const auto& pe : structures.plates) {
            if (pe.props.rho_A <= 0.0 && pe.props.rho_I <= 0.0) continue;
            plate::NodeCoords Xe;
            std::array<int, 9> eq;
            for (int k = 0; k < 3; ++k) {
                Xe(k, 0) = mesh.x[pe.nodes[k]];
                Xe(k, 1) = mesh.y[pe.nodes[k]];
                eq[3 * k + 0] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 0], pe.nodes[k], 0);
                eq[3 * k + 1] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 1], pe.nodes[k], 1);
                eq[3 * k + 2] = dofs.equation(pe.rot_dof[k]);
            }
            Eigen::Matrix<double, 9, 9> Me = Eigen::Matrix<double, 9, 9>::Zero();
            for (int q = 0; q < 3; ++q) {
                const auto e = plate::detail::edge_kin(Xe, xi3[q]);
                const double wds = w3[q] * e.J;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) {
                        const double NN = wds * e.N(i) * e.N(j);
                        Me(3 * i + 0, 3 * j + 0) += pe.props.rho_A * NN;   // u_x
                        Me(3 * i + 1, 3 * j + 1) += pe.props.rho_A * NN;   // u_y
                        Me(3 * i + 2, 3 * j + 2) += pe.props.rho_I * NN;   // phi (rotary inertia)
                    }
            }
            detail::scatter<9>(eq, Me, builder);
        }
    }
    // --- 5-node plate: N quartic -> N_iN_j degree 8; 5-point Gauss (exact to degree 9).
    {
        constexpr double a5 = 0.5384693101056831, b5 = 0.9061798459386640;
        const std::array<double, 5> xi5{-b5, -a5, 0.0, a5, b5};
        const std::array<double, 5> w5{0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                                       0.4786286704993665, 0.2369268850561891};
        for (const auto& pe : structures.plates5) {
            if (pe.props.rho_A <= 0.0 && pe.props.rho_I <= 0.0) continue;
            plate::NodeCoords5 Xe;
            std::array<int, 15> eq;
            for (int k = 0; k < 5; ++k) {
                Xe(k, 0) = mesh.x[pe.nodes[k]];
                Xe(k, 1) = mesh.y[pe.nodes[k]];
                eq[3 * k + 0] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 0], pe.nodes[k], 0);
                eq[3 * k + 1] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 1], pe.nodes[k], 1);
                eq[3 * k + 2] = dofs.equation(pe.rot_dof[k]);
            }
            Eigen::Matrix<double, 15, 15> Me = Eigen::Matrix<double, 15, 15>::Zero();
            for (int q = 0; q < 5; ++q) {
                const auto e = plate::detail::edge_kin5(Xe, xi5[q]);
                const double wds = w5[q] * e.J;
                for (int i = 0; i < 5; ++i)
                    for (int j = 0; j < 5; ++j) {
                        const double NN = wds * e.N(i) * e.N(j);
                        Me(3 * i + 0, 3 * j + 0) += pe.props.rho_A * NN;
                        Me(3 * i + 1, 3 * j + 1) += pe.props.rho_A * NN;
                        Me(3 * i + 2, 3 * j + 2) += pe.props.rho_I * NN;
                    }
            }
            detail::scatter<15>(eq, Me, builder);
        }
    }
    // --- Embedded beam: the beam's OWN node coordinates + OWN extra DOFs (not mesh nodes).
    // The same consistent mass as the 3-node plate; rho_A/rho_I from eb.props (the pile's weight
    // -> the driver). A massless pile is as wrong in seismics as a massless wall (stiffness
    // without inertia).
    {
        constexpr double g3 = 0.7745966692414834;
        const std::array<double, 3> xi3{-g3, 0.0, g3};
        const std::array<double, 3> w3{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
        for (const auto& eb : structures.embedded_beams) {
            if (eb.props.rho_A <= 0.0 && eb.props.rho_I <= 0.0) continue;
            for (const auto& el : eb.elements) {
                plate::NodeCoords Xe;
                std::array<int, 9> eq;
                for (int k = 0; k < 3; ++k) {
                    Xe(k, 0) = eb.node_x[el[k]]; Xe(k, 1) = eb.node_y[el[k]];
                    eq[3 * k + 0] = ebeam::trans_eq(eb, el[k], 0, dofs);
                    eq[3 * k + 1] = ebeam::trans_eq(eb, el[k], 1, dofs);
                    eq[3 * k + 2] = dofs.equation(eb.dof_phi[el[k]]);
                }
                Eigen::Matrix<double, 9, 9> Me = Eigen::Matrix<double, 9, 9>::Zero();
                for (int q = 0; q < 3; ++q) {
                    const auto e = plate::detail::edge_kin(Xe, xi3[q]);
                    const double wds = w3[q] * e.J;
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j) {
                            const double NN = wds * e.N(i) * e.N(j);
                            Me(3 * i + 0, 3 * j + 0) += eb.props.rho_A * NN;
                            Me(3 * i + 1, 3 * j + 1) += eb.props.rho_A * NN;
                            Me(3 * i + 2, 3 * j + 2) += eb.props.rho_I * NN;
                        }
                }
                detail::scatter<9>(eq, Me, builder);
            }
        }
    }
}

void assemble_structural_weight(const mesh::Mesh& mesh, const DofMap& dofs,
                                       const Structures& structures, double g,
                                       Eigen::VectorXd& f) {
    // --- Plate (3-node): int N_i ds, 3-point Gauss (more than exact for quadratic N x J).
    {
        constexpr double g3 = 0.7745966692414834;   // sqrt(3/5)
        const std::array<double, 3> xi3{-g3, 0.0, g3};
        const std::array<double, 3> w3{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
        for (const auto& pe : structures.plates) {
            if (pe.props.rho_A <= 0.0) continue;
            plate::NodeCoords Xe;
            std::array<int, 3> eqy;
            for (int k = 0; k < 3; ++k) {
                Xe(k, 0) = mesh.x[pe.nodes[k]];
                Xe(k, 1) = mesh.y[pe.nodes[k]];
                eqy[k] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 1], pe.nodes[k], 1);
            }
            for (int q = 0; q < 3; ++q) {
                const auto e = plate::detail::edge_kin(Xe, xi3[q]);
                const double wds = w3[q] * e.J * pe.props.rho_A * g;
                for (int i = 0; i < 3; ++i)
                    if (eqy[i] >= 0) f(eqy[i]) -= wds * e.N(i);
            }
        }
    }
    // --- 5-node plate: the same, 5-point Gauss.
    {
        constexpr double a5 = 0.5384693101056831, b5 = 0.9061798459386640;
        const std::array<double, 5> xi5{-b5, -a5, 0.0, a5, b5};
        const std::array<double, 5> w5{0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                                       0.4786286704993665, 0.2369268850561891};
        for (const auto& pe : structures.plates5) {
            if (pe.props.rho_A <= 0.0) continue;
            plate::NodeCoords5 Xe;
            std::array<int, 5> eqy;
            for (int k = 0; k < 5; ++k) {
                Xe(k, 0) = mesh.x[pe.nodes[k]];
                Xe(k, 1) = mesh.y[pe.nodes[k]];
                eqy[k] = detail::plate_trans_eq(dofs, pe.trans_dof[2 * k + 1], pe.nodes[k], 1);
            }
            for (int q = 0; q < 5; ++q) {
                const auto e = plate::detail::edge_kin5(Xe, xi5[q]);
                const double wds = w5[q] * e.J * pe.props.rho_A * g;
                for (int i = 0; i < 5; ++i)
                    if (eqy[i] >= 0) f(eqy[i]) -= wds * e.N(i);
            }
        }
    }
    // --- Embedded beam (pile row): the beam's OWN coordinates + OWN u_y extra DOFs.
    {
        constexpr double g3 = 0.7745966692414834;
        const std::array<double, 3> xi3{-g3, 0.0, g3};
        const std::array<double, 3> w3{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
        for (const auto& eb : structures.embedded_beams) {
            if (eb.props.rho_A <= 0.0) continue;
            for (const auto& el : eb.elements) {
                plate::NodeCoords Xe;
                std::array<int, 3> eqy;
                for (int k = 0; k < 3; ++k) {
                    Xe(k, 0) = eb.node_x[el[k]];
                    Xe(k, 1) = eb.node_y[el[k]];
                    eqy[k] = ebeam::trans_eq(eb, el[k], 1, dofs);
                }
                for (int q = 0; q < 3; ++q) {
                    const auto e = plate::detail::edge_kin(Xe, xi3[q]);
                    const double wds = w3[q] * e.J * eb.props.rho_A * g;
                    for (int i = 0; i < 3; ++i)
                        if (eqy[i] >= 0) f(eqy[i]) -= wds * e.N(i);
                }
            }
        }
    }
}

Eigen::VectorXd seismic_influence_x(const mesh::Mesh& mesh, const DofMap& dofs,
                                           const Structures& structures) {
    Eigen::VectorXd r = Eigen::VectorXd::Zero(dofs.equation_count());
    for (int n = 0; n < mesh.node_count; ++n) {                 // soil (and twin seam) nodes
        const int eq = dofs.equation(dofs.global_dof(n, 0));
        if (eq >= 0) r[eq] = 1.0;
    }
    auto mark = [&](int gdof) {                                 // independent wall translation DOF
        if (gdof < 0) return;
        const int eq = dofs.equation(gdof);
        if (eq >= 0) r[eq] = 1.0;
    };
    for (const auto& pe : structures.plates)  for (int k = 0; k < 3; ++k) mark(pe.trans_dof[2 * k + 0]);
    for (const auto& pe : structures.plates5) for (int k = 0; k < 5; ++k) mark(pe.trans_dof[2 * k + 0]);
    // The embedded beam's TRANSLATION DOFs are its own extra DOFs too (NOT mesh nodes) -> a
    // rigid base translation moves them as well. Skipping them means the pile takes NO inertial
    // force from the base = silently wrong SSI (the K.r=0 rigid-mode check catches it).
    for (const auto& eb : structures.embedded_beams)
        for (int d : eb.dof_x) mark(d);
    return r;
}

}  // namespace katai::core
