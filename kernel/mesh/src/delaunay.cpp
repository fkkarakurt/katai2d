#include <katai/mesh/delaunay.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <unordered_set>
#include <utility>

#include <katai/mesh/geom_predicates.hpp>

namespace katai::mesh {
namespace {

// Triangle with CCW vertices v[0..2] and neighbours n[k] across the edge
// opposite v[k] (i.e. edge (v[(k+1)%3], v[(k+2)%3])); -1 if none. Dead triangles
// (replaced during insertion) are flagged and skipped.
struct Tri {
    std::array<int, 3> v;
    std::array<int, 3> n;
    bool alive = true;
};

class Delaunay {
public:
    Delaunay(const std::vector<double>& px, const std::vector<double>& py)
        : x_(px), y_(py), n_input_(static_cast<int>(px.size())) {}

    Triangulation build() {
        triangulate_points();
        return extract();
    }

    Triangulation build_constrained(
        const std::vector<std::array<int, 2>>& segments) {
        triangulate_points();
        for (const auto& s : segments) recover_segment(s[0], s[1]);
        restore_delaunay();  // re-Delaunay everything except constrained edges
        return extract();
    }

private:
    void triangulate_points() {
        add_super_triangle();
        const std::vector<int> order = morton_order();
        int hint = 0;
        for (int idx : order) {
            const int located = locate(hint, idx);
            hint = insert_point(idx, located).hint;
        }
    }

    double orient(int a, int b, int c) const {
        return orient2d(x_[a], y_[a], x_[b], y_[b], x_[c], y_[c]);
    }
    double in_circle(int a, int b, int c, int d) const {
        return incircle(x_[a], y_[a], x_[b], y_[b], x_[c], y_[c], x_[d], y_[d]);
    }

    int add_tri(int a, int b, int c) {
        tris_.push_back(Tri{{a, b, c}, {-1, -1, -1}, true});
        return static_cast<int>(tris_.size()) - 1;
    }

    void add_super_triangle() {
        double minx = x_[0], maxx = x_[0], miny = y_[0], maxy = y_[0];
        for (int i = 1; i < n_input_; ++i) {
            minx = std::min(minx, x_[i]); maxx = std::max(maxx, x_[i]);
            miny = std::min(miny, y_[i]); maxy = std::max(maxy, y_[i]);
        }
        const double span = std::max({maxx - minx, maxy - miny, 1.0});
        const double cx = 0.5 * (minx + maxx), cy = 0.5 * (miny + maxy);
        // The super-triangle must lie outside every real Delaunay circumcircle;
        // thin near-hull triangles can have very large circumradii, so place it
        // extremely far. The robust (filtered + dd) predicates keep the large-
        // coordinate in-circle tests exact in sign.
        const double d = 1.0e8 * span;
        super_[0] = n_input_;
        super_[1] = n_input_ + 1;
        super_[2] = n_input_ + 2;
        x_.push_back(cx - 2.0 * d); y_.push_back(cy - d);
        x_.push_back(cx + 2.0 * d); y_.push_back(cy - d);
        x_.push_back(cx);           y_.push_back(cy + 2.0 * d);
        add_tri(super_[0], super_[1], super_[2]);  // CCW by construction
    }

    // The three super-triangle vertices occupy indices [n_input_, n_input_+3);
    // Steiner points (refinement) are appended afterwards and are NOT super.
    bool is_super(int vertex) const {
        return vertex >= n_input_ && vertex < n_input_ + 3;
    }

    // Walk from a hint triangle to the triangle containing point p (= vertex pi).
    int locate(int hint, int pi) const {
        int t = (hint >= 0 && tris_[hint].alive) ? hint : any_alive();
        const int cap = 2 * static_cast<int>(tris_.size()) + 16;
        for (int step = 0; step < cap; ++step) {
            const Tri& tr = tris_[t];
            bool inside = true;
            for (int k = 0; k < 3; ++k) {
                const int a = tr.v[(k + 1) % 3], b = tr.v[(k + 2) % 3];
                if (orient(a, b, pi) < 0.0) {  // p lies outside edge k
                    t = tr.n[k];
                    inside = false;
                    break;
                }
            }
            if (inside) return t;
            if (t < 0) break;  // should not happen inside the super-triangle
        }
        return linear_locate(pi);  // robustness fallback
    }

    int any_alive() const {
        for (int t = 0; t < static_cast<int>(tris_.size()); ++t)
            if (tris_[t].alive) return t;
        return 0;
    }

    int linear_locate(int pi) const {
        for (int t = 0; t < static_cast<int>(tris_.size()); ++t) {
            if (!tris_[t].alive) continue;
            const Tri& tr = tris_[t];
            if (orient(tr.v[0], tr.v[1], pi) >= 0.0 &&
                orient(tr.v[1], tr.v[2], pi) >= 0.0 &&
                orient(tr.v[2], tr.v[0], pi) >= 0.0)
                return t;
        }
        return any_alive();
    }

    int neighbour_slot(int t, int nb) const {
        for (int k = 0; k < 3; ++k)
            if (tris_[t].n[k] == nb) return k;
        return -1;
    }

    // Outcome of a point insertion. block_a >= 0 means the insertion was refused
    // because the point's Delaunay cavity would have crossed the constrained edge
    // (block_a, block_b) -- the caller should split that subsegment instead.
    struct InsertResult {
        int hint = -1;
        int block_a = -1;
        int block_b = -1;
    };

    // Bowyer-Watson cavity insertion: delete every triangle whose circumcircle
    // contains p, then re-triangulate the star-shaped cavity from p. With
    // detect_crossing, abort (without modifying the mesh) if the cavity would
    // cross a constrained edge, returning that edge for the caller to split.
    InsertResult insert_point(int pi, int seed, bool detect_crossing = false) {
        // 1. Flood-fill the cavity (triangles with p inside their circumcircle).
        //    A version stamp marks cavity membership without an O(T) clear.
        if (stamp_.size() < tris_.size()) stamp_.resize(tris_.size(), 0);
        const int mark = ++cur_stamp_;
        cavity_.clear();
        std::vector<int> stack = {seed};
        stamp_[seed] = mark;
        while (!stack.empty()) {
            const int cur = stack.back();
            stack.pop_back();
            cavity_.push_back(cur);
            for (int k = 0; k < 3; ++k) {
                const int nb = tris_[cur].n[k];
                if (nb < 0 || stamp_[nb] == mark) continue;
                const int va = tris_[cur].v[(k + 1) % 3];
                const int vb = tris_[cur].v[(k + 2) % 3];
                const Tri& t = tris_[nb];
                const bool inside = in_circle(t.v[0], t.v[1], t.v[2], pi) > 0.0;
                if (is_constrained(va, vb)) {
                    // A constrained edge bounds the cavity (keeps segments intact).
                    // If the cavity actually wants to cross it, the point is unsafe
                    // to insert -- report the subsegment so the caller can split it.
                    if (detect_crossing && inside) return {-1, va, vb};
                    continue;
                }
                if (inside) {
                    stamp_[nb] = mark;
                    stack.push_back(nb);
                }
            }
        }

        // 2. For each cavity boundary edge, spawn a new triangle (a, b, p).
        int first_new = -1;
        edge_lo_.clear();
        edge_hi_.clear();
        edge_tri_.clear();
        edge_slot_.clear();
        for (int cur : cavity_) {
            for (int k = 0; k < 3; ++k) {
                const int nb = tris_[cur].n[k];
                if (nb >= 0 && stamp_[nb] == mark) continue;  // interior cavity edge
                const int a = tris_[cur].v[(k + 1) % 3];
                const int b = tris_[cur].v[(k + 2) % 3];
                const int nt = add_tri(a, b, pi);  // CCW: p is interior side
                tris_[nt].n[2] = nb;               // neighbour across edge (a,b)
                if (nb >= 0) {
                    const int s = neighbour_slot(nb, cur);
                    if (s >= 0) tris_[nb].n[s] = nt;
                }
                // Stitch the two p-incident edges to adjacent new triangles:
                //   slot 0 -> directed edge (b, p);  slot 1 -> directed edge (p, a).
                stitch(b, pi, nt, 0);
                stitch(pi, a, nt, 1);
                if (first_new < 0) first_new = nt;
            }
        }

        // 3. Retire the cavity triangles.
        for (int cur : cavity_) tris_[cur].alive = false;
        return {first_new, -1, -1};
    }

    // Link directed edge (u, v) of (tri, slot) with its previously seen reverse.
    void stitch(int u, int v, int tri, int slot) {
        for (std::size_t i = 0; i < edge_lo_.size(); ++i) {
            if (edge_lo_[i] == v && edge_hi_[i] == u) {  // reverse edge found
                tris_[tri].n[slot] = edge_tri_[i];
                tris_[edge_tri_[i]].n[edge_slot_[i]] = tri;
                edge_lo_[i] = edge_lo_.back(); edge_hi_[i] = edge_hi_.back();
                edge_tri_[i] = edge_tri_.back(); edge_slot_[i] = edge_slot_.back();
                edge_lo_.pop_back(); edge_hi_.pop_back();
                edge_tri_.pop_back(); edge_slot_.pop_back();
                return;
            }
        }
        edge_lo_.push_back(u); edge_hi_.push_back(v);
        edge_tri_.push_back(tri); edge_slot_.push_back(slot);
    }

    // Sort input points along a Morton (Z-order) curve for spatial locality.
    std::vector<int> morton_order() const {
        double minx = x_[0], maxx = x_[0], miny = y_[0], maxy = y_[0];
        for (int i = 1; i < n_input_; ++i) {
            minx = std::min(minx, x_[i]); maxx = std::max(maxx, x_[i]);
            miny = std::min(miny, y_[i]); maxy = std::max(maxy, y_[i]);
        }
        const double sx = (maxx > minx) ? ((1 << 21) - 1) / (maxx - minx) : 0.0;
        const double sy = (maxy > miny) ? ((1 << 21) - 1) / (maxy - miny) : 0.0;
        std::vector<std::uint64_t> key(n_input_);
        for (int i = 0; i < n_input_; ++i) {
            const auto qx = static_cast<std::uint32_t>((x_[i] - minx) * sx);
            const auto qy = static_cast<std::uint32_t>((y_[i] - miny) * sy);
            key[i] = (spread(qx) << 1) | spread(qy);
        }
        std::vector<int> order(n_input_);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return key[a] < key[b]; });
        return order;
    }

    // Interleave the low 21 bits of v with zeros (Morton spread).
    static std::uint64_t spread(std::uint32_t v) {
        std::uint64_t x = v & 0x1fffff;
        x = (x | (x << 32)) & 0x1f00000000ffffULL;
        x = (x | (x << 16)) & 0x1f0000ff0000ffULL;
        x = (x | (x << 8)) & 0x100f00f00f00f00fULL;
        x = (x | (x << 4)) & 0x10c30c30c30c30c3ULL;
        x = (x | (x << 2)) & 0x1249249249249249ULL;
        return x;
    }

    // --- Constrained Delaunay (segment recovery) ----------------------------

    static std::uint64_t edge_key(int u, int w) {
        if (u > w) std::swap(u, w);
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u)) << 32) |
               static_cast<std::uint32_t>(w);
    }
    bool is_constrained(int u, int w) const {
        return constrained_.count(edge_key(u, w)) != 0;
    }
    static bool same_edge(int x, int y, int u, int w) {
        return (x == u && y == w) || (x == w && y == u);
    }

    // Proper intersection of open segments (p1,p2) and (p3,p4).
    bool seg_cross(int p1, int p2, int p3, int p4) const {
        const double d1 = orient(p3, p4, p1), d2 = orient(p3, p4, p2);
        const double d3 = orient(p1, p2, p3), d4 = orient(p1, p2, p4);
        return ((d1 > 0) != (d2 > 0)) && (d1 != 0) && (d2 != 0) &&
               ((d3 > 0) != (d4 > 0)) && (d3 != 0) && (d4 != 0);
    }

    bool edge_exists(int a, int b) const {
        for (const Tri& t : tris_) {
            if (!t.alive) continue;
            bool ha = false, hb = false;
            for (int k = 0; k < 3; ++k) { ha |= t.v[k] == a; hb |= t.v[k] == b; }
            if (ha && hb) return true;
        }
        return false;
    }

    // Locate the triangle/slot whose edge opposite v[slot] equals {u,w}.
    std::pair<int, int> find_edge(int u, int w) const {
        for (int t = 0; t < static_cast<int>(tris_.size()); ++t) {
            if (!tris_[t].alive) continue;
            for (int k = 0; k < 3; ++k)
                if (same_edge(tris_[t].v[(k + 1) % 3], tris_[t].v[(k + 2) % 3], u, w))
                    return {t, k};
        }
        return {-1, -1};
    }

    // Flip the edge opposite v[k] of triangle t (shared with neighbour n[k]).
    // The diagonal (b,c) of quad (a,b,d,c) is replaced by (a,d); both triangles
    // are rewritten in place and all neighbour links are repaired.
    void flip(int t, int k) {
        const int nb = tris_[t].n[k];
        const int k2 = neighbour_slot(nb, t);
        const int a = tris_[t].v[k];
        const int b = tris_[t].v[(k + 1) % 3];
        const int c = tris_[t].v[(k + 2) % 3];
        const int d = tris_[nb].v[k2];
        const int n_ab = tris_[t].n[(k + 2) % 3];   // across (a,b)
        const int n_ca = tris_[t].n[(k + 1) % 3];   // across (c,a)
        const int n_bd = tris_[nb].n[(k2 + 1) % 3]; // across (b,d)
        const int n_dc = tris_[nb].n[(k2 + 2) % 3]; // across (d,c)

        tris_[t].v = {a, b, d};
        tris_[t].n = {n_bd, nb, n_ab};
        tris_[nb].v = {a, d, c};
        tris_[nb].n = {n_dc, n_ca, t};

        if (n_ca >= 0) tris_[n_ca].n[neighbour_slot(n_ca, t)] = nb;
        if (n_bd >= 0) tris_[n_bd].n[neighbour_slot(n_bd, nb)] = t;
        // n_ab keeps pointing at t, n_dc keeps pointing at nb.
    }

    // Gather the edges of the current triangulation crossed by segment (a,b),
    // in order, by marching across triangles from a toward b.
    bool collect_crossing(int a, int b, std::vector<std::pair<int, int>>& out) {
        int start = -1, slot = -1, fu = -1, fw = -1;
        for (int t = 0; t < static_cast<int>(tris_.size()) && start < 0; ++t) {
            if (!tris_[t].alive) continue;
            for (int ka = 0; ka < 3; ++ka) {
                if (tris_[t].v[ka] != a) continue;
                const int u = tris_[t].v[(ka + 1) % 3], w = tris_[t].v[(ka + 2) % 3];
                if (seg_cross(a, b, u, w)) { start = t; slot = ka; fu = u; fw = w; break; }
            }
        }
        if (start < 0) return false;
        out.emplace_back(fu, fw);
        int pu = fu, pw = fw;
        int t = tris_[start].n[slot];
        int guard = 0;
        const int cap = static_cast<int>(tris_.size()) + 16;
        while (t >= 0 && guard++ < cap) {
            const Tri& tr = tris_[t];
            if (tr.v[0] == b || tr.v[1] == b || tr.v[2] == b) return true;
            int next = -1, eu = -1, ew = -1;
            for (int k = 0; k < 3; ++k) {
                const int x = tr.v[(k + 1) % 3], y = tr.v[(k + 2) % 3];
                if (same_edge(x, y, pu, pw)) continue;
                if (seg_cross(a, b, x, y)) { eu = x; ew = y; next = tr.n[k]; break; }
            }
            if (eu < 0) return false;  // degenerate (vertex on segment)
            out.emplace_back(eu, ew);
            pu = eu; pw = ew;
            t = next;
        }
        return false;
    }

    // Force segment (a,b) to be an edge (Anglada): repeatedly flip crossing edges
    // whose quadrilateral is convex; deferred edges become convex after others flip.
    void recover_segment(int a, int b) {
        if (a == b) return;
        if (edge_exists(a, b)) { constrained_.insert(edge_key(a, b)); return; }
        std::vector<std::pair<int, int>> cross;
        if (!collect_crossing(a, b, cross)) return;  // degenerate input: skip
        std::size_t qi = 0;
        int guard = 0;
        const int cap = 100 * static_cast<int>(cross.size()) + 1000;
        while (qi < cross.size() && guard++ < cap) {
            const int u = cross[qi].first, w = cross[qi].second;
            ++qi;
            const auto [t, k] = find_edge(u, w);
            if (t < 0 || tris_[t].n[k] < 0) continue;
            const int nb = tris_[t].n[k];
            const int A = tris_[t].v[k];
            const int B = tris_[nb].v[neighbour_slot(nb, t)];
            const double du = orient(A, B, u), dw = orient(A, B, w);
            if (!((du > 0) != (dw > 0)) || du == 0 || dw == 0) {
                cross.emplace_back(u, w);  // non-convex quad: defer
                continue;
            }
            flip(t, k);  // (u,w) -> (A,B)
            if (seg_cross(a, b, A, B)) cross.emplace_back(A, B);
        }
        constrained_.insert(edge_key(a, b));
    }

    // Restore the Delaunay property everywhere except across constrained edges.
    void restore_delaunay() {
        bool changed = true;
        int guard = 0;
        const int cap = 30 * static_cast<int>(tris_.size()) + 1000;
        while (changed && guard++ < cap) {
            changed = false;
            for (int t = 0; t < static_cast<int>(tris_.size()); ++t) {
                if (!tris_[t].alive) continue;
                bool flipped = false;
                for (int k = 0; k < 3 && !flipped; ++k) {
                    const int nb = tris_[t].n[k];
                    if (nb < 0) continue;
                    const int u = tris_[t].v[(k + 1) % 3], w = tris_[t].v[(k + 2) % 3];
                    if (is_constrained(u, w)) continue;
                    const int B = tris_[nb].v[neighbour_slot(nb, t)];
                    if (in_circle(tris_[t].v[0], tris_[t].v[1], tris_[t].v[2], B) <= 0.0)
                        continue;
                    const int A = tris_[t].v[k];
                    const double du = orient(A, B, u), dw = orient(A, B, w);
                    if (((du > 0) != (dw > 0)) && du != 0 && dw != 0) {
                        flip(t, k);
                        changed = flipped = true;
                    }
                }
            }
        }
    }

    Triangulation extract() const {
        Triangulation out;
        out.point_count = n_input_;
        out.x.assign(x_.begin(), x_.begin() + n_input_);  // input keeps its indices
        out.y.assign(y_.begin(), y_.begin() + n_input_);
        std::vector<int> remap(x_.size(), -1);
        for (int i = 0; i < n_input_; ++i) remap[i] = i;
        for (const Tri& t : tris_) {
            if (!t.alive) continue;
            if (is_super(t.v[0]) || is_super(t.v[1]) || is_super(t.v[2])) continue;
            // Drop triangles that fall outside a (possibly non-convex) PSLG
            // boundary -- e.g. those filling a concave notch of the convex hull.
            const double gx = (x_[t.v[0]] + x_[t.v[1]] + x_[t.v[2]]) / 3.0;
            const double gy = (y_[t.v[0]] + y_[t.v[1]] + y_[t.v[2]]) / 3.0;
            if (!in_domain(gx, gy)) continue;
            std::array<int, 3> tv;
            for (int k = 0; k < 3; ++k) {
                const int v = t.v[k];
                if (remap[v] < 0) {  // first sighting of a Steiner point
                    remap[v] = static_cast<int>(out.x.size());
                    out.x.push_back(x_[v]);
                    out.y.push_back(y_[v]);
                }
                tv[k] = remap[v];
            }
            out.triangles.push_back(tv);
        }
        return out;
    }

    // --- Ruppert quality refinement -----------------------------------------

    double dist2(int i, int j) const {
        const double dx = x_[i] - x_[j], dy = y_[i] - y_[j];
        return dx * dx + dy * dy;
    }

    std::pair<double, double> circumcenter(int a, int b, int c) const {
        const double ax = x_[a], ay = y_[a], bx = x_[b], by = y_[b];
        const double cx = x_[c], cy = y_[c];
        const double d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
        const double a2 = ax * ax + ay * ay, b2 = bx * bx + by * by,
                     c2 = cx * cx + cy * cy;
        const double ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
        const double uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
        return {ux, uy};
    }

    // A vertex p lies inside the diametral circle of segment (a,b) iff the angle
    // a-p-b is obtuse, i.e. (p-a).(p-b) < 0.
    bool inside_diametral(double px, double py, int a, int b) const {
        return (px - x_[a]) * (px - x_[b]) + (py - y_[a]) * (py - y_[b]) < 0.0;
    }

    bool segment_encroached(int a, int b) const {
        const auto [t, k] = find_edge(a, b);
        if (t < 0) return false;
        const int apex = tris_[t].v[k];
        if (!is_super(apex) && inside_diametral(x_[apex], y_[apex], a, b))
            return true;
        const int nb = tris_[t].n[k];
        if (nb >= 0) {
            const int apex2 = tris_[nb].v[neighbour_slot(nb, t)];
            if (!is_super(apex2) && inside_diametral(x_[apex2], y_[apex2], a, b))
                return true;
        }
        return false;
    }

    bool find_encroached_subsegment(int& sa, int& sb) const {
        for (const std::uint64_t key : constrained_) {
            const int a = static_cast<int>(key >> 32);
            const int b = static_cast<int>(key & 0xffffffffu);
            if (segment_encroached(a, b)) { sa = a; sb = b; return true; }
        }
        return false;
    }

    bool point_encroaches(double px, double py, int& sa, int& sb) const {
        for (const std::uint64_t key : constrained_) {
            const int a = static_cast<int>(key >> 32);
            const int b = static_cast<int>(key & 0xffffffffu);
            if (inside_diametral(px, py, a, b)) { sa = a; sb = b; return true; }
        }
        return false;
    }

    // Walk from triangle `start` toward point (px,py). Returns the containing
    // triangle, or -1 if a constrained (boundary) edge lies between -- i.e. the
    // point is outside the domain across (seg_a, seg_b), which should be split
    // instead of inserting the point there.
    int walk_to(int start, double px, double py, int& seg_a, int& seg_b) const {
        int t = (start >= 0 && tris_[start].alive) ? start : any_alive();
        const int cap = 2 * static_cast<int>(tris_.size()) + 16;
        for (int step = 0; step < cap; ++step) {
            const Tri& tr = tris_[t];
            bool inside = true;
            for (int k = 0; k < 3; ++k) {
                const int a = tr.v[(k + 1) % 3], b = tr.v[(k + 2) % 3];
                if (orient2d(x_[a], y_[a], x_[b], y_[b], px, py) < 0.0) {
                    if (is_constrained(a, b)) { seg_a = a; seg_b = b; return -1; }
                    t = tr.n[k];
                    inside = false;
                    break;
                }
            }
            if (inside) return t;
            if (t < 0) return -1;
        }
        return t;
    }

    // Flag input vertices where two boundary segments meet at < 60 deg. Splitting
    // a subsegment incident to such a corner at its midpoint cascades without
    // terminating (each split encroaches the matching subsegment on the other leg);
    // concentric-shell splitting (below) fixes it.
    void mark_acute_vertices() {
        acute_.clear();
        const double cos60 = 0.5;  // angle < 60 deg  <=>  cos(angle) > 0.5
        for (int v = 0; v < n_input_; ++v) {
            std::vector<std::pair<double, double>> dirs;
            for (const auto& s : boundary_) {
                int other = -1;
                if (s[0] == v) other = s[1];
                else if (s[1] == v) other = s[0];
                else continue;
                const double dx = x_[other] - x_[v], dy = y_[other] - y_[v];
                const double L = std::sqrt(dx * dx + dy * dy);
                if (L > 0.0) dirs.push_back({dx / L, dy / L});
            }
            for (std::size_t p = 0; p < dirs.size(); ++p)
                for (std::size_t q = p + 1; q < dirs.size(); ++q)
                    if (dirs[p].first * dirs[q].first +
                            dirs[p].second * dirs[q].second > cos60)
                        acute_.insert(v);
        }
    }

    void split_segment(int a, int b) {
        // Concentric-shell rule (Ruppert; Shewchuk 2000): if exactly one endpoint
        // is an acute input vertex v, split at a power-of-two distance from v
        // (the shell nearest the midpoint) instead of the midpoint. Successive
        // splits then halve onto the same shells on both legs of the corner, so
        // mutual encroachment balances and refinement terminates. Non-acute
        // subsegments keep the standard midpoint bisection.
        double t = 0.5;  // split fraction measured from a
        const bool aa = acute_.count(a) != 0, ab = acute_.count(b) != 0;
        if (aa != ab) {
            const double L = std::sqrt(dist2(a, b));
            const double d = std::ldexp(  // 2^round(log2(L/2)) in [0.354L, 0.707L]
                1.0, static_cast<int>(std::lround(std::log2(0.5 * L))));
            const double frac = d / L;             // fraction from the acute vertex
            t = aa ? frac : 1.0 - frac;            // acute endpoint is a or b
        }
        const double mx = x_[a] + t * (x_[b] - x_[a]);
        const double my = y_[a] + t * (y_[b] - y_[a]);
        constrained_.erase(edge_key(a, b));  // let the cavity merge across it
        x_.push_back(mx);
        y_.push_back(my);
        const int m = static_cast<int>(x_.size()) - 1;
        insert_point(m, locate(any_alive(), m));
        constrained_.insert(edge_key(a, m));
        constrained_.insert(edge_key(m, b));
    }

    // Point inside the meshing domain, tested against the actual PSLG boundary
    // polygon (boundary_) by ray casting (crossing-number / PNPOLY). Unlike a
    // convex-hull test this is correct for NON-CONVEX domains (e.g. a slope with a
    // re-entrant toe): triangles filling a concave notch -- inside the hull but
    // outside the polygon -- are correctly classified as exterior, so they are
    // neither refined (no escaped circumcentre runaway) nor emitted. For a convex
    // boundary it coincides with the hull test. When no boundary is registered
    // (raw/constrained Delaunay, not the quality mesher) every point is "inside".
    bool in_domain(double px, double py) const {
        if (outline_.empty()) return true;
        bool inside = false;
        for (const auto& s : outline_) {
            const double x0 = x_[s[0]], y0 = y_[s[0]];
            const double x1 = x_[s[1]], y1 = y_[s[1]];
            if ((y0 > py) != (y1 > py)) {
                const double xint = x0 + (py - y0) / (y1 - y0) * (x1 - x0);
                if (px < xint) inside = !inside;
            }
        }
        return inside;
    }

    // Return the first real triangle that is skinny (circumradius / shortest edge
    // exceeds the bound) or larger than the area cap -- constant max_area, or the
    // sizing field evaluated at the centroid when one is set (local refinement).
    int find_bad_triangle(double radius_edge_bound, double max_area) const {
        const double b2 = radius_edge_bound * radius_edge_bound;
        for (int t = 0; t < static_cast<int>(tris_.size()); ++t) {
            if (!tris_[t].alive) continue;
            const int a = tris_[t].v[0], b = tris_[t].v[1], c = tris_[t].v[2];
            if (is_super(a) || is_super(b) || is_super(c)) continue;
            // Only refine triangles inside the domain (exclude exterior ones, e.g.
            // those spawned by a circumcentre that escaped past the boundary).
            const double gx = (x_[a] + x_[b] + x_[c]) / 3.0;
            const double gy = (y_[a] + y_[b] + y_[c]) / 3.0;
            if (!in_domain(gx, gy)) continue;
            const double hmin2 = std::min({dist2(a, b), dist2(b, c), dist2(c, a)});
            const auto [ux, uy] = circumcenter(a, b, c);
            const double r2 = (ux - x_[a]) * (ux - x_[a]) + (uy - y_[a]) * (uy - y_[a]);
            if (r2 > b2 * hmin2) return t;
            const double cap = size_field_ ? size_field_(gx, gy) : max_area;
            if (cap > 0.0) {
                const double area = 0.5 * std::fabs(orient(a, b, c));
                if (area > cap) return t;
            }
        }
        return -1;
    }

    void refine(double min_angle_deg, double max_area) {
        const double pi = 3.14159265358979323846;
        const double bound = 1.0 / (2.0 * std::sin(min_angle_deg * pi / 180.0));
        // Refinement only touches in-domain triangles (find_bad_triangle), and a
        // circumcentre that would land outside the domain (walk_to crosses a
        // boundary) splits that boundary segment instead -- this prevents the
        // exterior runaway that an escaped circumcentre would otherwise cause. The
        // cap is a final safety valve.
        const int cap = 200000;
        for (int guard = 0; guard < cap; ++guard) {
            int sa, sb;
            if (find_encroached_subsegment(sa, sb)) {  // conform boundaries first
                split_segment(sa, sb);
                continue;
            }
            const int bad = find_bad_triangle(bound, max_area);
            if (bad < 0) return;  // quality achieved
            const auto [cx, cy] =
                circumcenter(tris_[bad].v[0], tris_[bad].v[1], tris_[bad].v[2]);
            int ea, eb;
            if (point_encroaches(cx, cy, ea, eb)) {  // would destroy a segment
                split_segment(ea, eb);
                continue;
            }
            // Reject a circumcentre that lies outside the domain: if reaching it
            // from the skinny triangle crosses a boundary segment, split that
            // segment instead (prevents exterior runaway).
            int wa = -1, wb = -1;
            const int seed = walk_to(bad, cx, cy, wa, wb);
            if (seed < 0) {
                if (wa >= 0) split_segment(wa, wb);
                continue;
            }
            x_.push_back(cx);
            y_.push_back(cy);
            const int ci = static_cast<int>(x_.size()) - 1;
            const InsertResult r = insert_point(ci, seed, true);
            if (r.block_a >= 0) {
                // The circumcentre sits across a constraint from the skinny
                // triangle; splitting that subsegment (rather than inserting the
                // circumcentre) is what removes the bad element and terminates.
                x_.pop_back();
                y_.pop_back();
                split_segment(r.block_a, r.block_b);
            }
        }
    }

public:
    Triangulation build_refined(const std::vector<std::array<int, 2>>& segments,
                                double min_angle_deg, double max_area,
                                const std::vector<std::array<int, 2>>& outline,
                                SizeField size_field = {}) {
        size_field_ = std::move(size_field);
        triangulate_points();
        for (const auto& s : segments) recover_segment(s[0], s[1]);
        restore_delaunay();
        boundary_ = segments;  // ALL constraints -> recovered + subdivided (encroachment)
        // The inside/outside test (in_domain) must use ONLY the domain OUTLINE, never the internal
        // constraint segments (plates/geogrids): an internal segment counted as a boundary would flip
        // the ray-casting parity and delete the mesh on one side of it. Empty -> fall back to all.
        outline_ = outline.empty() ? segments : outline;
        mark_acute_vertices();  // corners needing concentric-shell splitting
        refine(min_angle_deg, max_area);
        return extract();
    }

private:

    std::vector<double> x_, y_;
    int n_input_;
    std::array<int, 3> super_{};
    std::vector<Tri> tris_;

    // Scratch reused across insertions.
    std::vector<int> cavity_;
    std::vector<int> stamp_;  // version stamp for cavity membership
    int cur_stamp_ = 0;
    std::vector<int> edge_lo_, edge_hi_, edge_tri_, edge_slot_;
    std::unordered_set<std::uint64_t> constrained_;  // protected edges (CDT)
    std::vector<std::array<int, 2>> boundary_;       // ALL constraints (recovery + encroachment)
    std::vector<std::array<int, 2>> outline_;        // domain OUTLINE only (in_domain inside/outside)
    std::unordered_set<int> acute_;                  // input vertices at a < 60 deg corner
    SizeField size_field_;                           // local area cap (empty = constant max_area)
};

} // namespace

Triangulation delaunay_triangulate(const std::vector<double>& px,
                                   const std::vector<double>& py) {
    if (px.size() < 3 || px.size() != py.size()) return {};
    Delaunay d(px, py);
    return d.build();
}

Triangulation constrained_delaunay(
    const std::vector<double>& px, const std::vector<double>& py,
    const std::vector<std::array<int, 2>>& segments) {
    if (px.size() < 3 || px.size() != py.size()) return {};
    Delaunay d(px, py);
    return d.build_constrained(segments);
}

Triangulation quality_mesh(const std::vector<double>& px,
                           const std::vector<double>& py,
                           const std::vector<std::array<int, 2>>& segments,
                           double min_angle_deg, double max_area,
                           const std::vector<std::array<int, 2>>& outline) {
    if (px.size() < 3 || px.size() != py.size()) return {};
    Delaunay d(px, py);
    return d.build_refined(segments, min_angle_deg, max_area, outline);
}

Triangulation quality_mesh(const std::vector<double>& px,
                           const std::vector<double>& py,
                           const std::vector<std::array<int, 2>>& segments,
                           double min_angle_deg, const SizeField& max_area_at,
                           const std::vector<std::array<int, 2>>& outline) {
    if (px.size() < 3 || px.size() != py.size()) return {};
    Delaunay d(px, py);
    return d.build_refined(segments, min_angle_deg, 0.0, outline, max_area_at);
}

} // namespace katai::mesh
