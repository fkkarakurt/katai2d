#pragma once
// 2D geometry (minimal) — Phase 0's single primitive: a rectangular domain with a
// material. Real CAD (points/lines/regions, polygon booleans) arrives in Phase 1.

namespace katai::geometry {

// Axis-aligned rectangular domain. (x0, y0) is the bottom-left corner.
struct RectangularDomain {
    double x0 = 0.0;
    double y0 = 0.0;
    double width = 0.0;
    double height = 0.0;
    int material_id = 0;

    double x_max() const { return x0 + width; }
    double y_max() const { return y0 + height; }
    double area() const { return width * height; }
};

} // namespace katai::geometry
