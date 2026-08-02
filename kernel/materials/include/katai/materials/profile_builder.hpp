#pragma once
// Depth-varying stiffness / strength profile construction (Stage B6: extracted
// from the application driver). PLAXIS E'_inc / c'_inc about y_ref become a
// MaterialProfile, gated by what the resolved model actually reads -- the
// catalogue's profile_E / profile_c flags (registry.hpp), stated once with the
// model instead of re-derived at every caller. A model that cannot honour a
// gradient keeps a uniform() profile, and its callers take the bit-for-bit
// constant-E path.

#include <katai/materials/material_model.hpp>
#include <katai/materials/registry.hpp>

namespace katai::core {

inline MaterialProfile build_profile(const ModelEntry& entry,
                                     double E_inc, double c_inc, double y_ref) {
    MaterialProfile p;
    if (entry.profile_E) p.E_inc = E_inc;
    if (entry.profile_c) p.c_inc = c_inc;
    p.y_ref = y_ref;
    return p;
}

} // namespace katai::core
