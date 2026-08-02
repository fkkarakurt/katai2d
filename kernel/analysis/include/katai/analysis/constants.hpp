#pragma once
// Physical and reporting constants of the analysis layer (Stage B2: extracted
// from the application driver -- the engine cannot depend on the application
// for the unit weight of water).

namespace katai::core {

// Unit weight of water [kN/m^3].
inline constexpr double kGammaWater = 9.81;

// Gravity [m/s^2]. Weight -> mass: soil rho = gamma/g [Mg/m^3], plate rho_A = w/g [Mg/m].
inline constexpr double kGravity = 9.81;

// Cap for the interface Coulomb demand/capacity ratio. Where the joint goes into TENSION its shear
// capacity tau_max is exactly 0 (it has separated), so |tau|/tau_max is unbounded -- report a capped
// number rather than an infinity that would poison plots, tables and the .res file.
inline constexpr double kUtilCap = 999.0;

} // namespace katai::core
