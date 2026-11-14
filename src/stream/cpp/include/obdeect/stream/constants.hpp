#pragma once

// Physical and mathematical constants shared by all obdeect components.

#include <numbers>

namespace obdeect::stream {

inline constexpr double pi = std::numbers::pi;
inline constexpr double two_pi = 2.0 * std::numbers::pi;

// Speed of light in vacuum [m/ns].
inline constexpr double speed_of_light_m_per_ns = 0.299792458;

// Earth radius used by CORSIKA's CURVED option (6371.315 km) [m].
inline constexpr double earth_radius_m = 6371315.0;

constexpr double deg_to_rad(double degrees) { return degrees * pi / 180.0; }
constexpr double rad_to_deg(double radians) { return radians * 180.0 / pi; }

} // namespace obdeect::stream
