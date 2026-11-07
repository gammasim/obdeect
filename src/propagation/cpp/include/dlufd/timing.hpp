#pragma once

#include "dlufd/profile.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace dlufd {

// Speed of light in vacuum [m/ns], as obdeect vec.hpp kSpeedOfLightMPerNs.
inline constexpr double speed_of_light_m_per_ns = 0.299792458;

// obdeect propagation.ambient_group_index(): the refractive index n = 1 + (n-1) at the
// observation level. ln(n-1) is a natural cubic spline in height [m]; it is evaluated on 40000
// equidistant heights from the lowest to the highest level, linearly interpolated in ln(n-1).
// Requires at least 4 levels, strictly ascending heights and positive n-1.
inline double ambient_group_index(std::span<const double> altitude_m,
                                  std::span<const double> index_minus_one,
                                  double observation_level_m) {
  const std::size_t count = altitude_m.size();
  if (count != index_minus_one.size() || count < 4)
    throw std::invalid_argument("ambient index needs at least 4 matching levels");
  for (std::size_t i = 0; i < count; ++i) {
    if (!std::isfinite(altitude_m[i]) || !std::isfinite(index_minus_one[i]) ||
        !(index_minus_one[i] > 0.))
      throw std::invalid_argument("ambient index needs finite positive n-1");
    if (i > 0 && !(altitude_m[i] > altitude_m[i - 1]))
      throw std::invalid_argument("ambient index heights must be strictly ascending");
  }
  if (!std::isfinite(observation_level_m))
    throw std::invalid_argument("observation level must be finite");
  if (observation_level_m <= altitude_m.front())
    return 1. + index_minus_one.front();
  if (observation_level_m >= altitude_m.back())
    return 1.;
  if (observation_level_m > altitude_m[count - 2])
    throw std::invalid_argument("observation level lies in the extrapolated upper layer");

  std::vector<double> log_n_1(count);
  for (std::size_t i = 0; i < count; ++i)
    log_n_1[i] = std::log(index_minus_one[i]);
  const CubicSpline spline(std::vector<double>(altitude_m.begin(), altitude_m.end()), log_n_1);

  const double step = (altitude_m.back() - altitude_m.front()) / 39999.;
  const double coordinate = (observation_level_m - altitude_m.front()) / step;
  const auto lower = static_cast<std::size_t>(coordinate);
  const double fraction = coordinate - static_cast<double>(lower);
  const double first = spline.polynomial(altitude_m.front() + static_cast<double>(lower) * step);
  const double second =
      spline.polynomial(altitude_m.front() + static_cast<double>(lower + 1) * step);
  return 1. + std::exp(first + fraction * (second - first));
}

// Group index at the observation level for an atmospheric table.
inline double ambient_group_index(const AtmosphereProfile &profile, double observation_level_m) {
  return ambient_group_index(profile.table_altitude_m(), profile.table_index_minus_one(),
                             observation_level_m);
}

// Light travel time over a path length in a medium of group index n [ns] (obdeect trace.hpp:
// path * propagation_group_index / c).
// Throws std::invalid_argument for a non-finite path or index.
inline double propagation_time_ns(double path_m, double group_index) {
  require_finite(path_m, "path length");
  require_finite(group_index, "group index");
  return path_m * group_index / speed_of_light_m_per_ns;
}

// Extra travel time of light along a straight path in the atmosphere between two heights
// (vacuum-vs-air): (1/c) * integral of (n-1) ds, with ds = dh / cos(zenith) (plane-parallel).
// The integral is exact for the piecewise-linear n-1 of the profile table. A positive value
// means the light is slower than in vacuum over the same geometric path.
inline double straight_path_air_delay_ns(const AtmosphereProfile &profile, double height_low_m,
                                         double height_high_m, double zenith_rad) {
  require_finite(height_low_m, "height_low_m");
  require_finite(height_high_m, "height_high_m");
  require_finite(zenith_rad, "zenith angle");
  if (!(height_high_m >= height_low_m))
    throw std::invalid_argument("height_high_m must not be below height_low_m");
  if (!(zenith_rad >= 0. && zenith_rad < 1.5707963267948966))
    throw std::invalid_argument("zenith angle must be in [0, pi/2)");
  const double integral_cm = profile.index_integral_cm(100. * height_high_m) -
                             profile.index_integral_cm(100. * height_low_m);
  return integral_cm / 100. / (speed_of_light_m_per_ns * std::cos(zenith_rad));
}

} // namespace dlufd
