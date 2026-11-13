#pragma once

// Light travel time and the refractive index at the observation level. See docs/refraction.md.

#include "obdeect/propagation/profile.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace obdeect::propagation {

inline constexpr double speed_of_light_m_per_ns = 0.299792458;

namespace detail {

inline constexpr std::size_t ambient_grid_points = 40000;

} // namespace detail

// Refractive index n = 1 + (n-1) at the observation level from level heights [m] and n-1.
// (n-1) is interpolated as the logarithm with a natural cubic spline, evaluated on 40 000
// equidistant heights from the lowest to the highest level and linear in between. Below the
// first level the first value is used, at or above the last level n = 1. The height must lie
// below the second-highest level (the spline needs a level above it), otherwise it is rejected.
inline double ambient_group_index(const std::vector<double> &heights_m,
                                  const std::vector<double> &index_minus_one,
                                  double observation_m) {
  const std::size_t n = heights_m.size();
  if (n < 3 || index_minus_one.size() != n) {
    throw std::invalid_argument("ambient index: needs at least 3 levels of equal length");
  }
  if (!std::isfinite(observation_m)) {
    throw std::invalid_argument("ambient index: height must be finite");
  }
  std::vector<double> logs(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(heights_m[i]) || (i > 0 && !(heights_m[i] > heights_m[i - 1]))) {
      throw std::invalid_argument("ambient index: heights must be finite and strictly ascending");
    }
    if (!(index_minus_one[i] > 0.) || !std::isfinite(index_minus_one[i])) {
      throw std::invalid_argument("ambient index: n-1 must be positive and finite");
    }
    logs[i] = std::log(index_minus_one[i]);
  }
  if (observation_m >= heights_m.back()) {
    return 1.;
  }
  if (observation_m < heights_m.front()) {
    return 1. + index_minus_one.front();
  }
  if (observation_m >= heights_m[n - 2]) {
    throw std::invalid_argument("ambient index: height must lie below the second-highest level");
  }
  const double bottom = heights_m.front();
  const double step = (heights_m.back() - bottom) / double(detail::ambient_grid_points - 1);
  std::size_t k = static_cast<std::size_t>((observation_m - bottom) / step);
  k = std::min(k, detail::ambient_grid_points - 2);
  const double z0 = bottom + step * double(k), z1 = bottom + step * double(k + 1);
  const detail::NaturalSpline spline(heights_m, logs);
  const double y0 = spline(z0);
  const double y1 = spline(z1);
  const double fraction = (observation_m - z0) / (z1 - z0);
  return 1. + std::exp(y0 + fraction * (y1 - y0));
}

inline double ambient_group_index(const AtmosphereProfile &profile, double observation_m) {
  std::vector<double> heights(profile.levels());
  for (std::size_t i = 0; i < heights.size(); ++i) {
    heights[i] = profile.table_altitude_cm()[i] / 100.;
  }
  return ambient_group_index(heights, profile.table_index_minus_one(), observation_m);
}

// Time to travel path_m through air of index n.
inline double propagation_time_ns(double path_m, double index) {
  if (!std::isfinite(path_m) || !std::isfinite(index)) {
    throw std::invalid_argument("propagation time needs a finite path and index");
  }
  return path_m * index / speed_of_light_m_per_ns;
}

// Extra travel time of a straight path at zenith angle zenith [rad] between two heights [m]
// compared with vacuum: (1/c) integral of (n-1) dh / cos(zenith).
inline double straight_path_air_delay_ns(const AtmosphereProfile &profile, double low_m,
                                         double high_m, double zenith) {
  if (!std::isfinite(low_m) || !std::isfinite(high_m) || !std::isfinite(zenith)) {
    throw std::invalid_argument("air delay needs finite heights and zenith angle");
  }
  if (high_m < low_m) {
    throw std::invalid_argument("air delay needs low <= high");
  }
  if (zenith < 0. || zenith >= 1.5707963267948966) {
    throw std::invalid_argument("air delay needs 0 <= zenith < pi/2");
  }
  const double integral_cm =
      profile.index_integral_cm(100. * low_m) - profile.index_integral_cm(100. * high_m);
  return integral_cm / 100. / (speed_of_light_m_per_ns * std::cos(zenith));
}

} // namespace obdeect::propagation
