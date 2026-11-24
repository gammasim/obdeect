#pragma once

// Light travel time and the refractive index at the observation level. See
// docs/propagation/refraction.md.

#include "obdeect/propagation/profile.hpp"
#include "obdeect/stream/constants.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace obdeect::propagation {

inline constexpr double speed_of_light_m_per_ns = obdeect::stream::speed_of_light_m_per_ns;

namespace detail {

inline constexpr std::size_t ambient_grid_points = 40000;

// 8-point Gauss-Legendre rule on [0, 1]: nodes g = (1 + xi) / 2 and weights w / 2 of the rule on
// [-1, 1] (the same values as refraction_horizon.hpp).
inline constexpr std::array<double, 8> gauss8_xi = {
    -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
    0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
inline constexpr std::array<double, 8> gauss8_w = {
    0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
    0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
inline constexpr std::array<double, 8> gauss8_nodes = {
    0.5 * (1.0 + gauss8_xi[0]), 0.5 * (1.0 + gauss8_xi[1]), 0.5 * (1.0 + gauss8_xi[2]),
    0.5 * (1.0 + gauss8_xi[3]), 0.5 * (1.0 + gauss8_xi[4]), 0.5 * (1.0 + gauss8_xi[5]),
    0.5 * (1.0 + gauss8_xi[6]), 0.5 * (1.0 + gauss8_xi[7])};
inline constexpr std::array<double, 8> gauss8_weights = {
    0.5 * gauss8_w[0], 0.5 * gauss8_w[1], 0.5 * gauss8_w[2], 0.5 * gauss8_w[3],
    0.5 * gauss8_w[4], 0.5 * gauss8_w[5], 0.5 * gauss8_w[6], 0.5 * gauss8_w[7]};

// Largest piece [m] of a straight path in the Gauss-Legendre quadrature of the air delay.
inline constexpr double straight_delay_piece_m = 250.0;

// Mean of f(s) over the path s in [0, length_m] (m): 8-point Gauss-Legendre on pieces of at most
// straight_delay_piece_m. f is called at the nodes only.
template <class F> double straight_path_mean(double length_m, F &&f) {
  if (!(length_m > 0.0)) {
    return f(0.0);
  }
  const auto pieces = static_cast<std::size_t>(std::ceil(length_m / straight_delay_piece_m));
  const double piece = length_m / static_cast<double>(pieces);
  double sum = 0.0;
  for (std::size_t p = 0; p < pieces; ++p) {
    const double start = piece * static_cast<double>(p);
    for (std::size_t i = 0; i < gauss8_nodes.size(); ++i) {
      sum += gauss8_weights[i] * piece * f(start + gauss8_nodes[i] * piece);
    }
  }
  return sum / length_m;
}

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

// The same from cos(zenith) = |cz| > 0 (at most 1). Use this where the zenith angle itself would
// round to pi/2 (|cz| below about 1e-16, acos returns pi/2 and the zenith check fails): the delay
// is then finite and exact for any |cz| > 0. The result is large for nearly horizontal paths.
inline double straight_path_air_delay_from_cos_ns(const AtmosphereProfile &profile, double low_m,
                                                  double high_m, double cos_zenith) {
  if (!std::isfinite(low_m) || !std::isfinite(high_m)) {
    throw std::invalid_argument("air delay needs finite heights");
  }
  if (high_m < low_m) {
    throw std::invalid_argument("air delay needs low <= high");
  }
  if (!std::isfinite(cos_zenith) || !(cos_zenith > 0.) || cos_zenith > 1.) {
    throw std::invalid_argument("air delay needs 0 < cos(zenith) <= 1");
  }
  const double integral_cm =
      profile.index_integral_cm(100. * low_m) - profile.index_integral_cm(100. * high_m);
  return integral_cm / 100. / (speed_of_light_m_per_ns * cos_zenith);
}

// Extra travel time of a straight path between two heights [m] compared with vacuum:
// (1/c) integral of (n-1) dh / cos(zenith), for the zenith angle zenith [rad], 0 <= zenith < pi/2.
inline double straight_path_air_delay_ns(const AtmosphereProfile &profile, double low_m,
                                         double high_m, double zenith) {
  if (!std::isfinite(zenith)) {
    throw std::invalid_argument("air delay needs finite heights and zenith angle");
  }
  if (zenith < 0. || zenith >= 1.5707963267948966) {
    throw std::invalid_argument("air delay needs 0 <= zenith < pi/2");
  }
  return straight_path_air_delay_from_cos_ns(profile, low_m, high_m, std::cos(zenith));
}

// Extra travel time [ns] of a straight path of length_m [m] whose height changes linearly from h_a
// to h_b [m a.s.l.] (flat geometry): (length / c) times the mean of k_group (n-1) along the path,
// by 8-point Gauss-Legendre quadrature on pieces of at most 250 m. The profile's n-1 is evaluated
// at each node; for h_a == h_b the result is exact, k_group (n-1) at that height times length / c.
// k_group = 1 gives the air delay of the unscaled profile. See docs/propagation/refraction.md.
inline double straight_line_air_delay_ns(const AtmosphereProfile &profile, double h_a_m,
                                         double h_b_m, double length_m, double k_group) {
  if (!std::isfinite(h_a_m) || !std::isfinite(h_b_m) || !std::isfinite(length_m) ||
      !std::isfinite(k_group) || length_m < 0.0 || k_group < 0.0) {
    throw std::invalid_argument(
        "straight line air delay needs finite heights and length, k_group >= 0");
  }
  if (!(length_m > 0.0)) {
    return 0.0;
  }
  const double slope = (h_b_m - h_a_m) / length_m; // height per metre of path
  const double mean = detail::straight_path_mean(
      length_m, [&](double s) { return profile.refractive_index_minus_one(h_a_m + slope * s); });
  return length_m * k_group * mean / speed_of_light_m_per_ns;
}

} // namespace obdeect::propagation
