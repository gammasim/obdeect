#pragma once

// Ray bending and air delay in a spherical atmosphere: Earth radius R, refractive index
// n(r) = 1 + nu(h) with r = R + h and h the altitude a.s.l. of the profile.
//
// A photon is emitted at altitude h_e. Its vacuum straight line crosses the telescope level
// r_t = R + h_t at zenith angle theta_s (sin theta_s = b / r_t, b the impact parameter about
// the Earth centre). The real ray starts with the same direction at the emission point and
// follows Bouguer's law n(r) r sin(theta) = K = n_e b, theta measured from the local radial
// direction. With x = r - r_t the integrals from the telescope level to the emission level are
//
//   shift_m          = r_t Int (tan theta_ref - tan theta_vac) / r dx
//   excess_bent_ns   = (1/c) Int (n / cos theta_ref - 1 / cos theta_vac) dx
//   excess_straight_ns = (1/c) Int (nu / cos theta_vac) dx
//
// with s_ref = K / (n r), cos theta_ref = sqrt((1 - s_ref)(1 + s_ref)), sin theta_vac = b / r and
// cos theta_vac = sqrt(((r - r_t)(r + r_t) + r_t^2 cos^2 theta_s) / r^2). The shift is the
// lateral displacement of the arrival point along the horizontal direction of motion (negative:
// towards the point below the emission); the arrival direction needs no table (sin theta' /
// sin theta_s = n_e / n_t). In the flat limit R -> infinity these reduce to bend() and
// straight_path_air_delay_ns() of refraction.hpp and timing.hpp.
//
// The integrals are composite 8-point Gauss-Legendre sums on sub-intervals of at most 500 m of r.
// For cos theta_s < 0.2 the telescope-level singularity (cos theta_vac ~ sqrt(x) for theta_s = 90
// degrees) is removed by x = u^2 and geometrically refined first sub-intervals in x (1, 2, 4, ...
// m). Everything is double; no result depends on how callers batch photons.
//
// CurvedRefractionTables holds the integrals on a (theta_s, h_e) grid for one telescope altitude
// and interpolates them cubically (4 x 4 Lagrange). Heights above the profile top (or above 120 km,
// whichever is lower) take the values at that height.

#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/timing.hpp"
#include "obdeect/stream/constants.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace obdeect::propagation {

struct CurvedRayIntegrals {
  double shift_m{};
  double excess_bent_ns{};
  double excess_straight_ns{};
};

// Integrals of the profile with n-1 scaled by k (phase) and the group index by k_group
// (dispersion), from the integrals of the unscaled profile, first order in the scale:
//   shift_m            = k shift_m
//   excess_bent_ns     = k excess_bent_ns + (k_group - k) excess_straight_ns
//   excess_straight_ns = k_group excess_straight_ns
// Exact at k = k_group = 1. See CurvedRefractionTables::lookup_scaled.
inline void scale_ray_integrals(CurvedRayIntegrals &v, double k, double k_group) noexcept {
  const double straight = v.excess_straight_ns;
  v.shift_m = k * v.shift_m;
  v.excess_bent_ns = k * v.excess_bent_ns + (k_group - k) * straight;
  v.excess_straight_ns = k_group * straight;
}

namespace detail {

// 8-point Gauss-Legendre nodes and weights on [-1, 1].
inline constexpr std::array<double, 8> curved_gauss_nodes = {
    -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
    0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
inline constexpr std::array<double, 8> curved_gauss_weights = {
    0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
    0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};

inline constexpr double curved_piece_m = 500.;      // largest sub-interval in r [m]
inline constexpr double curved_first_piece_m = 1.;  // first sub-interval near the telescope [m]
inline constexpr double curved_singular_cos = 0.2;  // substitution x = u^2 below this cos theta_s
inline constexpr double curved_angle_slack = 1e-12; // tolerance on theta_s in [0, pi/2]
inline constexpr double curved_cos_floor = 1e-4;    // table floor in cos theta_s (horizon grid)
inline constexpr double curved_cos_ratio = 1.02;    // geometric step in cos theta_s near horizon

// Integrands of the three integrals at offset x = r - r_t [m] from the telescope level.
struct CurvedIntegrands {
  double shift;    // [m]
  double bent;     // dimensionless, times dr
  double straight; // dimensionless, times dr
};

// Common data of one ray family: telescope level, impact parameter, K = n_e b, cos theta_s.
struct CurvedRayFamily {
  double r_t;
  double h_t;
  double b;
  double cos_s;
  double K;
};

inline CurvedIntegrands curved_integrands(const AtmosphereProfile &profile,
                                          const CurvedRayFamily &family, double x) {
  const double r = family.r_t + x;
  const double nu = profile.index_minus_one_cm(100. * (family.h_t + x));
  const double n = 1. + nu;
  // cos theta_vac, in the stable form (r - r_t)(r + r_t) + r_t^2 cos^2 theta_s over r^2.
  const double cos_vac = std::sqrt(
      (x * (r + family.r_t) + family.r_t * family.r_t * family.cos_s * family.cos_s) / (r * r));
  const double s_vac = family.b / r;
  const double s_ref = family.K / (n * r);
  const double cos_ref = std::sqrt((1. - s_ref) * (1. + s_ref));
  return {family.r_t * (s_ref / cos_ref - s_vac / cos_vac) / r, n / cos_ref - 1. / cos_vac,
          nu / cos_vac};
}

} // namespace detail

// Integrals of the spherical ray from the telescope level (h_t) to the emission level (h_e) for
// the vacuum zenith angle theta_s at the telescope. Emission heights above the profile top are
// clamped to the top. Throws std::invalid_argument for non-finite input, R <= 0, theta_s outside
// [0, pi/2] (within 1e-12), a telescope altitude outside [bottom, top), and an emission below the
// telescope.
[[nodiscard]] inline CurvedRayIntegrals
curved_ray_integrals(const AtmosphereProfile &profile, double earth_radius_m,
                     double telescope_altitude_m, double theta_s_rad, double emission_altitude_m) {
  if (!std::isfinite(earth_radius_m) || !(earth_radius_m > 0.)) {
    throw std::invalid_argument("earth radius must be finite and positive");
  }
  const double h_t = telescope_altitude_m;
  if (!std::isfinite(h_t) || h_t < profile.bottom_m() || h_t >= profile.top_m()) {
    throw std::invalid_argument("telescope altitude must lie inside the atmospheric profile");
  }
  const double half_pi = obdeect::stream::pi / 2.;
  if (!std::isfinite(theta_s_rad) || theta_s_rad < -detail::curved_angle_slack ||
      theta_s_rad > half_pi + detail::curved_angle_slack) {
    throw std::invalid_argument("zenith angle at the telescope must lie in [0, pi/2]");
  }
  if (!std::isfinite(emission_altitude_m)) {
    throw std::invalid_argument("emission altitude must be finite");
  }
  if (emission_altitude_m < h_t) {
    throw std::invalid_argument("emission altitude must not lie below the telescope");
  }

  const double theta = std::clamp(theta_s_rad, 0., half_pi);
  const double h_e = std::min(emission_altitude_m, profile.top_m());
  const double r_t = earth_radius_m + h_t;
  const double sin_s = std::sin(theta);
  const double cos_s = std::cos(theta);
  const double b = r_t * sin_s;
  const double n_e = 1. + profile.index_minus_one_cm(100. * h_e);
  const detail::CurvedRayFamily family{r_t, h_t, b, cos_s, n_e * b};
  const double length = h_e - h_t; // [m], the integration range in x = r - r_t

  const bool near_telescope = cos_s < detail::curved_singular_cos;
  double shift = 0., bent = 0., straight = 0.;
  auto accumulate = [&](double weight, double x) {
    const detail::CurvedIntegrands f = detail::curved_integrands(profile, family, x);
    shift += weight * f.shift;
    bent += weight * f.bent;
    straight += weight * f.straight;
  };

  double x0 = 0.;
  double width = near_telescope ? detail::curved_first_piece_m : detail::curved_piece_m;
  while (x0 < length) {
    const double x1 = std::min(length, x0 + std::min(width, detail::curved_piece_m));
    if (near_telescope) {
      // x = u^2: Gauss on [sqrt(x0), sqrt(x1)] with the Jacobian 2u.
      const double u0 = std::sqrt(x0), u1 = std::sqrt(x1);
      const double mid = 0.5 * (u0 + u1), half = 0.5 * (u1 - u0);
      for (std::size_t k = 0; k < detail::curved_gauss_nodes.size(); ++k) {
        const double u = mid + half * detail::curved_gauss_nodes[k];
        accumulate(half * 2. * u * detail::curved_gauss_weights[k], u * u);
      }
    } else {
      const double mid = 0.5 * (x0 + x1), half = 0.5 * (x1 - x0);
      for (std::size_t k = 0; k < detail::curved_gauss_nodes.size(); ++k) {
        accumulate(half * detail::curved_gauss_weights[k],
                   mid + half * detail::curved_gauss_nodes[k]);
      }
    }
    x0 = x1;
    width = std::min(2. * width, detail::curved_piece_m);
  }

  const double c = speed_of_light_m_per_ns;
  return {shift, bent / c, straight / c};
}

// Grid of the integrals for one telescope altitude and Earth radius, built once with
// curved_ray_integrals at every node and interpolated cubically (4 x 4 Lagrange) in (theta_s, h_e).
class CurvedRefractionTables {
public:
  CurvedRefractionTables(AtmosphereProfile profile, double earth_radius_m,
                         double telescope_altitude_m)
      : profile_(std::move(profile)), earth_radius_m_(earth_radius_m),
        telescope_altitude_m_(telescope_altitude_m) {
    // Validates the arguments (the same checks as curved_ray_integrals at the telescope level).
    static_cast<void>(curved_ray_integrals(profile_, earth_radius_m_, telescope_altitude_m_, 0.,
                                           telescope_altitude_m_));
    build_theta_nodes();
    build_height_nodes();
    const std::size_t count = theta_.size() * heights_.size();
    shift_.resize(count);
    bent_.resize(count);
    straight_.resize(count);
    for (std::size_t i = 0; i < theta_.size(); ++i) {
      for (std::size_t j = 0; j < heights_.size(); ++j) {
        const CurvedRayIntegrals v = curved_ray_integrals(
            profile_, earth_radius_m_, telescope_altitude_m_, theta_[i], heights_[j]);
        const std::size_t k = i * heights_.size() + j;
        shift_[k] = v.shift_m;
        bent_[k] = v.excess_bent_ns;
        straight_[k] = v.excess_straight_ns;
      }
    }
  }

  // Cubic (4 x 4 Lagrange) interpolation on the nodes around the query point. False (out
  // unchanged) if an input is not finite, theta_s lies outside [0, pi/2] (within 1e-12) or h_e is
  // below the telescope. Heights above the grid take the value at its top. The stencil stays
  // inside the grid (no extrapolation). No allocation, no exceptions.
  [[nodiscard]] bool lookup(double theta_s_rad, double emission_altitude_m,
                            CurvedRayIntegrals &out) const noexcept {
    const double half_pi = obdeect::stream::pi / 2.;
    if (!std::isfinite(theta_s_rad) || !std::isfinite(emission_altitude_m) ||
        theta_s_rad < -detail::curved_angle_slack ||
        theta_s_rad > half_pi + detail::curved_angle_slack ||
        emission_altitude_m < telescope_altitude_m_) {
      return false;
    }
    const double theta = std::clamp(theta_s_rad, 0., half_pi);
    const double h = std::min(emission_altitude_m, heights_.back());
    const double g = zenith_coordinate(theta);

    // Stencil: four consecutive nodes around the bracketing interval [i0, i0 + 1].
    const std::size_t i1 =
        std::clamp(static_cast<std::size_t>(std::upper_bound(theta_.begin(), theta_.end(), theta) -
                                            theta_.begin()),
                   std::size_t{1}, theta_.size() - 1);
    const std::size_t j1 =
        std::clamp(static_cast<std::size_t>(std::upper_bound(heights_.begin(), heights_.end(), h) -
                                            heights_.begin()),
                   std::size_t{1}, heights_.size() - 1);
    const std::size_t i_start = std::min(i1 - 1, theta_.size() - 4);
    const std::size_t j_start = std::min(j1 - 1, heights_.size() - 4);
    std::array<double, 4> zenith_nodes{}, height_nodes{};
    for (std::size_t k = 0; k < 4; ++k) {
      zenith_nodes[k] = zenith_coordinate_[i_start + k];
      height_nodes[k] = heights_[j_start + k];
    }
    const std::array<double, 4> wi = lagrange_weights(zenith_nodes, g);
    const std::array<double, 4> wj = lagrange_weights(height_nodes, h);

    out.shift_m = interpolate(shift_, i_start, j_start, wi, wj);
    out.excess_bent_ns = interpolate(bent_, i_start, j_start, wi, wj);
    out.excess_straight_ns = interpolate(straight_, i_start, j_start, wi, wj);
    return true;
  }

  // Integrals for the profile with n-1 scaled by k (phase) and the group index by k_group, from the
  // tables of the unscaled profile, without rebuilding them. First order in the scale:
  //   shift_m          = k shift_m(1)
  //   excess_bent_ns   = k excess_bent_ns(1) + (k_group - k) excess_straight_ns(1)
  //   excess_straight_ns = k_group excess_straight_ns(1)
  // (geometry with k, time with k_group: k G + k_group S = k (G + S) + (k_group - k) S). Exact at
  // k = k_group = 1. The error against direct quadrature is measured in test_dispersion.cpp.
  // False (out unchanged) for non-finite or non-positive k or k_group, or when lookup() fails.
  [[nodiscard]] bool lookup_scaled(double theta_s_rad, double emission_altitude_m, double k,
                                   double k_group, CurvedRayIntegrals &out) const noexcept {
    if (!std::isfinite(k) || !(k > 0.) || !std::isfinite(k_group) || !(k_group > 0.)) {
      return false;
    }
    CurvedRayIntegrals base{};
    if (!lookup(theta_s_rad, emission_altitude_m, base)) {
      return false;
    }
    scale_ray_integrals(base, k, k_group);
    out = base;
    return true;
  }

  [[nodiscard]] double earth_radius_m() const noexcept { return earth_radius_m_; }
  [[nodiscard]] double telescope_altitude_m() const noexcept { return telescope_altitude_m_; }
  [[nodiscard]] const AtmosphereProfile &profile() const noexcept { return profile_; }
  // Grid sizes: zenith-angle nodes and emission-height nodes.
  [[nodiscard]] std::size_t theta_nodes() const noexcept { return theta_.size(); }
  [[nodiscard]] std::size_t height_nodes() const noexcept { return heights_.size(); }

private:
  // Lagrange weights of the four nodes p at x: exact (1 at its node, 0 at the others) at the nodes.
  static std::array<double, 4> lagrange_weights(const std::array<double, 4> &p, double x) noexcept {
    std::array<double, 4> w{};
    for (std::size_t k = 0; k < 4; ++k) {
      double value = 1.;
      for (std::size_t m = 0; m < 4; ++m) {
        if (m != k) {
          value *= (x - p[m]) / (p[k] - p[m]);
        }
      }
      w[k] = value;
    }
    return w;
  }

  // Cubic value from the 4 x 4 nodes starting at (i_start, j_start) with weights wi (zenith) and wj
  // (height).
  [[nodiscard]] double interpolate(const std::vector<double> &values, std::size_t i_start,
                                   std::size_t j_start, const std::array<double, 4> &wi,
                                   const std::array<double, 4> &wj) const noexcept {
    const std::size_t nh = heights_.size();
    double value = 0.;
    for (std::size_t k = 0; k < 4; ++k) {
      double row = 0.;
      for (std::size_t l = 0; l < 4; ++l) {
        row += wj[l] * values[(i_start + k) * nh + j_start + l];
      }
      value += wi[k] * row;
    }
    return value;
  }

  // Zenith-angle coordinate in which the table is interpolated: ln(cos theta + floor). The
  // integrals behave like cos(theta)^-2.75 near the horizon, so this coordinate is almost linear.
  static double zenith_coordinate(double theta) {
    return std::log(std::cos(theta) + detail::curved_cos_floor);
  }

  // 0.25 deg up to 80 deg, then cos theta decreasing geometrically by 2 % per node down to
  // detail::curved_cos_floor, then 90 deg.
  void build_theta_nodes() {
    const double deg = obdeect::stream::pi / 180.;
    for (int k = 0; k <= 320; ++k) {
      theta_.push_back(0.25 * k * deg);
    }
    double cos_node = std::cos(80. * deg);
    while (cos_node / detail::curved_cos_ratio > detail::curved_cos_floor) {
      cos_node /= detail::curved_cos_ratio;
      theta_.push_back(std::acos(cos_node));
    }
    theta_.push_back(obdeect::stream::pi / 2.);
    zenith_coordinate_.resize(theta_.size());
    for (std::size_t i = 0; i < theta_.size(); ++i) {
      zenith_coordinate_[i] = zenith_coordinate(theta_[i]);
    }
  }

  // Geometric steps of 5 % from 0.1 m to 20 km above the telescope (the integrals are quadratic in
  // the height near the telescope level for grazing rays), then 1 km steps up to
  // min(profile top, 120 km); the top node is the upper limit itself.
  void build_height_nodes() {
    const double h_t = telescope_altitude_m_;
    const double h_max = std::min(profile_.top_m(), std::max(120000., h_t + 1000.));
    heights_.push_back(h_t);
    for (double d = 0.1; d < 20000.; d *= 1.05) {
      if (h_t + d >= h_max) {
        break;
      }
      heights_.push_back(h_t + d);
    }
    for (int k = 1;; ++k) {
      const double h = h_t + 20000. + 1000. * k;
      if (h >= h_max) {
        break;
      }
      heights_.push_back(h);
    }
    heights_.push_back(h_max);
    // The 4-node stencil of lookup() needs at least four heights; a telescope within a fraction of
    // a metre of the profile top gets a short uniform grid instead.
    constexpr std::size_t min_heights = 8;
    if (heights_.size() < min_heights) {
      heights_.clear();
      for (std::size_t k = 0; k < min_heights; ++k) {
        heights_.push_back(h_t + (h_max - h_t) * static_cast<double>(k) /
                                     static_cast<double>(min_heights - 1));
      }
    }
  }

  AtmosphereProfile profile_;
  double earth_radius_m_{};
  double telescope_altitude_m_{};
  std::vector<double> theta_;             // zenith-angle nodes [rad], ascending
  std::vector<double> zenith_coordinate_; // ln(cos theta + floor) at the nodes
  std::vector<double> heights_;           // emission altitudes [m a.s.l.], ascending
  // Row-major [theta][height] arrays.
  std::vector<double> shift_, bent_, straight_;
};

} // namespace obdeect::propagation
