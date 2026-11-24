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

#include "obdeect/propagation/dispersion.hpp"
#include "obdeect/propagation/parallel.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/timing.hpp"
#include "obdeect/stream/constants.hpp"

#include "obdeect/propagation/quadrature.hpp"
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
  // Derivative of excess_bent_ns with respect to the group scale at k_group = k = 1 [ns]: the
  // group index integral along the bent ray, (1/c) Int nu / cos theta_ref dx. Only set by lookups
  // and curved_ray_integrals; used by scale_ray_integrals.
  double bent_group_ns{};
};

// Integrals of the profile with n-1 scaled by k (phase) and the group index by k_group
// (dispersion), from the integrals of the unscaled profile, first order in the scale:
//   shift_m            = k shift_m
//   excess_bent_ns     = k excess_bent_ns + (k_group - k) bent_group_ns
//   excess_straight_ns = k_group excess_straight_ns
// The group term is the derivative of the bent excess along the bent ray (nu / cos theta_ref), the
// same as the accurate-dispersion tables. Exact at k = k_group = 1. See
// CurvedRefractionTables::lookup_scaled.
inline void scale_ray_integrals(CurvedRayIntegrals &v, double k, double k_group) noexcept {
  const double straight = v.excess_straight_ns;
  v.shift_m = k * v.shift_m;
  v.excess_bent_ns = k * v.excess_bent_ns + (k_group - k) * v.bent_group_ns;
  v.excess_straight_ns = k_group * straight;
}

// First node of the four-node stencil i1 - 1 .. i1 + 2 for the bracket starting at node i1 - 1
// (1 <= i1 <= n - 1), clamped to the grid. The validation against an independent ray tracer
// (obdeect-validation propagation.refraction-curved-ray) prefers this stencil to the one centred
// on the bracket (i1 - 2 .. i1 + 1), which doubles the time and plane-limit errors.
inline std::size_t stencil_start(std::size_t i1, std::size_t n) noexcept {
  return std::min(i1 - 1, n - 4);
}

namespace detail {

inline constexpr double curved_piece_m = 500.;      // largest sub-interval in r [m]
inline constexpr double curved_first_piece_m = 1.;  // first sub-interval near the telescope [m]
inline constexpr double curved_singular_cos = 0.2;  // substitution x = u^2 below this cos theta_s
inline constexpr double curved_angle_slack = 1e-12; // tolerance on theta_s in [0, pi/2]
inline constexpr double curved_cos_floor = 1e-4;    // table floor in cos theta_s (horizon grid)
inline constexpr double curved_cos_ratio = 1.02;    // geometric step in cos theta_s near horizon

// Integrands of the three integrals at offset x = r - r_t [m] from the telescope level.
struct CurvedIntegrands {
  double shift;       // [m]
  double bent;        // dimensionless, times dr
  double straight;    // dimensionless, times dr
  double group_delta; // dimensionless, d(bent)/d(k_group) times dr
};

// Common data of one ray family: telescope level, impact parameter, K = n_e b, cos theta_s and
// the phase and group refractivity scales.
struct CurvedRayFamily {
  double r_t;
  double h_t;
  double b;
  double cos_s;
  double K;
  double k;
  double k_group;
};

inline CurvedIntegrands curved_integrands(const AtmosphereProfile &profile,
                                          const CurvedRayFamily &family, double x) {
  const double r = family.r_t + x;
  const double nu = profile.index_minus_one_cm(100. * (family.h_t + x));
  const double n = 1. + family.k * nu;
  // cos theta_vac, in the stable form (r - r_t)(r + r_t) + r_t^2 cos^2 theta_s over r^2.
  const double cos_vac = std::sqrt(
      (x * (r + family.r_t) + family.r_t * family.r_t * family.cos_s * family.cos_s) / (r * r));
  const double s_vac = family.b / r;
  const double s_ref = family.K / (n * r);
  const double cos_ref = std::sqrt((1. - s_ref) * (1. + s_ref));
  const double n_group = 1. + family.k_group * nu;
  return {family.r_t * (s_ref / cos_ref - s_vac / cos_vac) / r, n_group / cos_ref - 1. / cos_vac,
          family.k_group * nu / cos_vac, nu / cos_ref};
}

} // namespace detail

// Integrals of the spherical ray from the telescope level (h_t) to the emission level (h_e) for
// the vacuum zenith angle theta_s at the telescope. Emission heights above the profile top are
// clamped to the top. Throws std::invalid_argument for non-finite input, R <= 0, theta_s outside
// [0, pi/2] (within 1e-12), a telescope altitude outside [bottom, top), and an emission below the
// telescope.
[[nodiscard]] inline CurvedRayIntegrals
curved_ray_integrals(const AtmosphereProfile &profile, double earth_radius_m,
                     double telescope_altitude_m, double theta_s_rad, double emission_altitude_m,
                     double k = 1., double k_group = 1., double *bent_group_delta_ns = nullptr) {
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
  if (!std::isfinite(k) || !(k > 0.) || !std::isfinite(k_group) || !(k_group > 0.)) {
    throw std::invalid_argument("phase and group refractivity scales must be finite and positive");
  }

  const double theta = std::clamp(theta_s_rad, 0., half_pi);
  const double h_e = std::min(emission_altitude_m, profile.top_m());
  const double r_t = earth_radius_m + h_t;
  const double sin_s = std::sin(theta);
  const double cos_s = std::cos(theta);
  const double b = r_t * sin_s;
  const double n_e = 1. + k * profile.index_minus_one_cm(100. * h_e);
  const detail::CurvedRayFamily family{r_t, h_t, b, cos_s, n_e * b, k, k_group};
  const double length = h_e - h_t; // [m], the integration range in x = r - r_t

  const bool near_telescope = cos_s < detail::curved_singular_cos;
  double shift = 0., bent = 0., straight = 0., group_delta = 0.;
  // x = r - r_t in [0, length]: u^2 substitution near the telescope (the integrands have 1/sqrt
  // singularities there), with the first piece curved_first_piece_m wide.
  detail::gauss_pieces(
      length, near_telescope ? detail::curved_first_piece_m : detail::curved_piece_m,
      detail::curved_piece_m, near_telescope, [&](double x, double weight) {
        const detail::CurvedIntegrands f = detail::curved_integrands(profile, family, x);
        shift += weight * f.shift;
        bent += weight * f.bent;
        straight += weight * f.straight;
        group_delta += weight * f.group_delta;
        return true;
      });

  const double c = speed_of_light_m_per_ns;
  if (bent_group_delta_ns != nullptr) {
    *bent_group_delta_ns = group_delta / c;
  }
  return {shift, bent / c, straight / c, group_delta / c};
}

// Grid of the integrals for one telescope altitude and Earth radius, built once with
// curved_ray_integrals at every node and interpolated cubically (4 x 4 Lagrange) in (theta_s, h_e).
// The build runs on `threads` threads (1: in place); the tables do not depend on the count.
class CurvedRefractionTables {
public:
  CurvedRefractionTables(AtmosphereProfile profile, double earth_radius_m,
                         double telescope_altitude_m, bool enable_accurate_dispersion = false,
                         unsigned threads = 1)
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
    group_.resize(count);
    // Every node is independent and writes its own slot: the thread count does not change the
    // result.
    const std::size_t nh = heights_.size();
    for_row_blocks(count, threads, [&](std::size_t begin, std::size_t end) {
      for (std::size_t k = begin; k < end; ++k) {
        const CurvedRayIntegrals v = curved_ray_integrals(
            profile_, earth_radius_m_, telescope_altitude_m_, theta_[k / nh], heights_[k % nh]);
        shift_[k] = v.shift_m;
        bent_[k] = v.excess_bent_ns;
        straight_[k] = v.excess_straight_ns;
        group_[k] = v.bent_group_ns;
      }
    });
    if (enable_accurate_dispersion) {
      build_accurate_dispersion_tables(threads);
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
    const std::size_t i_start = stencil_start(i1, theta_.size());
    const std::size_t j_start = stencil_start(j1, heights_.size());
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
    out.bent_group_ns = interpolate(group_, i_start, j_start, wi, wj);
    return true;
  }

  // Integrals for the profile with n-1 scaled by k (phase) and the group index by k_group. The
  // default table uses a first-order scale; tables constructed with accurate dispersion enabled
  // use wavelength-dependent integrals above 75 degrees, blending to the accurate result by 80:
  //   shift_m          = k shift_m(1)
  //   excess_bent_ns   = k excess_bent_ns(1) + (k_group - k) excess_straight_ns(1)
  //   excess_straight_ns = k_group excess_straight_ns(1)
  // (geometry with k, time with k_group: k G + k_group S = k (G + S) + (k_group - k) S). Exact at
  // k = k_group = 1. The error against direct quadrature is measured in test_dispersion.cpp.
  // False (out unchanged) for non-finite or non-positive k or k_group, or when lookup() fails.
  [[nodiscard]] bool lookup_scaled(double theta_s_rad, double emission_altitude_m, double k,
                                   double k_group, CurvedRayIntegrals &out) const noexcept {
    const double half_pi = obdeect::stream::pi / 2.;
    if (!std::isfinite(k) || !(k > 0.) || !std::isfinite(k_group) || !(k_group > 0.) ||
        !std::isfinite(theta_s_rad) || theta_s_rad < -detail::curved_angle_slack ||
        theta_s_rad > half_pi + detail::curved_angle_slack || !std::isfinite(emission_altitude_m) ||
        emission_altitude_m < telescope_altitude_m_) {
      return false;
    }
    const double theta = std::clamp(theta_s_rad, 0., half_pi);
    if (!dispersion_k_.empty() && theta >= accurate_dispersion_blend_start_rad()) {
      CurvedRayIntegrals accurate{};
      if (!lookup_accurate(theta, emission_altitude_m, k, k_group, accurate)) {
        return false;
      }
      if (theta >= accurate_dispersion_blend_end_rad()) {
        out = accurate;
        return true;
      }
      CurvedRayIntegrals base{};
      if (!lookup(theta, emission_altitude_m, base)) {
        return false;
      }
      scale_ray_integrals(base, k, k_group);
      const double u =
          (theta - accurate_dispersion_blend_start_rad()) /
          (accurate_dispersion_blend_end_rad() - accurate_dispersion_blend_start_rad());
      const double weight = u * u * (3. - 2. * u);
      out = {base.shift_m + weight * (accurate.shift_m - base.shift_m),
             base.excess_bent_ns + weight * (accurate.excess_bent_ns - base.excess_bent_ns),
             base.excess_straight_ns +
                 weight * (accurate.excess_straight_ns - base.excess_straight_ns)};
      return true;
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
  [[nodiscard]] static double accurate_dispersion_blend_start_rad() noexcept {
    return 75. * obdeect::stream::pi / 180.;
  }
  [[nodiscard]] static double accurate_dispersion_blend_end_rad() noexcept {
    return 80. * obdeect::stream::pi / 180.;
  }

private:
  static constexpr std::size_t accurate_k_node_count = 9;
  static constexpr std::size_t accurate_theta_start = 296; // 74 degrees; includes 75-degree stencil

  void build_accurate_dispersion_tables(unsigned threads) {
    const double low = refractivity_e8(2000.) / refractivity_e8(200.);
    const double high = 1. / low;
    dispersion_k_.resize(accurate_k_node_count);
    const std::size_t nt = theta_.size() - accurate_theta_start;
    const std::size_t nh = heights_.size();
    const std::size_t plane_size = nt * nh;
    dispersion_shift_.resize(accurate_k_node_count * plane_size);
    dispersion_bent_.resize(accurate_k_node_count * plane_size);
    dispersion_bent_delta_.resize(accurate_k_node_count * plane_size);
    dispersion_straight_.resize(accurate_k_node_count * plane_size);
    for (std::size_t q = 0; q < accurate_k_node_count; ++q) {
      dispersion_k_[q] = low + (high - low) * static_cast<double>(q) /
                                   static_cast<double>(accurate_k_node_count - 1);
    }
    // Flat index f = q * plane_size + (i - accurate_theta_start) * nh + j; each node is
    // independent.
    for_row_blocks(accurate_k_node_count * plane_size, threads,
                   [&](std::size_t begin, std::size_t end) {
                     for (std::size_t index = begin; index < end; ++index) {
                       const std::size_t q = index / plane_size;
                       const std::size_t rest = index % plane_size;
                       const std::size_t i = accurate_theta_start + rest / nh;
                       const std::size_t j = rest % nh;
                       const double k = dispersion_k_[q];
                       double bent_group_delta = 0.;
                       const CurvedRayIntegrals equal =
                           curved_ray_integrals(profile_, earth_radius_m_, telescope_altitude_m_,
                                                theta_[i], heights_[j], k, k, &bent_group_delta);
                       dispersion_shift_[index] = equal.shift_m;
                       dispersion_bent_[index] = equal.excess_bent_ns;
                       dispersion_bent_delta_[index] = bent_group_delta;
                       dispersion_straight_[index] = equal.excess_straight_ns / k;
                     }
                   });
  }

  [[nodiscard]] bool lookup_accurate(double theta, double emission_altitude_m, double k,
                                     double k_group, CurvedRayIntegrals &out) const noexcept {
    if (k < dispersion_k_.front() || k > dispersion_k_.back()) {
      return false;
    }
    const double h = std::min(emission_altitude_m, heights_.back());
    const double g = zenith_coordinate(theta);
    const std::size_t i1 =
        std::clamp(static_cast<std::size_t>(std::upper_bound(theta_.begin(), theta_.end(), theta) -
                                            theta_.begin()),
                   std::size_t{1}, theta_.size() - 1);
    const std::size_t j1 =
        std::clamp(static_cast<std::size_t>(std::upper_bound(heights_.begin(), heights_.end(), h) -
                                            heights_.begin()),
                   std::size_t{1}, heights_.size() - 1);
    const std::size_t i_start = stencil_start(i1, theta_.size());
    const std::size_t j_start = stencil_start(j1, heights_.size());
    if (i_start < accurate_theta_start) {
      return false;
    }
    std::array<double, 4> zenith_nodes{}, height_nodes{};
    for (std::size_t q = 0; q < 4; ++q) {
      zenith_nodes[q] = zenith_coordinate_[i_start + q];
      height_nodes[q] = heights_[j_start + q];
    }
    const auto wi = lagrange_weights(zenith_nodes, g);
    const auto wj = lagrange_weights(height_nodes, h);
    const auto ki = static_cast<std::size_t>(
        std::upper_bound(dispersion_k_.begin(), dispersion_k_.end(), k) - dispersion_k_.begin());
    const std::size_t k0 =
        std::clamp(ki == 0 ? 0 : ki - 1, std::size_t{0}, dispersion_k_.size() - 4);
    std::array<double, 4> kn{};
    for (std::size_t q = 0; q < 4; ++q) {
      kn[q] = dispersion_k_[k0 + q];
    }
    const auto wk = lagrange_weights(kn, k);
    const std::size_t plane_size = (theta_.size() - accurate_theta_start) * heights_.size();
    auto interp = [&](const std::vector<double> &values) {
      double value = 0.;
      for (std::size_t q = 0; q < 4; ++q) {
        const std::size_t base = (k0 + q) * plane_size;
        double spatial = 0.;
        for (std::size_t a = 0; a < 4; ++a) {
          double row = 0.;
          for (std::size_t b = 0; b < 4; ++b) {
            const std::size_t index =
                base + (i_start + a - accurate_theta_start) * heights_.size() + j_start + b;
            row += wj[b] * values[index];
          }
          spatial += wi[a] * row;
        }
        value += wk[q] * spatial;
      }
      return value;
    };
    const double shift = interp(dispersion_shift_);
    const double bent_equal = interp(dispersion_bent_);
    const double bent_delta = interp(dispersion_bent_delta_);
    const double straight_unit = interp(dispersion_straight_);
    out = {shift, bent_equal + (k_group - k) * bent_delta, k_group * straight_unit};
    return true;
  }

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
  std::vector<double> shift_, bent_, straight_, group_;
  std::vector<double> dispersion_k_;
  std::vector<double> dispersion_shift_, dispersion_bent_, dispersion_bent_delta_,
      dispersion_straight_;
};

} // namespace obdeect::propagation
