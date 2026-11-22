#pragma once

// Exact plane-parallel bending near the horizon. The series of refraction.hpp (bend) is only
// valid up to about 87 degrees zenith; this header evaluates the same quantities directly from
// Snell's law, n(h) sin(theta(h)) = n_e sin(theta_e), with Gauss-Legendre quadrature. See
// docs/propagation/refraction.md ("Near the horizon").
//
// With h the height, nu = n - 1, s = sin(theta), and the path from the emission level h_e down
// to the observation level:
//   dx     = integral of (tan theta(h) - tan theta_e) dh            (displacement, negative)
//   T_air  = (1/c) integral of (n_group(h)/cos theta(h) - 1/cos theta_e) dh (air excess)
// The travel time is (h_e - obs) / (c cos theta_e) + T_air, as in bend().
//
// Near the horizon (cos theta_e < 0.2) the integrals are split at 1 m, 2 m, 4 m, ... and the
// sub-intervals are at most 250 m. The first sub-intervals use x = u^2 (Jacobian 2u) to resolve
// the square-root behaviour of 1/cos theta at the emission level.
//
// Public entry points use the same conventions as bend() (see refraction.hpp) and the same
// rejection rules, including the dispersion factors k and k_group (see bend()). The series limit
// bend_series_limit_exceeded() is defined in refraction.hpp; a photon that bend() rejects for that
// reason is bent by bend_quadrature() instead of being dropped.

#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/timing.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace obdeect::propagation {

namespace detail {

// 8-point Gauss-Legendre nodes and weights on [-1, 1].
inline constexpr std::array<double, 8> horizon_gauss_nodes = {
    -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
    0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
inline constexpr std::array<double, 8> horizon_gauss_weights = {
    0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
    0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};

inline constexpr double horizon_piece_m = 250.;     // default largest sub-interval [m]
inline constexpr double horizon_first_piece_m = 1.; // first sub-interval at the emission [m]
inline constexpr double horizon_singular_cos = 0.2; // x = u^2 substitution below this cos theta_e

// Integrands at height h [cm] for a ray with sin theta_e = s_e, cos theta_e = c_e, nu_e = k (n_e -
// 1) (k: dispersion factor of the profile, 1 without dispersion):
//   f_tan  = tan theta(h) - tan theta_e,
//   f_time = n_group(h) / cos theta(h) - 1 / cos theta_e.
// Returns false if sin theta(h) >= 1 (the ray cannot reach h).
inline bool horizon_integrands(const AtmosphereProfile &profile, double phase_scale, double h_cm,
                               double group_scale, double s_e, double c_e, double nu_e,
                               double &f_tan, double &f_time) {
  const double profile_nu = profile.index_minus_one_cm(h_cm);
  const double nu = phase_scale * profile_nu;
  const double s = s_e * (1. + nu_e) / (1. + nu); // Snell: n sin(theta) = n_e sin(theta_e)
  const double c2 = 1. - s * s;
  if (!(s < 1.) || !(c2 > 0.)) {
    return false;
  }
  const double c = std::sqrt(c2);
  f_tan = s / c - s_e / c_e;
  f_time = (1. + group_scale * profile_nu) / c - 1. / c_e;
  return true;
}

// Integrals over a path of length_cm below the emission level h_e_cm [cm]. Outputs: dx_cm and the
// air excess [cm] (time integral times c, without the vacuum part). Returns false if the ray
// cannot reach the lower end of the path.
inline bool horizon_path_integrals(const AtmosphereProfile &profile, double phase_scale,
                                   double group_scale, double h_e_cm, double length_cm, double s_e,
                                   double c_e, double nu_e, double piece_cm, double &dx_cm,
                                   double &excess_cm) {
  dx_cm = 0.;
  excess_cm = 0.;
  const bool near_horizontal = c_e < horizon_singular_cos;
  const double first_cm = std::min(100. * horizon_first_piece_m, piece_cm);
  double p0 = 0.; // distance below the emission level [cm]
  double width = first_cm;
  while (p0 < length_cm) {
    const double p1 = std::min(length_cm, p0 + width);
    for (std::size_t k = 0; k < horizon_gauss_nodes.size(); ++k) {
      double p = 0., weight = 0.;
      if (near_horizontal) {
        // p = u^2 on [sqrt(p0), sqrt(p1)], dp = 2 u du.
        const double u0 = std::sqrt(p0), u1 = std::sqrt(p1);
        const double u = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * horizon_gauss_nodes[k];
        p = u * u;
        weight = 0.5 * (u1 - u0) * 2. * u * horizon_gauss_weights[k];
      } else {
        p = 0.5 * (p0 + p1) + 0.5 * (p1 - p0) * horizon_gauss_nodes[k];
        weight = 0.5 * (p1 - p0) * horizon_gauss_weights[k];
      }
      double f_tan = 0., f_time = 0.;
      if (!horizon_integrands(profile, phase_scale, h_e_cm - p, group_scale, s_e, c_e, nu_e, f_tan,
                              f_time)) {
        return false;
      }
      dx_cm += weight * f_tan;
      excess_cm += weight * f_time;
    }
    p0 = p1;
    width = std::min(2. * width, piece_cm);
  }
  return true;
}

// Shared body of bend_quadrature and bend_quadrature_pieces. piece_m is the largest sub-interval.
// k and k_group as in bend(): the geometry uses k nu, the time uses k_group nu.
template <class T>
[[nodiscard]] bool bend_quadrature_impl(const RefractionTables &tables, double piece_m, double k,
                                        double k_group, T height_m, T &cx, T &cy, T &x, T &y,
                                        T &t) {
  const double h_m = static_cast<double>(height_m);
  const double ux = static_cast<double>(cx), uy = static_cast<double>(cy);
  if (!std::isfinite(h_m) || !std::isfinite(ux) || !std::isfinite(uy) ||
      !std::isfinite(static_cast<double>(x)) || !std::isfinite(static_cast<double>(y)) ||
      !std::isfinite(static_cast<double>(t)) || !std::isfinite(k) || !(k > 0.) ||
      !std::isfinite(k_group) || !(k_group > 0.)) {
    return false;
  }
  const double h = 100. * h_m;
  const double s2 = ux * ux + uy * uy;
  const double obs = tables.observation_level_cm();
  if (!std::isfinite(h) || !(s2 < 1.) || h < obs) {
    return false;
  }
  const double s_e = std::sqrt(s2);
  const double c_e = std::sqrt(1. - s2);
  const double nu_e = k * tables.profile().index_minus_one_cm(h);
  const double nu_obs = k * tables.index_minus_one_at_observation();
  const double s_obs = s_e * (1. + nu_e) / (1. + nu_obs);
  if (!(s_obs < 1.)) {
    return false; // total internal reflection
  }

  double dx_cm = 0., excess_cm = 0.;
  if (!horizon_path_integrals(tables.profile(), k, k_group, h, h - obs, s_e, c_e, nu_e,
                              100. * piece_m, dx_cm, excess_cm)) {
    return false;
  }

  // Travel time: vacuum path plus the group-index excess along the bent path.
  const double travel_cm = (h - obs) / c_e + excess_cm;
  const double travel_ns = travel_cm / (100. * speed_of_light_m_per_ns);
  if (s_e > 0.) {
    const double shift_m = dx_cm / 100.;
    x = static_cast<T>(static_cast<double>(x) + shift_m * ux / s_e);
    y = static_cast<T>(static_cast<double>(y) + shift_m * uy / s_e);
    const double scale = s_obs / s_e;
    cx = static_cast<T>(ux * scale);
    cy = static_cast<T>(uy * scale);
  }
  t = static_cast<T>(static_cast<double>(t) + travel_ns);
  return true;
}

} // namespace detail

// Same as bend_quadrature() with an explicit largest sub-interval piece_m [m] (> 0). Used by the
// tests to check convergence; production code calls bend_quadrature().
template <class T>
[[nodiscard]] bool bend_quadrature_pieces(const RefractionTables &tables, double piece_m,
                                          T height_m, T &cx, T &cy, T &x, T &y, T &t) {
  if (!std::isfinite(piece_m) || !(piece_m > 0.)) {
    throw std::invalid_argument("quadrature piece size must be finite and positive");
  }
  return detail::bend_quadrature_impl(tables, piece_m, 1., 1., height_m, cx, cy, x, y, t);
}

// Exact plane-parallel bending by quadrature, for photons that bend() rejects at the series limit.
// Same arguments, outputs and rejection rules as bend() in refraction.hpp, with the dispersion
// factors k and k_group of bend(). Returns false and leaves the photon unchanged if an input is not
// finite, the direction is not strictly inside the unit circle, the emission is below the
// observation level, or the ray cannot reach the observation level (sin theta >= 1). Allocates
// nothing; about 8 profile lookups per 250 m of path.
template <class T>
[[nodiscard]] bool bend_quadrature(const RefractionTables &tables, T height_m, T &cx, T &cy, T &x,
                                   T &y, T &t, double k, double k_group) {
  return detail::bend_quadrature_impl(tables, detail::horizon_piece_m, k, k_group, height_m, cx, cy,
                                      x, y, t);
}

// bend_quadrature() without dispersion: k = k_group = 1.
template <class T>
[[nodiscard]] bool bend_quadrature(const RefractionTables &tables, T height_m, T &cx, T &cy, T &x,
                                   T &y, T &t) {
  return bend_quadrature(tables, height_m, cx, cy, x, y, t, 1., 1.);
}

} // namespace obdeect::propagation
