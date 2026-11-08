#pragma once

#include "dlufd/profile.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace dlufd {

// Refraction tables for one profile and observation level: air-time excess and ray-bending
// coefficients from planar ray integration (see docs/refraction.md).
// Units inside: cm and ns. Public heights are m a.s.l.
class RefractionTables {
public:
  static constexpr double vacuum_speed_cm_per_ns = 29.9792458; // speed of light in vacuum
  static constexpr double pi = 3.14159265358979323846;

  RefractionTables(AtmosphereProfile profile, double observation_level_m)
      : profile_(std::move(profile)), obs_cm_(100. * observation_level_m) {
    require_finite(obs_cm_, "observation level");
    if (!(obs_cm_ >= profile_.bottom_cm() && obs_cm_ < profile_.top_cm()))
      throw std::invalid_argument("observation level must lie inside the atmospheric table");
    const auto &alt = profile_.alt_cm_;
    const auto &rho = profile_.rho_;
    const auto &n_1 = profile_.n_1_;
    const std::size_t n = alt.size();

    // etadsn: n-1 per density at the first level more than 1.5 km above the observation level.
    etadsn_ = AtmosphereProfile::default_etadsn;
    for (std::size_t i = 0; i < n; ++i) {
      if (alt[i] > obs_cm_ + 1.5e5) {
        etadsn_ = n_1[i] / rho[i];
        break;
      }
    }
    obs_refidx_ = 1. + profile_.index_minus_one_cm(obs_cm_);
    obs_thick_ = profile_.thickness_cm(obs_cm_);

    // Bending coefficients per table level, from a 45 degree ray integrated numerically.
    const double theta0 = 45. * (pi / 180.);
    std::vector<double> hori_a(n), time_a(n);
    table_time0_.resize(n);
    table_alt_ = alt;
    for (std::size_t ialt = 0; ialt < n; ++ialt) {
      const double n_sint0 = (1. + profile_.index_minus_one_cm(alt[ialt])) * std::sin(theta0);
      double x0 = 0., t_vt = 0., s0 = 0.;
      double x_45 = 0., t_45 = 0., s45 = 0.;
      const double theta1 = std::asin(n_sint0 / obs_refidx_);
      const double c = std::cos(theta0 + 0.28 * (theta1 - theta0));
      const double s = std::sin(theta0 + 0.28 * (theta1 - theta0));
      const double vci = 1. / vacuum_speed_cm_per_ns;
      const double airmass = 1. / std::cos(theta0);
      const double zem = alt[ialt];
      const double dz = zem - obs_cm_;
      const double de = etadsn_ * (obs_thick_ - profile_.thickness_cm(zem));
      const double t0_vt = (dz + de) * vci;
      const double t0_45 = t0_vt * airmass;
      trace_ray_planar(zem, obs_cm_, 0., 2.5e2, x0, t_vt, s0);
      trace_ray_planar(zem, obs_cm_, theta0, 2.5e2, x_45, t_45, s45);
      if (zem < obs_cm_) {
        t_vt *= -1.;
        t_45 *= -1.;
        x_45 *= -1.;
      }
      table_time0_[ialt] = t_vt - t0_vt;
      // Density-indexed coefficients are stored in order of increasing density.
      const std::size_t rev = n - 1 - ialt;
      hori_a[rev] =
          x_45 >= 0. ? std::sqrt(x_45 * (c * c * c) / s) : -std::sqrt(-x_45 * (c * c * c) / s);
      if (((t_45 - t0_45) - (t_vt - t0_vt)) < 0.)
        time_a[rev] = std::sqrt(((t0_45 - t_45) - (t0_vt - t_vt)) * (c * c * c) / (s * s));
      else
        time_a[rev] = std::sqrt(-((t0_45 - t_45) - (t0_vt - t_vt)) * (c * c * c) / (s * s));
    }

    // Cubic splines of the coefficients, resampled on the equidistant tables.
    std::vector<double> rho_rev(n);
    for (std::size_t i = 0; i < n; ++i)
      rho_rev[i] = rho[n - 1 - i];
    const CubicSpline spline_hori(rho_rev, hori_a);
    const CubicSpline spline_time_a(rho_rev, time_a);
    const CubicSpline spline_time0(alt, table_time0_);
    fast_rho_fac_ = static_cast<double>(AtmosphereProfile::fast_points - 1) / rho.front();
    fast_rho_rev_.resize(AtmosphereProfile::fast_points);
    fast_hori_a_.resize(AtmosphereProfile::fast_points);
    fast_time_a_.resize(AtmosphereProfile::fast_points);
    fast_time0_.resize(AtmosphereProfile::fast_points);
    const auto fast_alt = profile_.fast_altitude_cm();
    for (std::size_t i = 0; i < AtmosphereProfile::fast_points; ++i) {
      const double r = i < AtmosphereProfile::fast_points - 1
                           ? static_cast<double>(i) / fast_rho_fac_
                           : rho.front();
      fast_rho_rev_[i] = r;
      fast_hori_a_[i] = spline_hori.value(r);
      fast_time_a_[i] = spline_time_a.value(r);
      fast_time0_[i] = spline_time0.value(fast_alt[i]);
    }
  }

  const AtmosphereProfile &profile() const { return profile_; }
  double observation_level_cm() const { return obs_cm_; }
  double etadsn() const { return etadsn_; }
  double observation_refractive_index() const { return obs_refidx_; }
  double observation_thickness() const { return obs_thick_; }
  double fast_density_factor() const { return fast_rho_fac_; }
  // Density at the equidistant grid of fast_density_factor(), increasing.
  std::span<const double> fast_density() const { return fast_rho_rev_; }
  std::span<const double> fast_horizontal() const { return fast_hori_a_; }
  std::span<const double> fast_time_coefficient() const { return fast_time_a_; }
  std::span<const double> fast_time0() const { return fast_time0_; }

  // Linear interpolation in altitude on the table levels, used for vertical rays.
  double vertical_time0_cm(double zem) const {
    require_finite(zem, "height");
    const std::size_t n = table_alt_.size();
    if (zem <= table_alt_.front())
      return table_time0_.front();
    if (zem >= table_alt_.back())
      return table_time0_.back();
    std::size_t k = 0;
    while (k + 2 < n && table_alt_[k + 1] <= zem)
      ++k;
    const double f = (zem - table_alt_[k]) / (table_alt_[k + 1] - table_alt_[k]);
    return table_time0_[k] * (1. - f) + table_time0_[k + 1] * f;
  }

private:
  // Planar ray from emlev down to olev, step in cm.
  // Outputs: horizontal displacement xo [cm], time to via [ns] and path length so [cm].
  void trace_ray_planar(double emlev, double olev, double za, double step, double &xo, double &to,
                        double &so) const {
    double ss = 0., t = 0.;
    const double vci = 1. / vacuum_speed_cm_per_ns;
    double cz = std::cos(za);
    double sz = std::sin(za);
    if (emlev < olev)
      std::swap(emlev, olev);
    double h = emlev >= olev ? emlev - olev : 0.;
    double y = olev + h;
    double x = h * std::tan(za);
    double n1 = 1. + profile_.index_minus_one_cm(y);
    const double p = n1 * sz; // Snell invariant
    while (h > 0.) {
      const double nstep = (h > step * cz) ? step : h / cz;
      const double y2 = y - nstep * cz;
      const double h2 = y2 - olev;
      const double n2 = 1. + profile_.index_minus_one_cm(y2);
      const double sz2 = p / n2;
      const double cz2 = std::sqrt(1. - sz2 * sz2);
      const double dx = nstep * 0.5 * (sz + sz2);
      const double x2 = x - dx;
      ss += nstep;
      const double dt = nstep * 0.5 * (n1 + n2) * vci;
      t += dt;
      sz = sz2;
      cz = cz2;
      x = x2;
      y = y2;
      h = h2;
      n1 = n2;
    }
    xo = x;
    to = t;
    so = ss;
  }

  AtmosphereProfile profile_;
  double obs_cm_ = 0., etadsn_ = 0., obs_refidx_ = 1., obs_thick_ = 0.;
  std::vector<double> table_alt_, table_time0_;
  double fast_rho_fac_ = 0.;
  std::vector<double> fast_rho_rev_, fast_hori_a_, fast_time_a_, fast_time0_;
};

// Per-photon ray bending for a photon emitted at height
// h_emission_m (m a.s.l.) with direction cosines (cx, cy) (downward, cz > 0), the straight-line
// position (x_m, y_m) at the observation level, and time t_ns.
//
// On return the direction is the observed one, the position is displaced horizontally, and
// the time includes the light travel time to the observation level. Returns false and leaves
// the photon unchanged when the direction is outside the unit circle (upward or horizontal
// rays, which are not bent) or when any input (or the height in cm) is not finite.
// Rejection is a return value, not an exception, so the per-photon kernel stays allocation-free.
//
// Template on the scalar type of the photon state; table lookups are evaluated in double.
template <class Scalar>
bool bend(const RefractionTables &tables, Scalar h_emission_m, Scalar &cx, Scalar &cy, Scalar &x_m,
          Scalar &y_m, Scalar &t_ns) {
  using std::asin;
  using std::cos;
  using std::sin;
  using std::sqrt;
  const auto finite = [](Scalar value) { return std::isfinite(static_cast<double>(value)); };
  if (!(finite(h_emission_m) && finite(cx) && finite(cy) && finite(x_m) && finite(y_m) &&
        finite(t_ns)))
    return false;
  const AtmosphereProfile &profile = tables.profile();
  const double zem = 100. * static_cast<double>(h_emission_m);
  if (!std::isfinite(zem))
    return false;
  const double rho = profile.density_cm(zem);
  const double de =
      tables.etadsn() * std::abs(tables.observation_thickness() - profile.thickness_cm(zem));
  const double dz = std::abs(zem - tables.observation_level_cm());
  const Scalar sin2 = cx * cx + cy * cy;
  if (!(sin2 < Scalar(1)))
    return false;
  const Scalar airmass = Scalar(1) / sqrt(Scalar(1) - sin2);
  const Scalar straight_ns =
      Scalar(dz + de) * airmass * Scalar(1. / RefractionTables::vacuum_speed_cm_per_ns);

  const auto last = static_cast<long long>(AtmosphereProfile::fast_points) - 2;
  const auto alt = tables.profile().fast_altitude_cm();
  const double fac_h = profile.fast_height_factor();
  const auto ialt = static_cast<std::size_t>(
      std::clamp(static_cast<long long>(fac_h * (zem - profile.bottom_cm())), 0LL, last));
  const double ralt = fac_h * (zem - alt[ialt]);
  const auto time0 = tables.fast_time0();
  const double time0_fast = ralt * (time0[ialt + 1] - time0[ialt]) + time0[ialt];

  const Scalar sin_em = sqrt(sin2);
  if (sin_em <= Scalar(0)) {
    // Vertical: no bending, travel time only.
    t_ns += straight_ns + Scalar(tables.vertical_time0_cm(zem));
    return true;
  }

  const double fac_rho = tables.fast_density_factor();
  const auto irho =
      static_cast<std::size_t>(std::clamp(static_cast<long long>(fac_rho * rho), 0LL, last));
  const double rrho = fac_rho * (rho - tables.fast_density()[irho]);
  const auto hori = tables.fast_horizontal();
  const auto tca = tables.fast_time_coefficient();
  const double hori_cm = rrho * (hori[irho + 1] - hori[irho]) + hori[irho];
  const double time_cm = rrho * (tca[irho + 1] - tca[irho]) + tca[irho];

  const Scalar sin_obs = sin_em * Scalar(1. + profile.index_minus_one_cm(zem)) /
                         Scalar(tables.observation_refractive_index());
  if (sin_obs > Scalar(1))
    return false;
  const Scalar theta_em = asin(sin_em);
  const Scalar theta_obs = asin(sin_obs);
  const Scalar c = cos(theta_em + Scalar(0.28) * (theta_obs - theta_em));
  const Scalar s = sin(theta_em + Scalar(0.28) * (theta_obs - theta_em));
  const Scalar c3 = c * c * c;
  const Scalar hori_off_cm = -(Scalar(hori_cm) * Scalar(hori_cm)) * s / c3;
  const Scalar travel_ns =
      Scalar(time0_fast) - (Scalar(time_cm) * Scalar(time_cm)) * (s * s) / c3 + straight_ns;

  // Displacement along the horizontal direction of the photon.
  x_m += hori_off_cm / Scalar(100) * cx / sin_em;
  y_m += hori_off_cm / Scalar(100) * cy / sin_em;
  cx = cx * sin_obs / sin_em;
  cy = cy * sin_obs / sin_em;
  t_ns += travel_ns;
  return true;
}

} // namespace dlufd
