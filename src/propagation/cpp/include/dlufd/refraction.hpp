#pragma once

// Ray bending and travel time in a plane-parallel atmosphere with n(h) = 1 + nu(h). See
// docs/refraction.md for the derivation.
//
// Snell's law n(h) sin(theta(h)) = n_e sin(theta_e) fixes the ray. With delta = n_e / n - 1,
// sin(theta) = s_e (1 + delta) and, expanded about the emission angle,
//   tan(theta) - tan(theta_e) = sum_k f_k (s_e delta)^k,
//   n / cos(theta)            = (1 + nu) sum_k g_k (s_e delta)^k,
// with Taylor coefficients f_k, g_k of s/sqrt(1-s^2) and 1/sqrt(1-s^2). delta^k is a
// polynomial in nu_e and nu, so each term reduces to products of nu_e^a and the profile
// moments M_b(h) = integral of nu^b dh, which are tabulated once per profile.

#include "dlufd/profile.hpp"
#include "dlufd/timing.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace dlufd {

class RefractionTables;
template <class T>
[[nodiscard]] bool bend(const RefractionTables &tables, T height_m, T &cx, T &cy, T &x, T &y, T &t);

class RefractionTables {
public:
  static constexpr int terms = 6; // order of the expansion in s_e * delta

  RefractionTables(AtmosphereProfile profile, double observation_level_m)
      : profile_(std::move(profile)) {
    if (!std::isfinite(observation_level_m) || observation_level_m < profile_.bottom_m() ||
        observation_level_m >= profile_.top_m()) {
      throw std::invalid_argument("observation level must lie inside the atmospheric profile");
    }
    obs_cm_ = 100. * observation_level_m;
    nu_obs_ = profile_.index_minus_one_cm(obs_cm_);
    for (int power = 0; power <= AtmosphereProfile::max_moment; ++power) {
      moment_obs_[static_cast<std::size_t>(power)] = profile_.index_moment_cm(obs_cm_, power);
    }
    // Gladstone-Dale ratio at the first level more than 1.5 km above the observation level.
    const auto &z = profile_.table_altitude_cm();
    // Without such a level the default is kept.
    etadsn_ = AtmosphereProfile::default_etadsn;
    for (std::size_t i = 0; i < z.size(); ++i) {
      if (z[i] > obs_cm_ + 1.5e5) {
        etadsn_ = profile_.table_index_minus_one()[i] / profile_.table_density()[i];
        break;
      }
    }
    build_coefficients();
  }

  [[nodiscard]] double etadsn() const noexcept { return etadsn_; }
  [[nodiscard]] double observation_level_cm() const noexcept { return obs_cm_; }
  [[nodiscard]] const AtmosphereProfile &profile() const noexcept { return profile_; }
  [[nodiscard]] double index_minus_one_at_observation() const noexcept { return nu_obs_; }

  // Air travel time of a vertical ray from the observation level to height h [cm], in
  // excess of the vacuum time [ns].
  double vertical_time0_cm(double h) const {
    if (!std::isfinite(h)) {
      throw std::invalid_argument("height must be finite");
    }
    return (moment(h, 1)) / (100. * speed_of_light_m_per_ns);
  }

  // Integral of nu^power from the observation level to h [cm].
  [[nodiscard]] double moment(double h_cm, int power) const {
    if (power < 0 || power > AtmosphereProfile::max_moment) {
      throw std::invalid_argument("moment order must be in [0, max_moment]");
    }
    return profile_.index_moment_cm(h_cm, power) - moment_obs_[static_cast<std::size_t>(power)];
  }

private:
  template <class T> friend bool bend(const RefractionTables &, T, T &, T &, T &, T &, T &);

  // Coefficient of nu_e^a * M_b in the expansion of (nu_e - nu)^k (1 + nu)^m, m = -k (position)
  // or 1 - k (time), truncated at total degree `terms`.
  [[nodiscard]] double position_coefficient(int k, int a, int b) const {
    return position_[static_cast<std::size_t>(k)][static_cast<std::size_t>(a)]
                    [static_cast<std::size_t>(b)];
  }
  [[nodiscard]] double time_coefficient(int k, int a, int b) const {
    return time_[static_cast<std::size_t>(k)][static_cast<std::size_t>(a)]
                [static_cast<std::size_t>(b)];
  }

  static constexpr std::size_t dim = terms + 1;
  using Table = std::array<std::array<std::array<double, dim>, dim>, dim>;

  // Binomial coefficient of x^i in (1 + x)^m for integer m (also negative).
  static double binomial(int m, int i) {
    double value = 1.;
    for (int j = 0; j < i; ++j) {
      value *= double(m - j) / double(j + 1);
    }
    return value;
  }

  static void fill(Table &table, int shift) {
    for (auto &a : table)
      for (auto &b : a)
        b.fill(0.);
    for (int k = 1; k <= terms; ++k) {
      const int m = shift - k;
      for (int j = 0; j <= k; ++j) {
        // (nu_e - nu)^k contributes C(k, j) nu_e^(k-j) (-nu)^j
        double c = 1.;
        for (int q = 0; q < j; ++q) {
          c *= double(k - q) / double(q + 1);
        }
        c *= (j % 2 == 0) ? 1. : -1.;
        for (int i = 0; k + i <= terms; ++i) {
          table[static_cast<std::size_t>(k)][static_cast<std::size_t>(k - j)]
               [static_cast<std::size_t>(j + i)] += c * binomial(m, i);
        }
      }
    }
  }

  void build_coefficients() {
    fill(position_, 0);
    fill(time_, 1);
  }

  AtmosphereProfile profile_;
  double obs_cm_{};
  double nu_obs_{};
  double etadsn_{};
  std::array<double, AtmosphereProfile::max_moment + 1> moment_obs_{};
  Table position_{};
  Table time_{};
};

// Bends a photon emitted at height_m [m a.s.l.] with direction cosines (cx, cy) (downward)
// to the observation level of the tables. On input (x, y) is the straight-line position at the
// observation level and t the time at emission; on output the position is displaced by the
// bending, (cx, cy) is the direction at the observation level and t has the travel time added
// (the straight vacuum path along the emission direction plus the air excess).
//
// Returns false and leaves every argument unchanged if an input is not finite, the direction
// is not strictly inside the unit circle, the emission is below the observation level, or the
// ray is too close to the horizon for the expansion (s_e (nu_obs - nu_e) / cos^2 > 0.1).
// Arithmetic is done in double.
template <class T>
[[nodiscard]] bool bend(const RefractionTables &tables, T height_m, T &cx, T &cy, T &x, T &y,
                        T &t) {
  const double h_m = static_cast<double>(height_m);
  const double ux = static_cast<double>(cx), uy = static_cast<double>(cy);
  if (!std::isfinite(h_m) || !std::isfinite(ux) || !std::isfinite(uy) ||
      !std::isfinite(static_cast<double>(x)) || !std::isfinite(static_cast<double>(y)) ||
      !std::isfinite(static_cast<double>(t))) {
    return false;
  }
  const double h = 100. * h_m;
  const double s2 = ux * ux + uy * uy;
  if (!std::isfinite(h) || !(s2 < 1.) || h < tables.observation_level_cm()) {
    return false;
  }
  const double s_e = std::sqrt(s2);
  const double c_e = std::sqrt(1. - s2);
  const double nu_e = tables.profile().index_minus_one_cm(h);
  const double nu_obs = tables.index_minus_one_at_observation();
  const double height = h - tables.observation_level_cm(); // cm
  if (s_e * std::fabs(nu_obs - nu_e) > 0.1 * c_e * c_e) {
    return false;
  }
  const double s_obs = s_e * (1. + nu_e) / (1. + nu_obs);
  if (!(s_obs < 1.)) {
    return false; // total internal reflection
  }

  // Moments M_b = integral of nu^b over [obs, h], b = 0 is the height.
  std::array<double, AtmosphereProfile::max_moment + 1> m{};
  m[0] = height;
  for (int b = 1; b <= AtmosphereProfile::max_moment; ++b) {
    m[static_cast<std::size_t>(b)] = tables.moment(h, b);
  }
  std::array<double, RefractionTables::terms + 1> nu_power{};
  nu_power[0] = 1.;
  for (std::size_t a = 1; a < nu_power.size(); ++a) {
    nu_power[a] = nu_power[a - 1] * nu_e;
  }

  // Taylor coefficients g_k (of 1/sqrt(1-s^2)) and f_k = s_e g_k + g_(k-1) at s_e.
  std::array<double, RefractionTables::terms + 1> g{};
  g[0] = 1. / c_e;
  double previous = 0.;
  for (int k = 0; k < RefractionTables::terms; ++k) {
    const double next =
        (s_e * double(2 * k + 1) * g[static_cast<std::size_t>(k)] + double(k) * previous) /
        (double(k + 1) * c_e * c_e);
    previous = g[static_cast<std::size_t>(k)];
    g[static_cast<std::size_t>(k + 1)] = next;
  }

  double displacement = 0.;    // integral of (tan(theta) - tan(theta_e)) dh [cm]
  double excess = m[1] * g[0]; // integral of nu g_0 dh [cm]
  double s_power = 1.;
  for (int k = 1; k <= RefractionTables::terms; ++k) {
    s_power *= s_e;
    const double f_k = s_e * g[static_cast<std::size_t>(k)] + g[static_cast<std::size_t>(k - 1)];
    double position_sum = 0., time_sum = 0.;
    for (int a = 0; a <= k; ++a) {
      for (int b = 0; a + b <= RefractionTables::terms; ++b) {
        const double mb = m[static_cast<std::size_t>(b)];
        position_sum +=
            tables.position_coefficient(k, a, b) * nu_power[static_cast<std::size_t>(a)] * mb;
        time_sum += tables.time_coefficient(k, a, b) * nu_power[static_cast<std::size_t>(a)] * mb;
      }
    }
    displacement += f_k * s_power * position_sum;
    excess += g[static_cast<std::size_t>(k)] * s_power * time_sum;
  }

  const double travel_ns = (height * g[0] + excess) / (100. * speed_of_light_m_per_ns);
  if (s_e > 0.) {
    const double shift_m = displacement / 100.;
    x = static_cast<T>(static_cast<double>(x) + shift_m * ux / s_e);
    y = static_cast<T>(static_cast<double>(y) + shift_m * uy / s_e);
    const double scale = s_obs / s_e;
    cx = static_cast<T>(ux * scale);
    cy = static_cast<T>(uy * scale);
  }
  t = static_cast<T>(static_cast<double>(t) + travel_ns);
  return true;
}

} // namespace dlufd
