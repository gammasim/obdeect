#include "check.hpp"

#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/refraction_horizon.hpp"
#include "obdeect/stream/constants.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using obdeect::propagation::bend;
using obdeect::propagation::bend_quadrature;
using obdeect::propagation::bend_quadrature_pieces;
using obdeect::propagation::bend_series_limit_exceeded;
using obdeect::propagation::RefractionTables;

constexpr double kPi = obdeect::stream::pi;
constexpr double kObservationM = 2180.;

struct Photon {
  double cx, cy, x, y, t;
};

// Photon with zenith angle theta (degrees) in the x-z plane, at the origin, time zero.
Photon start(double zenith_deg) { return {std::sin(zenith_deg * kPi / 180.), 0., 0., 0., 0.}; }

bool throws_invalid(auto &&action) {
  try {
    action();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

} // namespace

int main() {
  return run_checks([] {
    const auto profile = obdeect::propagation::AtmosphereProfile::from_file(
        std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/atmprof_ecmwf_north_winter_fixed.dat");
    const RefractionTables tables(profile, kObservationM);
    const std::array<double, 3> emissions{10000., 20000., 40000.};

    // Agreement with the series bend() for zenith 0..82 degrees.
    double max_dx = 0., max_dt = 0., max_dir = 0.;
    for (double h : emissions) {
      for (int zenith = 0; zenith <= 82; zenith += 2) {
        Photon a = start(zenith), b = start(zenith);
        CHECK(bend(tables, h, a.cx, a.cy, a.x, a.y, a.t));
        CHECK(bend_quadrature(tables, h, b.cx, b.cy, b.x, b.y, b.t));
        max_dx = std::max(max_dx, std::fabs(a.x - b.x));
        max_dt = std::max(max_dt, std::fabs(a.t - b.t));
        max_dir = std::max({max_dir, std::fabs(a.cx - b.cx), std::fabs(a.cy - b.cy)});
        CHECK_NEAR(a.x, b.x, 1e-4);
        CHECK_NEAR(a.t, b.t, 1e-4);
        CHECK_NEAR(a.cx, b.cx, 1e-12);
        CHECK_NEAR(a.cy, b.cy, 1e-12);
      }
    }
    std::printf("bend vs quadrature, zenith 0..82, 10/20/40 km: max |dx| = %.3g m, "
                "max |dt| = %.3g ns, max |dcos| = %.3g\n",
                max_dx, max_dt, max_dir);

    // Dispersion factors: bend_quadrature(k, k_group) agrees with bend(k, k_group) where the series
    // is valid (zenith 0..82 degrees), and the scaled quadrature is what the stage uses beyond it.
    {
      const double k = 1.03, k_group = 1.02;
      double dispersed_dx = 0., dispersed_dt = 0.;
      for (double h : emissions) {
        for (int zenith = 0; zenith <= 82; zenith += 2) {
          Photon a = start(zenith), b = start(zenith);
          CHECK(bend(tables, h, a.cx, a.cy, a.x, a.y, a.t, k, k_group));
          CHECK(bend_quadrature(tables, h, b.cx, b.cy, b.x, b.y, b.t, k, k_group));
          dispersed_dx = std::max(dispersed_dx, std::fabs(a.x - b.x));
          dispersed_dt = std::max(dispersed_dt, std::fabs(a.t - b.t));
          CHECK_NEAR(a.x, b.x, 1e-4);
          // The series uses a straight-path group correction; quadrature follows the bent ray.
          constexpr int steps = 10000;
          const double step = (h - kObservationM) / steps;
          const double sin_e = start(zenith).cx;
          const double nu_e = profile.index_minus_one_cm(100. * h);
          double group_integral = 0.;
          for (int i = 0; i <= steps; ++i) {
            const double nu = profile.index_minus_one_cm(100. * (kObservationM + i * step));
            const double sin_h = sin_e * (1. + k * nu_e) / (1. + k * nu);
            const double weight = i == 0 || i == steps ? 1. : (i % 2 == 0 ? 2. : 4.);
            group_integral += weight * nu / std::sqrt(1. - sin_h * sin_h);
          }
          const double bent_m = group_integral * step / 3.;
          const double straight_m =
              tables.moment(100. * h, 1) / (100. * std::cos(zenith * kPi / 180.));
          const double correction =
              (k_group - k) * (bent_m - straight_m) / obdeect::stream::speed_of_light_m_per_ns;
          CHECK_NEAR(a.t + correction, b.t, 1e-4);
        }
      }
      // Beyond the series limit the scaled quadrature still gives the photon's bent position.
      Photon far = start(88.);
      CHECK(!bend(tables, 20000., far.cx, far.cy, far.x, far.y, far.t, k, k_group));
      CHECK(bend_series_limit_exceeded(tables, 20000., far.cx, far.cy, k));
      CHECK(bend_quadrature(tables, 20000., far.cx, far.cy, far.x, far.y, far.t, k, k_group));
      CHECK(far.x < 0. && std::isfinite(far.t));
      // k = 1 is the dispersion-free quadrature bit for bit.
      Photon plain = start(88.), unit = start(88.);
      CHECK(bend_quadrature(tables, 20000., plain.cx, plain.cy, plain.x, plain.y, plain.t));
      CHECK(bend_quadrature(tables, 20000., unit.cx, unit.cy, unit.x, unit.y, unit.t, 1., 1.));
      CHECK(plain.x == unit.x && plain.t == unit.t && plain.cx == unit.cx);
      std::printf("dispersion (k 1.03, k_group 1.02): quadrature vs series max |dx| = %.3g m, "
                  "|dt| = %.3g ns\n",
                  dispersed_dx, dispersed_dt);
    }

    // Group delay follows the bent ray, including near the horizon.
    {
      constexpr double h = 20000., k = 1.03, kg = 1.02;
      Photon phase = start(89.), group = start(89.);
      CHECK(bend_quadrature(tables, h, phase.cx, phase.cy, phase.x, phase.y, phase.t, k, k));
      CHECK(bend_quadrature(tables, h, group.cx, group.cy, group.x, group.y, group.t, k, kg));
      constexpr int intervals = 200000;
      const double step = (h - kObservationM) / intervals;
      const double sin_e = start(89.).cx;
      const double nu_e = profile.index_minus_one_cm(100. * h);
      double integral = 0.;
      for (int i = 0; i <= intervals; ++i) {
        const double nu = profile.index_minus_one_cm(100. * (kObservationM + i * step));
        const double sin_h = sin_e * (1. + k * nu_e) / (1. + k * nu);
        const double f = nu / std::sqrt(1. - sin_h * sin_h);
        integral += (i == 0 || i == intervals ? 1. : (i % 2 == 0 ? 2. : 4.)) * f;
      }
      const double expected =
          (kg - k) * integral * step / (3. * obdeect::stream::speed_of_light_m_per_ns);
      CHECK_NEAR(group.t - phase.t, expected, 0.001);
    }

    // Convergence: halving the largest sub-interval (250 m -> 125 m) changes the result by less
    // than 1e-4 m and 1e-3 ns near the horizon.
    double conv_dx = 0., conv_dt = 0., conv_dx_fine = 0., conv_dt_fine = 0.;
    for (double h : emissions) {
      for (double zenith : {87., 88., 89.}) {
        Photon coarse = start(zenith), mid = start(zenith), fine = start(zenith);
        CHECK(bend_quadrature_pieces(tables, 250., h, coarse.cx, coarse.cy, coarse.x, coarse.y,
                                     coarse.t));
        CHECK(bend_quadrature_pieces(tables, 125., h, mid.cx, mid.cy, mid.x, mid.y, mid.t));
        CHECK(bend_quadrature_pieces(tables, 62.5, h, fine.cx, fine.cy, fine.x, fine.y, fine.t));
        conv_dx = std::max(conv_dx, std::fabs(coarse.x - mid.x));
        conv_dt = std::max(conv_dt, std::fabs(coarse.t - mid.t));
        conv_dx_fine = std::max(conv_dx_fine, std::fabs(mid.x - fine.x));
        conv_dt_fine = std::max(conv_dt_fine, std::fabs(mid.t - fine.t));
        CHECK(std::fabs(coarse.x - mid.x) < 1e-4);
        CHECK(std::fabs(coarse.t - mid.t) < 1e-3);
        CHECK(coarse.cx == mid.cx && coarse.cy == mid.cy);
      }
    }
    std::printf("convergence 87..89 deg: 250->125 m max |dx| = %.3g m, |dt| = %.3g ns; "
                "125->62.5 m max |dx| = %.3g m, |dt| = %.3g ns\n",
                conv_dx, conv_dt, conv_dx_fine, conv_dt_fine);

    // Photons beyond the series limit: bend() rejects them, bend_quadrature() bends them. Zenith
    // 84..89 degrees; the predicate decides which case applies, and bend() must reject exactly
    // the photons it flags.
    int beyond = 0;
    for (double h : emissions) {
      for (int zenith = 84; zenith <= 89; ++zenith) {
        Photon s = start(zenith);
        const bool series = bend(tables, h, s.cx, s.cy, s.x, s.y, s.t);
        const bool limit =
            bend_series_limit_exceeded(tables, h, start(zenith).cx, start(zenith).cy);
        CHECK(series == !limit);
        if (series) {
          continue;
        }
        ++beyond;
        Photon q = start(zenith);
        const double cx_in = q.cx;
        CHECK(bend_quadrature(tables, h, q.cx, q.cy, q.x, q.y, q.t));
        CHECK(std::isfinite(q.x) && std::isfinite(q.y) && std::isfinite(q.t));
        CHECK(q.x < 0.);
        CHECK(q.cx > 0. && q.cx < cx_in);
        CHECK(q.cy == 0.);
        CHECK_NEAR(q.cx * q.cx + q.cy * q.cy + (1. - q.cx * q.cx - q.cy * q.cy), 1., 1e-15);
      }
    }
    CHECK(beyond > 0);
    std::printf("photons beyond the series limit (84..89 deg, 10/20/40 km): %d cases\n", beyond);

    // A float photon gives the same result as double to single precision.
    {
      // Both calls start from the same float-rounded direction; the float result is then
      // accurate to float precision (x ~ 3e4 m, t ~ 1.6e6 ns).
      float fcx = static_cast<float>(std::sin(88. * kPi / 180.)), fcy = 0.f, fx = 0.f, fy = 0.f,
            ft = 0.f;
      Photon d{static_cast<double>(fcx), 0., 0., 0., 0.};
      CHECK(bend_quadrature(tables, 20000., d.cx, d.cy, d.x, d.y, d.t));
      CHECK(bend_quadrature(tables, 20000.f, fcx, fcy, fx, fy, ft));
      CHECK_NEAR(static_cast<double>(fx), d.x, 1e-2);
      CHECK_NEAR(static_cast<double>(ft), d.t, 0.25);
      CHECK_NEAR(static_cast<double>(fcx), d.cx, 1e-7);
    }

    // Invalid inputs are rejected, and the photon is left unchanged, as in bend().
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (double bad : {nan, inf, -inf}) {
      for (std::size_t field = 0; field < 6; ++field) {
        std::array<double, 6> v{20000., std::sin(0.7), 0., 1., 2., 3.}; // h, cx, cy, x, y, t
        v[field] = bad;
        const std::array<double, 6> before = v;
        CHECK(!bend_quadrature(tables, v[0], v[1], v[2], v[3], v[4], v[5]));
        for (std::size_t k = 0; k < v.size(); ++k) {
          CHECK(v[k] == before[k] || (std::isnan(v[k]) && std::isnan(before[k])));
        }
      }
    }
    {
      // Height whose conversion to cm overflows.
      std::array<double, 6> v{std::numeric_limits<double>::max(), std::sin(0.7), 0., 1., 2., 3.};
      const std::array<double, 6> before = v;
      CHECK(!bend_quadrature(tables, v[0], v[1], v[2], v[3], v[4], v[5]));
      CHECK(v == before);
    }
    {
      // Direction upward or horizontal, and emission below the observation level.
      Photon upward{1.1, 0., 0., 0., 0.};
      CHECK(!bend_quadrature(tables, 20000., upward.cx, upward.cy, upward.x, upward.y, upward.t));
      Photon horizontal{1., 0., 0., 0., 0.};
      CHECK(!bend_quadrature(tables, 20000., horizontal.cx, horizontal.cy, horizontal.x,
                             horizontal.y, horizontal.t));
      Photon below = start(88.);
      const Photon below_before = below;
      CHECK(!bend_quadrature(tables, 1000., below.cx, below.cy, below.x, below.y, below.t));
      CHECK(below.cx == below_before.cx && below.x == below_before.x && below.t == below_before.t);
      CHECK(!bend_series_limit_exceeded(tables, 1000., below.cx, below.cy));
      CHECK(!bend_series_limit_exceeded(tables, nan, below.cx, below.cy));
    }
    // Invalid piece sizes are rejected.
    {
      Photon p = start(88.);
      CHECK(throws_invalid([&] {
        static_cast<void>(bend_quadrature_pieces(tables, 0., 20000., p.cx, p.cy, p.x, p.y, p.t));
      }));
      CHECK(throws_invalid([&] {
        static_cast<void>(bend_quadrature_pieces(tables, nan, 20000., p.cx, p.cy, p.x, p.y, p.t));
      }));
    }
  });
}
