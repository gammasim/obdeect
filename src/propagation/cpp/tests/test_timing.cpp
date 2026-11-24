// Light travel time and the refractive index at the observation level (timing.hpp): analytic
// cases (constant and exponential n-1, closed-form air delay), regression values on the ECMWF
// fixture, the edge rules and the rejections.
#include "check.hpp"
#include "stream_fixtures.hpp"

#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/timing.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using obdeect_test::kInf;
using obdeect_test::kNan;

// Exponential atmosphere, levels every 1 km from 0 to 10 km (as in test_profile.cpp):
// rho = rho0 exp(-h / H) g/cm^3, H = 7 km, thickness = int from h to 10 km of rho (g/cm^2).
constexpr double exponential_rho0 = 1.2e-3;
constexpr double exponential_scale_m = 7000.;

std::string exponential_profile() {
  std::ostringstream text;
  text << "# exponential test profile\n";
  const double rho_top = exponential_rho0 * std::exp(-10000. / exponential_scale_m);
  for (int km = 0; km <= 10; ++km) {
    const double rho = exponential_rho0 * std::exp(-km * 1000. / exponential_scale_m);
    text << km << ' ' << rho << ' ' << 7e5 * (rho - rho_top) << ' '
         << obdeect::propagation::AtmosphereProfile::default_etadsn * rho << '\n';
  }
  return text.str();
}

// Levels of an exponential index n-1 = a exp(-h / scale) at heights 0, 1, ..., 6 km.
std::vector<double> exponential_heights() { return {0., 1000., 2000., 3000., 4000., 5000., 6000.}; }

std::vector<double> exponential_index(double a, double scale_m) {
  std::vector<double> index;
  for (const double h : exponential_heights()) {
    index.push_back(a * std::exp(-h / scale_m));
  }
  return index;
}

} // namespace

int main() {
  return run_checks([] {
    // Constant n-1: the log spline is constant, so n_g = 1 + (n-1) at every height in range.
    const std::vector<double> levels{0., 1000., 2000., 3000.};
    const std::vector<double> constant(levels.size(), 3e-4);
    for (const double h : {0., 1.5, 1000., 1999.}) {
      CHECK_NEAR(obdeect::propagation::ambient_group_index(levels, constant, h) / (1. + 3e-4), 1.,
                 1e-14);
    }

    // Exponential n-1 = a exp(-h / H): the log is linear, which the natural spline reproduces
    // exactly, so n_g = 1 + a exp(-h / H) at any height below the second-highest level.
    const double a = 3e-4, scale = 8000.;
    for (const double h : {0., 350., 1234.5, 2500., 4999.}) {
      CHECK_NEAR(obdeect::propagation::ambient_group_index(exponential_heights(),
                                                           exponential_index(a, scale), h) /
                     (1. + a * std::exp(-h / scale)),
                 1., 1e-12);
    }

    // Edge rules: below the lowest level the lowest value, at or above the highest level 1.
    CHECK(obdeect::propagation::ambient_group_index(levels, constant, -10.) == 1. + 3e-4);
    CHECK(obdeect::propagation::ambient_group_index(levels, constant, 3000.) == 1.);

    // Closed form of the straight-path air delay in the exponential atmosphere of
    // exponential_profile() (n-1 = etadsn rho, rho = rho0 exp(-h / H), H = 7 km):
    // (1/c) int (n-1) dh / cos(zenith) = etadsn rho0 H (exp(-low/H) - exp(-high/H)) / (c cos).
    const double low = 2000., high = 6000., zenith = 0.6;
    const auto exponential = obdeect::propagation::AtmosphereProfile(exponential_profile());
    const double closed =
        obdeect::propagation::AtmosphereProfile::default_etadsn * exponential_rho0 *
        exponential_scale_m *
        (std::exp(-low / exponential_scale_m) - std::exp(-high / exponential_scale_m)) /
        (obdeect::propagation::speed_of_light_m_per_ns * std::cos(zenith));
    // The profile interpolates its thickness column (exponential minus a constant) in log space,
    // which leaves a relative error of about 1e-6 at these heights.
    CHECK_NEAR(obdeect::propagation::straight_path_air_delay_ns(exponential, low, high, zenith) /
                   closed,
               1., 1e-5);

    // Regression values on the ECMWF fixture (obdeect output, python3.14 with the altitude in km
    // and the refractive index as n-1; not an independent reference).
    const auto profile = obdeect::propagation::AtmosphereProfile::from_file(
        std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/atmprof_ecmwf_north_winter_fixed.dat");
    CHECK_NEAR(obdeect::propagation::ambient_group_index(profile, 1000.), 1.000255944992766, 1e-12);
    CHECK_NEAR(obdeect::propagation::ambient_group_index(profile, 2180.), 1.000225191893913, 1e-12);
    CHECK_NEAR(obdeect::propagation::ambient_group_index(profile, 4500.), 1.000177406919158, 1e-12);
    CHECK_NEAR(obdeect::propagation::ambient_group_index(profile, 60000.), 1.0000000583180002,
               1e-12);

    // Below the table the ground value, above the table 1; the profile height limit is rejected.
    CHECK_NEAR(obdeect::propagation::ambient_group_index(profile, -10.),
               1. + profile.table_index_minus_one()[0], 1e-15);
    CHECK(obdeect::propagation::ambient_group_index(profile, 130000.) == 1.);
    CHECK(throws_with([&] { obdeect::propagation::ambient_group_index(profile, 119000.); },
                      "height must lie below the second-highest level"));

    // Propagation time: path * n / c.
    CHECK_NEAR(obdeect::propagation::propagation_time_ns(1000., 1.0002),
               1000. * 1.0002 / 0.299792458, 1e-9);

    // Rejected levels, heights and zenith angles.
    const std::vector<double> index_ok{3e-4, 2e-4, 1e-4, 5e-5};
    CHECK(throws_with(
        [&] {
          obdeect::propagation::ambient_group_index(
              levels, std::vector<double>{3e-4, 2e-4, 1e-4, 0.}, 100.);
        },
        "n-1 must be positive and finite"));
    CHECK(throws_with([&] { obdeect::propagation::ambient_group_index(levels, index_ok, 2500.); },
                      "height must lie below the second-highest level"));
    CHECK(throws_with(
        [&] {
          obdeect::propagation::ambient_group_index(std::vector<double>{0., 0., 2000., 3000.},
                                                    index_ok, 100.);
        },
        "heights must be finite and strictly ascending"));
    CHECK(throws_with(
        [&] { obdeect::propagation::straight_path_air_delay_ns(profile, high, low, 0.); },
        "air delay needs low <= high"));
    CHECK(throws_with(
        [&] { obdeect::propagation::straight_path_air_delay_ns(profile, low, high, 1.6); },
        "air delay needs 0 <= zenith < pi/2"));

    // Non-finite inputs throw (heights, zenith angle, group index, path length).
    for (const double bad : {kNan, kInf, -kInf}) {
      CHECK(throws_with([&] { obdeect::propagation::ambient_group_index(profile, bad); },
                        "height must be finite"));
      CHECK(throws_with(
          [&] { obdeect::propagation::straight_path_air_delay_ns(profile, bad, high, zenith); },
          "needs finite heights"));
      CHECK(throws_with(
          [&] { obdeect::propagation::straight_path_air_delay_ns(profile, low, bad, zenith); },
          "needs finite heights"));
      CHECK(throws_with(
          [&] { obdeect::propagation::straight_path_air_delay_ns(profile, low, high, bad); },
          "finite heights and zenith angle"));
      CHECK(throws_with([&] { obdeect::propagation::propagation_time_ns(bad, 1.0002); },
                        "finite path and index"));
      CHECK(throws_with([&] { obdeect::propagation::propagation_time_ns(1000., bad); },
                        "finite path and index"));
      CHECK(throws_with([&] { obdeect::propagation::ambient_group_index(levels, index_ok, bad); },
                        "height must be finite"));
    }
  });
}
