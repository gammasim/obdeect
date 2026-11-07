#include "check.hpp"

#include "dlufd/profile.hpp"
#include "dlufd/timing.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef DLUFD_TEST_DATA_DIR
#error "DLUFD_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

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
    const auto profile = dlufd::AtmosphereProfile::from_file(
        std::string(DLUFD_TEST_DATA_DIR) + "/atmprof_ecmwf_north_winter_fixed.dat");

    // Reference values: obdeect propagation.ambient_group_index() on the same table
    // (python3.14, obdeect checkout, altitude in km and refractive_index = n-1 columns).
    CHECK_NEAR(dlufd::ambient_group_index(profile, 1000.), 1.000255944992766, 1e-12);
    CHECK_NEAR(dlufd::ambient_group_index(profile, 2180.), 1.000225191893913, 1e-12);
    CHECK_NEAR(dlufd::ambient_group_index(profile, 4500.), 1.000177406919158, 1e-12);
    CHECK_NEAR(dlufd::ambient_group_index(profile, 60000.), 1.0000000583180002, 1e-12);

    // obdeect edge rules: below the table the ground value, above the table 1.
    CHECK_NEAR(dlufd::ambient_group_index(profile, -10.), 1. + profile.table_index_minus_one()[0],
               1e-15);
    CHECK(dlufd::ambient_group_index(profile, 130000.) == 1.);
    CHECK(throws_invalid([&] { dlufd::ambient_group_index(profile, 119000.); }));

    // Arrays with invalid content are rejected as in obdeect.
    const std::vector<double> heights{0., 1000., 2000., 3000.};
    const std::vector<double> index{3e-4, 2e-4, 1e-4, 0.};
    CHECK(throws_invalid([&] { dlufd::ambient_group_index(heights, index, 100.); }));
    const std::vector<double> index_ok{3e-4, 2e-4, 1e-4, 5e-5};
    CHECK(throws_invalid([&] { dlufd::ambient_group_index(heights, index_ok, 2500.); }));
    CHECK(throws_invalid([&] {
      dlufd::ambient_group_index(std::vector<double>{0., 0., 2000., 3000.}, index_ok, 100.);
    }));

    // Propagation time: path * n / c.
    CHECK_NEAR(dlufd::propagation_time_ns(1000., 1.0002), 1000. * 1.0002 / 0.299792458, 1e-9);

    // Straight-path air delay against a fine midpoint quadrature of (n-1) / c / cos(zenith).
    const double zenith = 0.6;
    const double low = 2180., high = 20000.;
    const int steps = 200000;
    const double dh = (high - low) / steps;
    double sum = 0.;
    for (int i = 0; i < steps; ++i)
      sum += profile.refractive_index_minus_one(low + (i + 0.5) * dh) * dh;
    const double reference = sum / (dlufd::speed_of_light_m_per_ns * std::cos(zenith));
    CHECK_NEAR(dlufd::straight_path_air_delay_ns(profile, low, high, zenith) / reference, 1., 1e-6);
    CHECK_NEAR(dlufd::straight_path_air_delay_ns(profile, low, low, zenith), 0., 1e-15);
    CHECK(throws_invalid([&] { dlufd::straight_path_air_delay_ns(profile, high, low, 0.); }));
    CHECK(throws_invalid([&] { dlufd::straight_path_air_delay_ns(profile, low, high, 1.6); }));

    // Non-finite inputs throw (heights, zenith angle, group index, path length).
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (double bad : {nan, inf, -inf}) {
      CHECK(throws_invalid([&] { dlufd::ambient_group_index(profile, bad); }));
      CHECK(throws_invalid([&] { dlufd::straight_path_air_delay_ns(profile, bad, high, zenith); }));
      CHECK(throws_invalid([&] { dlufd::straight_path_air_delay_ns(profile, low, bad, zenith); }));
      CHECK(throws_invalid([&] { dlufd::straight_path_air_delay_ns(profile, low, high, bad); }));
      CHECK(throws_invalid([&] { dlufd::propagation_time_ns(bad, 1.0002); }));
      CHECK(throws_invalid([&] { dlufd::propagation_time_ns(1000., bad); }));
    }
    CHECK(throws_invalid([&] { dlufd::ambient_group_index(heights, index_ok, nan); }));
    CHECK(throws_invalid([&] { dlufd::ambient_group_index(heights, index_ok, inf); }));
  });
}
