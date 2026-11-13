// Analytic and first-principles tests of the profile, ray bending, travel time and extinction.
// The ray integrals are checked against a direct numerical integration of the exact integrands
// tan(theta(h)) and n(h)/cos(theta(h)) with Snell's law, using only the profile lookups.
#include "check.hpp"

#include "obdeect/propagation/extinction.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/timing.hpp"
#include "obdeect/propagation/transmission.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kScaleHeightKm = 8.0;
constexpr double kNu0 = 2.9e-4;
constexpr double kRho0 = 1.25e-3;

// Exponential atmosphere to 80 km in 2 km levels; n-1 = kNu0 exp(-h/H), thickness = integral
// of the density from h to the top.
std::string exponential_text(double nu_scale = 1.0) {
  std::ostringstream text;
  text.precision(12);
  const double top = 80.;
  const double rho_top = kRho0 * std::exp(-top / kScaleHeightKm);
  for (int km = 0; km <= 80; km += 2) {
    const double rho = kRho0 * std::exp(-km / kScaleHeightKm);
    text << km << ' ' << rho << ' ' << kScaleHeightKm * 1e5 * (rho - rho_top) << ' '
         << nu_scale * kNu0 * std::exp(-km / kScaleHeightKm) << '\n';
  }
  return text.str();
}

// Constant n-1 with decreasing density.
std::string constant_index_text(double nu) {
  std::ostringstream text;
  text.precision(12);
  for (int km = 0; km <= 20; km += 2) {
    const double rho = kRho0 * std::exp(-km / kScaleHeightKm);
    const double rho_top = kRho0 * std::exp(-20. / kScaleHeightKm);
    text << km << ' ' << rho << ' ' << kScaleHeightKm * 1e5 * (rho - rho_top) << ' ' << nu << '\n';
  }
  return text.str();
}

double simpson(const std::function<double(double)> &f, double a, double b, int n) {
  const double h = (b - a) / n;
  double sum = f(a) + f(b);
  for (int i = 1; i < n; ++i) {
    sum += f(a + i * h) * (i % 2 == 0 ? 2. : 4.);
  }
  return sum * h / 3.;
}

struct Result {
  double x, y, cx, cy, t;
};

Result bend_from(const obdeect::propagation::RefractionTables &tables, double height_m, double zenith,
                 double azimuth) {
  const double s = std::sin(zenith);
  double cx = s * std::cos(azimuth), cy = s * std::sin(azimuth), x = 0., y = 0., t = 0.;
  CHECK(obdeect::propagation::bend(tables, height_m, cx, cy, x, y, t));
  return {x, y, cx, cy, t};
}

// Direct integration of Snell's law: displacement relative to the straight line [m] and
// travel time minus the straight vacuum path [ns].
void direct_integrals(const obdeect::propagation::AtmosphereProfile &profile, double obs, double h_e,
                      double zenith, double &displacement, double &excess_ns) {
  const double nu_e = profile.refractive_index_minus_one(h_e);
  const double invariant = (1. + nu_e) * std::sin(zenith);
  auto sine = [&](double h) { return invariant / (1. + profile.refractive_index_minus_one(h)); };
  auto tangent = [&](double h) {
    const double s = sine(h);
    return s / std::sqrt(1. - s * s);
  };
  auto slowness = [&](double h) {
    const double s = sine(h);
    return (1. + profile.refractive_index_minus_one(h)) / std::sqrt(1. - s * s);
  };
  displacement = simpson(tangent, obs, h_e, 20000) - (h_e - obs) * std::tan(zenith);
  excess_ns = (simpson(slowness, obs, h_e, 20000) - (h_e - obs) / std::cos(zenith)) /
              obdeect::propagation::speed_of_light_m_per_ns;
}

void test_profile_integrals() {
  const obdeect::propagation::AtmosphereProfile profile(exponential_text());
  const double scale_cm = kScaleHeightKm * 1e5;
  for (double km : {0., 0.7, 3.3, 10., 55.5, 79.9}) {
    const double h = km * 1e5;
    const double top = 80. * 1e5;
    const double nu_integral =
        kNu0 * scale_cm * (std::exp(-h / scale_cm) - std::exp(-top / scale_cm));
    CHECK_NEAR(profile.index_integral_cm(h) / nu_integral, 1., 1e-9);
    const double rho_top = kRho0 * std::exp(-80. / kScaleHeightKm);
    const double rho = kRho0 * std::exp(-km / kScaleHeightKm);
    CHECK_NEAR(profile.thickness_cm(h) / (scale_cm * (rho - rho_top)), 1., 1e-9);
    CHECK_NEAR(profile.density_cm(h) / rho, 1., 1e-9);
  }
  // Monotonic in height.
  double previous = profile.index_integral_cm(0.);
  for (double h = 1000.; h < 80e5; h += 7919.) {
    const double current = profile.index_integral_cm(h);
    CHECK(current < previous);
    previous = current;
  }
}

void test_constant_index_is_straight() {
  const double nu = 2e-4;
  const obdeect::propagation::AtmosphereProfile profile(constant_index_text(nu));
  const double obs = 1000.;
  const obdeect::propagation::RefractionTables tables(profile, obs);
  for (double zenith_deg : {0., 10., 45., 70.}) {
    const double zenith = zenith_deg * kPi / 180.;
    const auto r = bend_from(tables, 15000., zenith, 0.4);
    CHECK_NEAR(r.x, 0., 1e-9);
    CHECK_NEAR(r.y, 0., 1e-9);
    CHECK_NEAR(r.cx, std::sin(zenith) * std::cos(0.4), 1e-15);
    CHECK_NEAR(r.cy, std::sin(zenith) * std::sin(0.4), 1e-15);
    CHECK_NEAR(r.t, 14000. * (1. + nu) / (obdeect::propagation::speed_of_light_m_per_ns * std::cos(zenith)), 1e-9);
  }
}

void test_vertical_ray() {
  const obdeect::propagation::AtmosphereProfile profile(exponential_text());
  const obdeect::propagation::RefractionTables tables(profile, 2000.);
  const auto r = bend_from(tables, 30000., 0., 0.);
  const double nu_integral = kNu0 * kScaleHeightKm * 1e3 *
                             (std::exp(-2. / kScaleHeightKm) - std::exp(-30. / kScaleHeightKm));
  CHECK_NEAR(r.x, 0., 0.);
  CHECK_NEAR(r.t - 28000. / obdeect::propagation::speed_of_light_m_per_ns,
             nu_integral / obdeect::propagation::speed_of_light_m_per_ns, 1e-9);
}

void test_against_direct_integration() {
  const obdeect::propagation::AtmosphereProfile profile(exponential_text());
  const double obs = 2000.;
  const obdeect::propagation::RefractionTables tables(profile, obs);
  for (double h_e : {5000., 12000., 30000., 60000.}) {
    for (double zenith_deg : {2., 20., 45., 60., 75., 82.}) {
      const double zenith = zenith_deg * kPi / 180.;
      double displacement = 0., excess = 0.;
      direct_integrals(profile, obs, h_e, zenith, displacement, excess);
      const auto r = bend_from(tables, h_e, zenith, 0.);
      CHECK_NEAR(r.x, displacement, 1e-7 + 1e-8 * std::fabs(displacement));
      const double straight_ns = (h_e - obs) / (obdeect::propagation::speed_of_light_m_per_ns * std::cos(zenith));
      CHECK_NEAR(r.t - straight_ns, excess, 1e-7 + 1e-8 * std::fabs(excess));
      // Snell's law at the observation level.
      const double nu_e = profile.refractive_index_minus_one(h_e);
      const double nu_obs = profile.refractive_index_minus_one(obs);
      CHECK_NEAR(std::hypot(r.cx, r.cy), std::sin(zenith) * (1. + nu_e) / (1. + nu_obs), 1e-14);
    }
  }
}

void test_real_profile_against_direct_integration() {
  const auto profile = obdeect::propagation::AtmosphereProfile::from_file(std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) +
                                                           "/atmprof_ecmwf_north_winter_fixed.dat");
  const double obs = 2180.;
  const obdeect::propagation::RefractionTables tables(profile, obs);
  for (double h_e : {10000., 20000., 40000.}) {
    for (double zenith_deg : {20., 60.}) {
      const double zenith = zenith_deg * kPi / 180.;
      double displacement = 0., excess = 0.;
      direct_integrals(profile, obs, h_e, zenith, displacement, excess);
      const auto r = bend_from(tables, h_e, zenith, 0.);
      CHECK_NEAR(r.x, displacement, 1e-6);
      CHECK_NEAR(r.t - (h_e - obs) / (obdeect::propagation::speed_of_light_m_per_ns * std::cos(zenith)), excess,
                 1e-6);
    }
  }
}

void test_small_angle_limit() {
  // First order in the zenith angle: displacement = -theta * integral (nu - nu_e) dh.
  const obdeect::propagation::AtmosphereProfile profile(exponential_text());
  const double obs = 1000., h_e = 20000., theta = 1e-3;
  const obdeect::propagation::RefractionTables tables(profile, obs);
  const auto r = bend_from(tables, h_e, theta, 0.);
  const double nu_e = kNu0 * std::exp(-h_e / (kScaleHeightKm * 1e3));
  const double integral =
      kNu0 * kScaleHeightKm * 1e3 *
          (std::exp(-obs / (kScaleHeightKm * 1e3)) - std::exp(-h_e / (kScaleHeightKm * 1e3))) -
      nu_e * (h_e - obs);
  CHECK_NEAR(r.x / (-theta * integral), 1., 1e-3); // neglects terms of order n-1
  CHECK(r.x < 0.);
}

void test_symmetry_and_monotonicity() {
  const obdeect::propagation::AtmosphereProfile profile(exponential_text());
  const obdeect::propagation::RefractionTables tables(profile, 2000.);
  const auto along_x = bend_from(tables, 15000., 0.9, 0.);
  for (double azimuth : {0.3, 1.2, 2.5, -1.0}) {
    const auto r = bend_from(tables, 15000., 0.9, azimuth);
    CHECK_NEAR(r.x, along_x.x * std::cos(azimuth), 1e-9);
    CHECK_NEAR(r.y, along_x.x * std::sin(azimuth), 1e-9);
    CHECK_NEAR(r.t, along_x.t, 1e-9);
    CHECK_NEAR(std::hypot(r.cx, r.cy), std::hypot(along_x.cx, along_x.cy), 1e-14);
  }
  // |displacement| grows with zenith angle and emission height; the travel time with height.
  double last_x = 0.;
  for (double zenith = 0.1; zenith < 1.4; zenith += 0.1) {
    const double x = -bend_from(tables, 15000., zenith, 0.).x;
    CHECK(x > last_x);
    last_x = x;
  }
  double last_height_x = 0., last_time = 0.;
  for (double h = 3000.; h < 79000.; h += 4000.) {
    const auto r = bend_from(tables, h, 0.8, 0.);
    CHECK(-r.x > last_height_x);
    CHECK(r.t > last_time);
    last_height_x = -r.x;
    last_time = r.t;
  }
}

void test_zero_index_limit() {
  // A tenth of the refractivity: the displacement scales with it to first order.
  const obdeect::propagation::AtmosphereProfile full(exponential_text(1.0));
  const obdeect::propagation::AtmosphereProfile weak(exponential_text(1e-3));
  const obdeect::propagation::RefractionTables tables_full(full, 2000.), tables_weak(weak, 2000.);
  const double ratio =
      bend_from(tables_weak, 20000., 0.7, 0.).x / bend_from(tables_full, 20000., 0.7, 0.).x;
  CHECK_NEAR(ratio, 1e-3, 2e-6);
}

void test_extinction_gradient() {
  // tau is linear in log10(height) between the table heights and 0 at H2.
  const obdeect::propagation::TransmissionTable table(
      "# H2= 2.000, H1= 3.000 4.000 5.000\n400 0.1 0.2 0.3\n401 0.1 0.2 0.3\n", 2000.0);
  const double cz = 0.6, cx = 0.8;
  const double airmass = obdeect::propagation::airmass<double>(cx, 0.0);
  CHECK_NEAR(airmass, 1. / cz, 1e-14);
  CHECK_NEAR(table.transmission<double>(400., 4000., airmass), std::exp(-0.2 / cz), 1e-14);
  CHECK_NEAR(table.transmission<double>(400., 5000., airmass), std::exp(-0.3 / cz), 1e-14);
  // Between the nodes tau is linear in log10(h).
  const double f =
      (std::log10(4500.) - std::log10(4000.)) / (std::log10(5000.) - std::log10(4000.));
  CHECK_NEAR(table.optical_depth<double>(400., 4500.), 0.2 + 0.1 * f, 1e-14);
  // In-telescope path: telescope at 3 km with 100 m inside.
  const auto depths = table.telescope_depths(3000., 100.);
  CHECK_NEAR(depths.od_tel[400], 0.1, 1e-14);
  CHECK_NEAR(depths.od_focus[400], table.optical_depth<double>(400., 3100.) - 0.1, 1e-14);
  const double expected =
      std::exp(-(std::fabs(0.2 - depths.od_tel[400]) / cz + depths.od_focus[400]));
  CHECK_NEAR(
      table.transmission2<double>(400., 4000., airmass, depths.od_tel[400], depths.od_focus[400]),
      expected, 1e-14);
  // Transmission falls with airmass and with emission height.
  double previous = 1.;
  for (double h = 2100.; h < 5000.; h += 100.) {
    const double value = table.transmission<double>(400., h, 1.0);
    CHECK(value <= previous);
    previous = value;
  }
  CHECK(table.transmission<double>(400., 4000., 2.0) <
        table.transmission<double>(400., 4000., 1.0));
}

bool throws_invalid(const std::function<void()> &action) {
  try {
    action();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

void test_overflow_and_reflection() {
  // Absurd altitudes are rejected, not converted to a cell count.
  CHECK(throws_invalid([] {
    obdeect::propagation::AtmosphereProfile("0 1e-3 100 2e-4\n1e304 9e-4 90 2e-4\n1e305 8e-4 80 2e-4\n"
                             "1e306 7e-4 70 2e-4\n1e307 6e-4 60 2e-4\n");
  }));
  // Altitudes beyond 1000 km would make the fine grid unbounded.
  CHECK(throws_invalid([] {
    obdeect::propagation::AtmosphereProfile("0 1e-3 100 2e-4\n1 9e-4 90 2e-4\n2 8e-4 80 2e-4\n"
                             "3 7e-4 70 2e-4\n2000 6e-4 60 2e-4\n");
  }));
  // n-1 rising with height (valid table, unphysical): a ray from the layer above the denser
  // ground would need sin(theta_obs) > 1 and is not bent.
  std::ostringstream text;
  for (int km = 0; km <= 20; km += 2) {
    const double rho = kRho0 * std::exp(-km / kScaleHeightKm);
    const double rho_top = kRho0 * std::exp(-20. / kScaleHeightKm);
    text << km << ' ' << rho << ' ' << kScaleHeightKm * 1e5 * (rho - rho_top) << ' '
         << 1e-4 * (1 + km / 4.0) << '\n';
  }
  const obdeect::propagation::AtmosphereProfile rising(text.str());
  const obdeect::propagation::RefractionTables tables(rising, 0.);
  double cx = std::sin(1.55), cy = 0., x = 0., y = 0., t = 0.;
  CHECK(!obdeect::propagation::bend(tables, 15000., cx, cy, x, y, t));
  CHECK(cx == std::sin(1.55) && x == 0. && t == 0.);
}

void test_public_argument_validation() {
  const obdeect::propagation::AtmosphereProfile profile(exponential_text());
  const obdeect::propagation::RefractionTables tables(profile, 2000.);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK(throws_invalid([&] { (void)profile.index_moment_cm(nan, 1); }));
  CHECK(throws_invalid([&] { (void)profile.index_moment_cm(1e5, -1); }));
  CHECK(throws_invalid([&] { (void)profile.index_moment_cm(1e5, 7); }));
  CHECK(throws_invalid([&] { (void)tables.moment(1e5, -1); }));
  CHECK(throws_invalid([&] { (void)tables.moment(1e5, 7); }));
  CHECK(throws_invalid([&] { (void)tables.moment(nan, 1); }));
  CHECK(throws_invalid([&] { (void)obdeect::propagation::airmass<double>(1.1, 0.); }));
  CHECK(throws_invalid([&] { (void)obdeect::propagation::airmass<double>(nan, 0.); }));
  CHECK(obdeect::propagation::airmass<double>(1., 0.) == obdeect::propagation::unlimited_airmass);
  const obdeect::propagation::TransmissionTable open("# H2= 2.000, H1= 3.000\n400 0.5\n", 2000.);
  CHECK(throws_invalid([&] { (void)open.transmission<double>(1500., nan, 1.); }));
  CHECK(throws_invalid([&] { (void)obdeect::propagation::table_wavelength(400.5); }));
  CHECK(throws_invalid([&] { (void)obdeect::propagation::table_wavelength(-1.); }));
  CHECK(obdeect::propagation::table_wavelength(1e10) == -1);
  CHECK(obdeect::propagation::table_wavelength(400.) == 400);
}

void test_etadsn() {
  const obdeect::propagation::AtmosphereProfile profile(exponential_text());
  // Level more than 1.5 km above the observation level: n-1 / density of that level.
  const obdeect::propagation::RefractionTables low(profile, 2000.);
  CHECK_NEAR(low.etadsn(), kNu0 / kRho0, 1e-6);
  // No such level (observation level within 1.5 km of the top): the default.
  const obdeect::propagation::RefractionTables high(profile, 79000.);
  CHECK(high.etadsn() == obdeect::propagation::AtmosphereProfile::default_etadsn);
}

void test_table_header_comments() {
  const obdeect::propagation::TransmissionTable table(
      "# H2= 2.000, H1= 3.000 4.000 # heights\n400 0.1 0.2 # first row\n401 0.1 0.2\n", 2000.0);
  CHECK(table.last_wavelength_nm() == 401);
  // At most 100 heights including H2.
  std::string header = "# H2= 2.000, H1=", row = "400";
  for (int i = 1; i <= 100; ++i) { // 100 H1 plus H2
    header += " " + std::to_string(2.0 + 0.1 * i);
    row += " 0.5";
  }
  CHECK(throws_invalid([&] { obdeect::propagation::TransmissionTable(header + "\n" + row + "\n", 2000.0); }));
  CHECK_NEAR(table.optical_depth<double>(400., 4000.), 0.2, 1e-14);
}

} // namespace

int main() {
  return run_checks([] {
    test_profile_integrals();
    test_constant_index_is_straight();
    test_vertical_ray();
    test_against_direct_integration();
    test_real_profile_against_direct_integration();
    test_small_angle_limit();
    test_symmetry_and_monotonicity();
    test_zero_index_limit();
    test_extinction_gradient();
    test_public_argument_validation();
    test_overflow_and_reflection();
    test_etadsn();
    test_table_header_comments();
  });
}
