#include "check.hpp"

#include "dlufd/profile.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef DLUFD_TEST_DATA_DIR
#error "DLUFD_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

const std::string kProfilePath =
    std::string(DLUFD_TEST_DATA_DIR) + "/atmprof_ecmwf_north_winter_fixed.dat";

// Exponential atmosphere, 1 km levels to 10 km: rho = 1.2e-3 exp(-h/7 km) g/cm^3,
// thickness = integral from h to 10 km = 7e5 cm * (rho - rho(10 km)), n-1 = etadsn * rho
// (the relation atmo.c uses in the top layer).
std::string exponential_profile() {
  std::ostringstream text;
  text << "# exponential test profile\n";
  for (int km = 0; km <= 10; ++km) {
    const double rho = 1.2e-3 * std::exp(-km / 7.);
    const double rho_top = 1.2e-3 * std::exp(-10. / 7.);
    text << km << ' ' << rho << ' ' << 7e5 * (rho - rho_top) << ' '
         << dlufd::AtmosphereProfile::default_etadsn * rho << '\n';
  }
  return text.str();
}

bool rejects(const std::string &text) {
  try {
    dlufd::AtmosphereProfile profile(text);
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

} // namespace

int main() {
  return run_checks([] {
    // The exponential profile reproduces its analytic density, thickness and n-1 between levels.
    const dlufd::AtmosphereProfile exp_profile(exponential_profile());
    CHECK(exp_profile.levels() == 11);
    // Density and n-1 agree to 3e-6 (linear interpolation on the 3 m grid). Thickness is the
    // spline of log(thickness) through the truncated levels; it is within 0.5 % everywhere
    // and within 1e-3 below 7.5 km, 5e-3 above (coarse 1 km table).
    for (double km = 0.25; km < 9.9; km += 0.5) {
      const double h = 1000. * km;
      const double rho = 1.2e-3 * std::exp(-km / 7.);
      const double rho_top = 1.2e-3 * std::exp(-10. / 7.);
      CHECK_NEAR(exp_profile.density(h) / rho, 1., 1e-5);
      CHECK_NEAR(exp_profile.refractive_index_minus_one(h) /
                     (dlufd::AtmosphereProfile::default_etadsn * rho),
                 1., 1e-5);
      // The thickness spline is least accurate in the top 2 km of this coarse table.
      CHECK_NEAR(exp_profile.thickness(h) / (7e5 * (rho - rho_top)), 1., km < 7.5 ? 1e-3 : 5e-3);
    }

    // The ECMWF table (CORSIKA format, 55 levels, comment header skipped).
    const auto profile = dlufd::AtmosphereProfile::from_file(kProfilePath);
    CHECK(profile.levels() == 55);
    CHECK_NEAR(profile.bottom_m(), 0., 1e-9);
    CHECK_NEAR(profile.top_m(), 120000., 1e-6);

    // Table nodes: CORSIKA interpolates linearly on a 3 m grid of cubic-spline values, so
    // the nodes are reproduced to about 2e-5 relative (density), not to rounding.
    std::ifstream file(kProfilePath);
    std::string line;
    double worst_rho = 0., worst_thick = 0., worst_n1 = 0.;
    while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      std::istringstream row(line);
      double km = 0., rho = 0., thick = 0., n_1 = 0.;
      row >> km >> rho >> thick >> n_1;
      if (km >= 120.)
        continue; // The top level is outside the lookup range (CORSIKA returns 0 there).
      const double h = 1000. * km;
      worst_rho = std::fmax(worst_rho, std::fabs(profile.density(h) / rho - 1.));
      if (thick > 0.)
        worst_thick = std::fmax(worst_thick, std::fabs(profile.thickness(h) / thick - 1.));
      // Above 110 km the top layer uses n-1 = etadsn * rho (atmo.c); the 115 km node then
      // differs by 1e-3 relative, at n-1 ~ 8e-12.
      if (km < 110.)
        worst_n1 = std::fmax(worst_n1, std::fabs(profile.refractive_index_minus_one(h) / n_1 - 1.));
    }
    CHECK(worst_rho < 3e-5);
    CHECK(worst_thick < 1e-5);
    CHECK(worst_n1 < 3e-5);

    // Thickness decreases strictly with height over the whole table.
    double previous = profile.thickness(0.);
    for (double h = 10.; h <= 120000.; h += 10.) {
      const double current = profile.thickness(h);
      CHECK(current < previous);
      previous = current;
    }
    CHECK(profile.thickness(120000.) == 0.);
    CHECK(profile.density(130000.) == 0.);

    // Below the table the bottom values are returned.
    CHECK(profile.density(-5.) == profile.density(0.));

    // Input rejection.
    CHECK(rejects("0 1e-3 100\n1 1e-3 90\n2 1e-3 80\n3 1e-3 70\n4 1e-3 60\n")); // 3 columns
    CHECK(rejects(
        "0 1e-3 100 2e-4 x\n1 1e-3 90 2e-4\n2 1e-3 80 2e-4\n3 1e-3 70 2e-4\n4 1e-3 60 2e-4\n"));
    CHECK(rejects(exponential_profile().substr(0, 60))); // too few levels
    CHECK(rejects(
        "0 1e-3 100 2e-4\n1 1e-3 90 2e-4\n1 1e-3 80 2e-4\n3 1e-3 70 2e-4\n4 1e-3 60 2e-4\n"));
    CHECK(rejects(
        "0 1e-3 100 2e-4\n1 2e-3 90 2e-4\n2 1e-3 80 2e-4\n3 1e-3 70 2e-4\n4 5e-4 60 2e-4\n"));
    CHECK(
        rejects("0 1e-3 100 2e-4\n1 0 90 2e-4\n2 1e-3 80 2e-4\n3 1e-3 70 2e-4\n4 1e-3 60 2e-4\n"));
    CHECK(
        rejects("0 1e-3 100 0\n1 1e-3 90 2e-4\n2 1e-3 80 2e-4\n3 1e-3 70 2e-4\n4 1e-3 60 2e-4\n"));
    CHECK(rejects(
        "0 1e-3 100 2e-4\n1 9e-4 90 2e-4\n2 8e-4 80 2e-4\n3 7e-4 70 2e-4\n4 6e-4 9000 2e-4\n"));
  });
}
