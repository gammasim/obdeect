#include "check.hpp"

#include "obdeect/propagation/profile.hpp"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

const std::string kProfilePath =
    std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/atmprof_ecmwf_north_winter_fixed.dat";

// Exponential atmosphere, 1 km levels to 10 km: rho = 1.2e-3 exp(-h/7 km) g/cm^3,
// thickness = integral from h to 10 km = 7e5 cm * (rho - rho(10 km)), n-1 = etadsn * rho
// (the relation of the top layer).
std::string exponential_row(int km) {
  const double rho = 1.2e-3 * std::exp(-km / 7.);
  const double rho_top = 1.2e-3 * std::exp(-10. / 7.);
  std::ostringstream row;
  row << km << ' ' << rho << ' ' << 7e5 * (rho - rho_top) << ' '
      << obdeect::propagation::AtmosphereProfile::default_etadsn * rho;
  return row.str();
}

// Exponential profile. `trailing` is appended to every row; row `replaced_km` (if >= 0) is
// written as `replacement` instead.
std::string exponential_profile(const std::string &trailing = "", int replaced_km = -1,
                                const std::string &replacement = "") {
  std::ostringstream text;
  text << "# exponential test profile\n";
  for (int km = 0; km <= 10; ++km) {
    if (km == replaced_km)
      text << replacement << '\n';
    else
      text << exponential_row(km) << trailing << '\n';
  }
  return text.str();
}

// exponential_row(km) with its column `column` (0 = altitude .. 3 = n-1) replaced by `token`.
std::string with_column(int km, std::size_t column, const std::string &token) {
  std::istringstream row(exponential_row(km));
  std::vector<std::string> columns;
  for (std::string column_text; row >> column_text;)
    columns.push_back(column_text);
  columns.at(column) = token;
  std::string joined;
  for (const auto &column_text : columns)
    joined += (joined.empty() ? "" : " ") + column_text;
  return joined;
}

// The error message of a rejected profile, or an empty string if it is accepted.
std::string parse_error(const std::string &text) {
  try {
    obdeect::propagation::AtmosphereProfile profile(text);
  } catch (const std::invalid_argument &error) {
    return error.what();
  }
  return "";
}

bool rejects(const std::string &text) { return !parse_error(text).empty(); }

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
    // The exponential profile reproduces its analytic density, thickness and n-1 between levels.
    const obdeect::propagation::AtmosphereProfile exp_profile(exponential_profile());
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
                     (obdeect::propagation::AtmosphereProfile::default_etadsn * rho),
                 1., 1e-5);
      // The thickness spline is least accurate in the top 2 km of this coarse table.
      CHECK_NEAR(exp_profile.thickness(h) / (7e5 * (rho - rho_top)), 1., km < 7.5 ? 1e-3 : 5e-3);
    }

    // The ECMWF table (atmprof format, 55 levels, comment header skipped).
    const auto profile = obdeect::propagation::AtmosphereProfile::from_file(kProfilePath);
    CHECK(profile.levels() == 55);
    CHECK_NEAR(profile.bottom_m(), 0., 1e-9);
    CHECK_NEAR(profile.top_m(), 120000., 1e-6);

    // Table nodes: lookups interpolate linearly on a 3 m grid of cubic-spline values, so
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
        continue; // The top level is outside the lookup range (the lookup returns 0 there).
      const double h = 1000. * km;
      worst_rho = std::fmax(worst_rho, std::fabs(profile.density(h) / rho - 1.));
      if (thick > 0.)
        worst_thick = std::fmax(worst_thick, std::fabs(profile.thickness(h) / thick - 1.));
      // Above 110 km the top layer uses n-1 = etadsn * rho; the 115 km node then
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

    // Non-finite heights are rejected by every public height lookup: NaN fails all boundary
    // comparisons and would reach an integer conversion, infinities would be clamped.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (double bad : {nan, inf, -inf}) {
      CHECK(throws_invalid([&] { profile.density(bad); }));
      CHECK(throws_invalid([&] { profile.thickness(bad); }));
      CHECK(throws_invalid([&] { profile.refractive_index_minus_one(bad); }));
      CHECK(throws_invalid([&] { profile.density_cm(bad); }));
      CHECK(throws_invalid([&] { profile.thickness_cm(bad); }));
      CHECK(throws_invalid([&] { profile.index_minus_one_cm(bad); }));
      CHECK(throws_invalid([&] { profile.index_integral_cm(bad); }));
    }
    // A finite height whose conversion to cm overflows is non-finite too.
    CHECK(throws_invalid([&] { profile.density(std::numeric_limits<double>::max()); }));
    CHECK(throws_invalid([&] { profile.thickness(-std::numeric_limits<double>::max()); }));
    // Finite heights outside the table keep their defined values (not an error).
    CHECK(profile.density(1e9) == 0.);

    // Input rejection.
    CHECK(rejects("0 1e-3 100\n1 1e-3 90\n2 1e-3 80\n3 1e-3 70\n4 1e-3 60\n")); // 3 columns
    // A malformed token in any of the four parsed columns is rejected, and the error names the
    // line (row km = 5 is line 7: one comment line, then rows from km = 0).
    for (std::size_t column = 0; column < 4; ++column) {
      for (const std::string token : {"x", "2e-4x"}) {
        const std::string error =
            parse_error(exponential_profile("", 5, with_column(5, column, token)));
        CHECK(error.find("not a finite number") != std::string::npos);
        CHECK(error.find("line 7") != std::string::npos);
      }
    }
    // Trailing columns beyond the four are ignored: they are accepted
    // and give the same table as the four-column text, whatever they contain.
    const obdeect::propagation::AtmosphereProfile four_columns(exponential_profile());
    for (const std::string trailing : {" x", " 1e-3 2 abc", " 0.5"}) {
      const obdeect::propagation::AtmosphereProfile extra(exponential_profile(trailing));
      CHECK(extra.levels() == four_columns.levels());
      for (double h = 0.; h <= 10000.; h += 333.) {
        CHECK(extra.density(h) == four_columns.density(h));
        CHECK(extra.thickness(h) == four_columns.thickness(h));
        CHECK(extra.refractive_index_minus_one(h) == four_columns.refractive_index_minus_one(h));
      }
    }
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
