#include "check.hpp"
#include "sky_stream_check.hpp"

#include "dlamb/programs/xyzls.hpp"

#include "../tools/telpos_file.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

using namespace dlamb;
using namespace dlamb::programs;
using gstrickd::Column;

namespace {

constexpr double deg = pi / 180.0;

// Default: 1e7 photons into 4 pi, source 1000 m above the origin, telescope radius 15 m.
// The original emits only into the cone that contains the 1.02 R sphere:
// frac = 0.5 (1 - cos(asin(1.02 R / d))).
void default_count_matches_original() {
  const XyzlsOptions o;
  const auto plan = plan_xyzls(o);
  CHECK(plan.cases.size() == 1);
  const auto &em = plan.cases[0][0];
  CHECK(std::holds_alternative<SphereTarget>(em.source.target));
  const double cos_max = std::sqrt(1.0 - std::pow(1.02 * 15.0 / 1000.0, 2));
  const double expected = 1e7 * 0.5 * (1.0 - cos_max);
  const int n = 4000;
  double sum = 0.0;
  for (std::uint32_t e = 1; e <= static_cast<std::uint32_t>(n); ++e)
    sum += static_cast<double>(em.count(e));
  // Each count is floor(mean) plus a Bernoulli term: variance <= 1/4.
  CHECK_NEAR(sum / n, expected, 4.0 * 0.5 / std::sqrt(static_cast<double>(n)));
}

// Default run: stream layout, frames, cone of directions, delta pulse, 400 nm.
void default_stream() {
  XyzlsOptions o;
  o.events = 2;
  const std::string bytes = write_to_string([&](std::ostream &out) { write_xyzls(out, o); });
  const auto s = read_stream(bytes);
  CHECK(s.ended);
  CHECK(!s.types.empty() && s.types.front() == 1 && s.types.back() == 8);
  CHECK(s.events == 2);
  CHECK(std::get<std::int64_t>(s.run.at("n_telescopes")) == 1);
  CHECK(std::get<std::int64_t>(s.run.at("refraction_applied")) == 0);

  const auto plan = plan_xyzls(o);
  const std::uint64_t expected = plan.cases[0][0].count(1) + plan.cases[0][0].count(2);
  CHECK(s.photons.at(1) == expected);
  const double cos_max = std::sqrt(1.0 - std::pow(1.02 * 15.0 / 1000.0, 2));
  for (const auto &c : s.chunks) {
    CHECK(c.frame == gstrickd::Frame::sky);
    for (std::size_t i = 0; i < c.rows(); ++i) {
      CHECK_NEAR(c[Column::w][i], 1.0, 1e-6);
      CHECK_NEAR(c[Column::wavelength][i], 400.0, 1e-4);
      CHECK_NEAR(c[Column::t][i], 0.0, 1e-6);
      CHECK_NEAR(c[Column::z][i], 1000.0, 1e-3);
      CHECK(c[Column::cz][i] <= -cos_max + 1e-6);
    }
  }
}

// Non-default: cone of 30 deg (no fiducial cut for fixed angles), 2 photons per bunch,
// exponential pulse, 300 nm, two events. Counts, directions, weights and moments are exact.
void non_default_moments() {
  XyzlsOptions o;
  o.photons = {2e5};
  o.events = 2;
  o.bunch_size = 2.0;
  o.z_cm = 20000.0;
  o.angular = FixedAngle{1.0 - std::cos(30.0 * deg)};
  o.pulse = Exponential{0.0, 2.0};
  o.spectrum = Fixed{300.0};
  const auto s = read_stream(write_to_string([&](std::ostream &out) { write_xyzls(out, o); }));
  CHECK(s.photons.at(1) == 2 * 100000U);
  double t1 = 0.0, t2 = 0.0;
  std::uint64_t n = 0;
  for (const auto &c : s.chunks) {
    for (std::size_t i = 0; i < c.rows(); ++i) {
      CHECK_NEAR(c[Column::w][i], 2.0, 1e-6);
      CHECK_NEAR(c[Column::wavelength][i], 300.0, 1e-4);
      CHECK_NEAR(c[Column::cz][i], -std::cos(30.0 * deg), 1e-5);
      const double t = c[Column::t][i];
      t1 += t;
      t2 += t * t;
      ++n;
    }
  }
  const double mean = t1 / static_cast<double>(n);
  const double var = (t2 - n * mean * mean) / static_cast<double>(n - 1);
  // Exponential with decay 2 ns and mean 0: variance 4, kurtosis 9.
  CHECK_NEAR(mean, 0.0, 4.0 * 2.0 / std::sqrt(static_cast<double>(n)));
  CHECK_NEAR(var, 4.0, 4.0 * 4.0 * std::sqrt(8.0 / static_cast<double>(n)));
}

// Non-isotropic profiles beyond the fiducial sphere are rejected, not silently changed.
void rejects_unported_profile() {
  XyzlsOptions o;
  o.angular = Lambertian{};
  bool threw = false;
  try {
    (void)plan_xyzls(o);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
}

// Emission altitude: zem = observation level + z for every sky photon (format.md, zem).
void emission_altitude() {
  XyzlsOptions o;
  o.altitude_m = 1234.5;
  o.z_cm = 20000.0;
  o.photons = {2e4};
  const auto s = read_stream(write_to_string([&](std::ostream &out) { write_xyzls(out, o); }));
  CHECK(!s.chunks.empty());
  for (const auto &c : s.chunks) {
    CHECK((c.mask() & gstrickd::column_bit(Column::zem)) != 0U);
    for (std::size_t i = 0; i < c.rows(); ++i)
      CHECK_NEAR(c[Column::zem][i] - c[Column::z][i], 1234.5, 1e-2);
  }
}

// Negative photon counts are rejected, not clamped to zero.
void rejects_negative_photons() {
  XyzlsOptions o;
  o.photons = {1e6, -1.0};
  bool threw = false;
  try {
    (void)plan_xyzls(o);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
}

// Telescope lists of --telpos-file / --corsika-telpos-file (see cpp/tools/telpos_file.hpp).
std::vector<TelposCm> telpos(const std::string &text, bool corsika) {
  std::istringstream in(text);
  return cli::read_telpos(in, "--telpos-file", corsika);
}

bool telpos_rejects(const std::string &text, bool corsika) {
  try {
    (void)telpos(text, corsika);
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

// Valid files: blank lines and comments are skipped; each line gives one telescope.
void telpos_valid() {
  const auto list = telpos("# two telescopes\n\n0 0 0 1500\n 100 -50 0 1200 # comment\n", false);
  CHECK(list.size() == 2);
  CHECK(list[1].x == 100.0 && list[1].y == -50.0 && list[1].z == 0.0 && list[1].r == 1200.0);
  const auto corsika =
      telpos("SHOWER 1\nTELESCOPE 0 0 0 1500\nOBSLEV 2150\nTELESCOPE 10 20 0 500\n", true);
  CHECK(corsika.size() == 2);
  CHECK(corsika[1].x == 10.0 && corsika[1].y == 20.0 && corsika[1].r == 500.0);
}

// Empty files, or files without telescope records, are errors (no silent default telescope).
void telpos_rejects_empty() {
  CHECK(telpos_rejects("", false));
  CHECK(telpos_rejects("# nothing here\n\n", false));
  CHECK(telpos_rejects("OBSLEV 2150\n", true));
  CHECK(telpos_rejects("", true));
}

// Malformed lines are errors, even when other lines are valid.
void telpos_rejects_malformed() {
  CHECK(telpos_rejects("0 0 0 1500\n1 2\n", false));
  CHECK(telpos_rejects("0 0 0\n", false));
  CHECK(telpos_rejects("0 0 0 1500 9\n", false));
  CHECK(telpos_rejects("0 0 x 1500\n", false));
  CHECK(telpos_rejects("0 0 0 0\n", false));
  CHECK(telpos_rejects("0 0 0 -5\n", false));
  CHECK(telpos_rejects("0 0 0 nan\n", false));
  CHECK(telpos_rejects("0 0 0 inf\n", false));
  CHECK(telpos_rejects("TELESCOPE 0 0 0\n", true));
  CHECK(telpos_rejects("TELESCOPE 0 0 0 1500 9\n", true));
}

} // namespace

int main() {
  return run_checks([] {
    default_count_matches_original();
    default_stream();
    non_default_moments();
    rejects_unported_profile();
    emission_altitude();
    rejects_negative_photons();
    telpos_valid();
    telpos_rejects_empty();
    telpos_rejects_malformed();
  });
}
