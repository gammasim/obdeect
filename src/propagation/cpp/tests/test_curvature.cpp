// Spherical-atmosphere extinction: slant optical depth, its flat-atmosphere limit, a brute-force
// midpoint integration of the same piecewise integrand, and the per-row kernels.
#include "check.hpp"

#include "obdeect/propagation/extinction.hpp"
#include "obdeect/propagation/pipeline.hpp"
#include "obdeect/propagation/thinning.hpp"
#include "obdeect/propagation/transmission.hpp"
#include "obdeect/propagation/wavelength.hpp"

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using obdeect::stream::Column;

constexpr double kEarthRadius = obdeect::stream::earth_radius_m;
constexpr double kSiteAltitude = 2000.0; // H2 of the fixture and of the synthetic table
constexpr double kTelescope = 2030.0;    // telescope altitude [m]
constexpr std::uint64_t kSeed = 20261114;
const double kPi = obdeect::stream::pi;

const std::string kFixture =
    std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/atm_trans_2000_1_10_0_0_2000.dat";

template <class Function> bool throws_invalid(Function function) {
  try {
    function();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

double relative(double value, double reference) {
  return std::abs(value - reference) / std::abs(reference);
}

// Synthetic table: 80 heights from 2.1 km to 100 km (geometric), one row at 400 nm with
// tau(h) = 1.5 (1 - exp(-(h - H2) / 6 km)) sampled at the heights.
struct SyntheticTable {
  std::vector<double> height_m; // including H2 first
  std::vector<double> tau;      // tau at height_m, tau[0] = 0
};

SyntheticTable synthetic_table() {
  SyntheticTable table;
  table.height_m.push_back(kSiteAltitude);
  table.tau.push_back(0.0);
  for (int k = 0; k < 80; ++k) {
    const double h_km = 2.0 * std::pow(50.0, (k + 1) / 80.0);
    table.height_m.push_back(1000.0 * h_km);
    table.tau.push_back(1.5 * (1.0 - std::exp(-(h_km - 2.0) / 6.0)));
  }
  return table;
}

std::string table_text(const SyntheticTable &table) {
  std::ostringstream text;
  text.precision(17);
  text << "# H2= 2.000, H1=";
  for (std::size_t k = 1; k < table.height_m.size(); ++k) {
    text << ' ' << table.height_m[k] / 1000.0;
  }
  text << "\n400";
  for (std::size_t k = 1; k < table.tau.size(); ++k) {
    text << ' ' << table.tau[k];
  }
  text << '\n';
  return text.str();
}

// Brute-force slant depth of the piecewise integrand: midpoint rule with n points per table
// segment, written independently of the kernel (no Gauss nodes, direct cos theta(h) formula).
// The first segment at h_t uses h = h_t + u^2 when cos_t < 0.2 (the same singularity handling).
double brute_slant(const SyntheticTable &table, double h_t, double zem, double cos_t, double R,
                   int n) {
  const double ln10 = std::log(10.0);
  const double p = (R + h_t) * std::sqrt(1.0 - cos_t * cos_t);
  auto cos_theta = [&](double h) {
    const double r = R + h;
    return std::sqrt(1.0 - p * p / (r * r));
  };
  double total = 0.0;
  for (std::size_t i = 0; i + 1 < table.height_m.size(); ++i) {
    const double lo = std::max(table.height_m[i], h_t);
    const double hi = std::min(table.height_m[i + 1], zem);
    if (!(hi > lo)) {
      continue;
    }
    const double slope = (table.tau[i + 1] - table.tau[i]) /
                         (std::log10(table.height_m[i + 1]) - std::log10(table.height_m[i]));
    auto integrand = [&](double h) { return slope / (h * ln10) / cos_theta(h); };
    if (cos_t < 0.2 && lo == h_t) {
      const double top = std::sqrt(hi - h_t);
      const double du = top / n;
      for (int k = 0; k < n; ++k) {
        const double u = (k + 0.5) * du;
        total += integrand(h_t + u * u) * 2.0 * u * du;
      }
    } else {
      const double dh = (hi - lo) / n;
      for (int k = 0; k < n; ++k) {
        total += integrand(lo + (k + 0.5) * dh) * dh;
      }
    }
  }
  return total;
}

// Plane-parallel slant depth of the table: |tau(zem) - tau(h_t)| / cos_t.
double flat_slant(const obdeect::propagation::TransmissionTable &table, double h_t, double zem,
                  double cos_t) {
  return std::abs(table.optical_depth<double>(400.0, zem) -
                  table.optical_depth<double>(400.0, h_t)) /
         cos_t;
}

void test_flat_limit_fixture() {
  const auto table = obdeect::propagation::TransmissionTable::from_file(kFixture, kSiteAltitude);
  for (const double cos_t : {1.0, 0.7, 0.3}) {
    for (const double zem : {2500.0, 5000.0, 20000.0}) {
      const double curved = table.slant_optical_depth<double>(400.0, kTelescope, zem, cos_t, 1e12);
      CHECK(relative(curved, flat_slant(table, kTelescope, zem, cos_t)) <= 1e-6);
    }
  }
}

void test_flat_limit_and_brute_force_synthetic() {
  const auto synthetic = synthetic_table();
  const obdeect::propagation::TransmissionTable table(table_text(synthetic), kSiteAltitude);
  // Flat-atmosphere limit on the synthetic grid.
  for (const double cos_t : {1.0, 0.7, 0.3}) {
    for (const double zem : {3000.0, 30000.0, 100000.0}) {
      const double curved = table.slant_optical_depth<double>(400.0, kTelescope, zem, cos_t, 1e12);
      CHECK(relative(curved, flat_slant(table, kTelescope, zem, cos_t)) <= 1e-6);
    }
  }
  // Spherical geometry against the brute-force midpoint integral.
  for (const double degrees : {0.0, 45.0, 70.0, 85.0, 89.9}) {
    const double cos_t = std::cos(degrees * kPi / 180.0);
    for (const double zem : {30000.0, 60000.0}) {
      const double kernel =
          table.slant_optical_depth<double>(400.0, kTelescope, zem, cos_t, kEarthRadius);
      const double brute = brute_slant(synthetic, kTelescope, zem, cos_t, kEarthRadius, 100000);
      CHECK(brute > 0.0);
      CHECK(relative(kernel, brute) <= 1e-5);
    }
  }
}

void test_curvature_at_large_angles() {
  const auto table = obdeect::propagation::TransmissionTable::from_file(kFixture, kSiteAltitude);
  const double zem = 30000.0;
  const double cos80 = std::cos(80.0 * kPi / 180.0);
  const double curved80 =
      table.slant_optical_depth<double>(400.0, kTelescope, zem, cos80, kEarthRadius);
  CHECK(curved80 > 0.0);
  CHECK(curved80 < flat_slant(table, kTelescope, zem, cos80));
  // The curvature reduction at 80 degrees is about 3 % (ratio 0.972 for this table).
  CHECK(curved80 < 0.99 * flat_slant(table, kTelescope, zem, cos80));

  const double cos8999 = std::cos(89.99 * kPi / 180.0);
  const double curved8999 =
      table.slant_optical_depth<double>(400.0, kTelescope, zem, cos8999, kEarthRadius);
  const double plane8999 = flat_slant(table, kTelescope, zem, cos8999);
  CHECK(std::isfinite(curved8999) && curved8999 > 0.0);
  CHECK(plane8999 > 10.0 * curved8999);

  // Grazing photons: finite, and the slant depth decreases with the cosine.
  const double grazing =
      table.slant_optical_depth<double>(400.0, kTelescope, zem, 0.0, kEarthRadius);
  CHECK(std::isfinite(grazing));
  CHECK(grazing >= curved8999);
  double previous = grazing;
  for (const double cos_t : {0.001, 0.01, 0.1, 0.3, 0.7, 1.0}) {
    const double value =
        table.slant_optical_depth<double>(400.0, kTelescope, zem, cos_t, kEarthRadius);
    CHECK(std::isfinite(value));
    CHECK(value <= previous);
    previous = value;
  }
  // Vertical photons: curvature is irrelevant.
  CHECK(relative(table.slant_optical_depth<double>(400.0, kTelescope, zem, 1.0, kEarthRadius),
                 flat_slant(table, kTelescope, zem, 1.0)) <= 1e-9);
}

void test_monotone_zero_and_validation() {
  const auto table = obdeect::propagation::TransmissionTable::from_file(kFixture, kSiteAltitude);
  const double cos_t = std::cos(70.0 * kPi / 180.0);
  double previous = 0.0;
  for (double zem = kTelescope - 500.0; zem <= 60000.0; zem += 250.0) {
    const double value =
        table.slant_optical_depth<double>(400.0, kTelescope, zem, cos_t, kEarthRadius);
    CHECK(value >= previous);
    CHECK(value >= 0.0);
    previous = value;
  }
  CHECK(table.slant_optical_depth<double>(400.0, kTelescope, kTelescope, cos_t, kEarthRadius) ==
        0.0);
  CHECK(table.slant_optical_depth<double>(400.0, kTelescope, kTelescope - 1000.0, cos_t,
                                          kEarthRadius) == 0.0);

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK(throws_invalid(
      [&] { (void)table.slant_optical_depth<double>(400.0, kTelescope, 30000.0, cos_t, 0.0); }));
  CHECK(throws_invalid([&] {
    (void)table.slant_optical_depth<double>(400.0, kTelescope, 30000.0, cos_t, -kEarthRadius);
  }));
  CHECK(throws_invalid(
      [&] { (void)table.slant_optical_depth<double>(400.0, kTelescope, 30000.0, cos_t, nan); }));
  CHECK(throws_invalid(
      [&] { (void)table.slant_optical_depth<double>(400.0, kTelescope, 30000.0, cos_t, inf); }));
  CHECK(throws_invalid([&] {
    (void)table.slant_optical_depth<double>(400.0, kTelescope, 30000.0, -0.1, kEarthRadius);
  }));
  CHECK(throws_invalid([&] {
    (void)table.slant_optical_depth<double>(400.0, kTelescope, 30000.0, 1.0 + 1e-6, kEarthRadius);
  }));
  CHECK(throws_invalid([&] {
    (void)table.slant_optical_depth<double>(400.0, kTelescope, 30000.0, nan, kEarthRadius);
  }));
  CHECK(throws_invalid(
      [&] { (void)table.slant_optical_depth<double>(400.0, nan, 30000.0, cos_t, kEarthRadius); }));
  CHECK(throws_invalid([&] {
    (void)table.slant_optical_depth<double>(400.0, kTelescope, nan, cos_t, kEarthRadius);
  }));
  CHECK(throws_invalid([&] {
    (void)table.slant_optical_depth<double>(nan, kTelescope, 30000.0, cos_t, kEarthRadius);
  }));
  // Rounding-size excess over 1 is accepted.
  CHECK(!throws_invalid([&] {
    (void)table.slant_optical_depth<double>(400.0, kTelescope, 30000.0, 1.0 + 1e-12, kEarthRadius);
  }));
  CHECK(throws_invalid([&] {
    (void)table.transmission_curved<double>(400.0, 30000.0, cos_t, kTelescope, -1.0, kEarthRadius);
  }));
  CHECK(throws_invalid([&] {
    (void)table.transmission_curved<double>(400.0, 30000.0, cos_t, kTelescope, nan, kEarthRadius);
  }));
}

void test_transmission_curved() {
  const auto table = obdeect::propagation::TransmissionTable::from_file(kFixture, kSiteAltitude);
  const double cos_t = 0.5;
  const double zem = 20000.0;
  const double focus = 0.02;
  const double slant =
      table.slant_optical_depth<double>(400.0, kTelescope, zem, cos_t, kEarthRadius);
  CHECK_NEAR(table.transmission_curved<double>(400.0, zem, cos_t, kTelescope, focus, kEarthRadius),
             std::exp(-(slant + focus)), 1e-14);
  // Below the observation level, and for a telescope depth above max_telescope_depth: zero.
  CHECK(table.transmission_curved<double>(400.0, 1500.0, cos_t, kTelescope, focus, kEarthRadius) ==
        0.0);
  // Zem below the telescope: no slant depth, transmission = exp(-od_focus).
  CHECK_NEAR(
      table.transmission_curved<double>(400.0, 2010.0, cos_t, kTelescope, focus, kEarthRadius),
      std::exp(-focus), 1e-14);
  // A total depth above max_total_depth gives zero: tau rises to 300 at 3 km.
  const obdeect::propagation::TransmissionTable opaque("# H2= 2.000, H1= 3.000\n400 300\n",
                                                       kSiteAltitude);
  CHECK(opaque.slant_optical_depth<double>(400.0, kTelescope, 3000.0, 0.5, kEarthRadius) > 100.0);
  CHECK(opaque.transmission_curved<double>(400.0, 3000.0, 0.5, kTelescope, focus, kEarthRadius) ==
        0.0);
  // Float evaluation agrees with double to single precision.
  const float as_float =
      table.transmission_curved<float>(400.0F, 20000.0F, 0.5F, 2030.0F, 0.02F, kEarthRadius);
  CHECK_NEAR(as_float, std::exp(-(slant + focus)), 1e-6);
}

// Chunk of rows generated from their global index g (direction, emission altitude, wavelength),
// so that chunks split or reversed in order carry the same photons.
obdeect::stream::PhotonChunk make_chunk(std::size_t first, std::size_t count, bool reversed,
                                        std::size_t total) {
  obdeect::stream::PhotonChunk chunk;
  chunk.resize(count, obdeect::stream::column_bit(Column::cx) |
                          obdeect::stream::column_bit(Column::cy) |
                          obdeect::stream::column_bit(Column::zem) |
                          obdeect::stream::column_bit(Column::wavelength) |
                          obdeect::stream::column_bit(Column::w));
  for (std::size_t r = 0; r < count; ++r) {
    const std::size_t g = reversed ? total - 1 - (first + r) : first + r;
    const double zenith = 1.5 * static_cast<double>(g) / static_cast<double>(total);
    const double azimuth = 0.7 * static_cast<double>(g);
    chunk[Column::cx][r] = static_cast<float>(std::sin(zenith) * std::cos(azimuth));
    chunk[Column::cy][r] = static_cast<float>(std::sin(zenith) * std::sin(azimuth));
    chunk[Column::zem][r] = static_cast<float>(2500.0 + 300.0 * static_cast<double>(g % 120));
    const double wavelengths[] = {300.0, 350.0, 400.0, 550.0, 700.0, 1000.0};
    chunk[Column::wavelength][r] = static_cast<float>(wavelengths[g % 6]);
    chunk[Column::w][r] = 1.0F;
  }
  return chunk;
}

void test_acceptance_kernel() {
  const auto table = obdeect::propagation::TransmissionTable::from_file(kFixture, kSiteAltitude);
  const auto depths = table.telescope_depths(kTelescope, 28.0);
  const obdeect::propagation::CurvedPath path{kEarthRadius, kTelescope, 0.0, 0.0, 1.0};
  const std::size_t n = 240;

  // Equals transmission_curved row by row; the cosine is -(d . n) = sqrt(1 - cx^2 - cy^2).
  const auto chunk = make_chunk(0, n, false, n);
  std::vector<double> p(n);
  obdeect::propagation::atmospheric_acceptance_curved<double>(chunk, table, depths, path, p);
  for (std::size_t row = 0; row < n; ++row) {
    const double cx = chunk[Column::cx][row], cy = chunk[Column::cy][row];
    const double cos_t = std::sqrt(1.0 - cx * cx - cy * cy);
    const double wavelength = chunk[Column::wavelength][row];
    const auto iwl = static_cast<std::size_t>(wavelength);
    if (wavelength >= 1000.0F) {
      CHECK(p[row] == 0.0);
      continue;
    }
    const double expected = table.transmission_curved<double>(
        wavelength, chunk[Column::zem][row], cos_t, kTelescope, depths.od_focus[iwl], kEarthRadius);
    CHECK_NEAR(p[row], expected, 1e-14);
  }

  // Row order: reversed rows carry the same photons.
  const auto reversed = make_chunk(0, n, true, n);
  std::vector<double> q(n);
  obdeect::propagation::atmospheric_acceptance_curved<double>(reversed, table, depths, path, q);
  for (std::size_t row = 0; row < n; ++row) {
    CHECK_NEAR(q[row], p[n - 1 - row], 0.0);
  }

  // Chunking: halves give the same values as the whole chunk.
  std::vector<double> first_half(n / 2), second_half(n - n / 2);
  const auto half_a = make_chunk(0, n / 2, false, n);
  const auto half_b = make_chunk(n / 2, n - n / 2, false, n);
  obdeect::propagation::atmospheric_acceptance_curved<double>(half_a, table, depths, path,
                                                              first_half);
  obdeect::propagation::atmospheric_acceptance_curved<double>(half_b, table, depths, path,
                                                              second_half);
  for (std::size_t row = 0; row < n / 2; ++row) {
    CHECK_NEAR(first_half[row], p[row], 0.0);
  }
  for (std::size_t row = n / 2; row < n; ++row) {
    CHECK_NEAR(second_half[row - n / 2], p[row], 0.0);
  }

  // Tilted vertical: cos_zenith = sqrt(1 - cx^2 - cy^2) cos a - cx sin a.
  const double a = 0.3;
  const obdeect::propagation::CurvedPath tilted{kEarthRadius, kTelescope, std::sin(a), 0.0,
                                                std::cos(a)};
  std::vector<double> t(n);
  obdeect::propagation::atmospheric_acceptance_curved<double>(chunk, table, depths, tilted, t);
  std::size_t changed = 0;
  for (std::size_t row = 0; row < n; ++row) {
    const double cx = chunk[Column::cx][row], cy = chunk[Column::cy][row];
    const double wavelength = chunk[Column::wavelength][row];
    if (wavelength >= 1000.0F) {
      continue;
    }
    const double cos_t = std::sqrt(1.0 - cx * cx - cy * cy) * std::cos(a) - cx * std::sin(a);
    if (cos_t < 0.0) {
      CHECK(t[row] == 0.0); // moves away from the ground
      continue;
    }
    const auto iwl = static_cast<std::size_t>(wavelength);
    const double expected = table.transmission_curved<double>(
        wavelength, chunk[Column::zem][row], cos_t, kTelescope, depths.od_focus[iwl], kEarthRadius);
    CHECK_NEAR(t[row], expected, 1e-14);
    if (std::abs(t[row] - p[row]) > 1e-9) {
      ++changed;
    }
  }
  CHECK(changed > n / 4);

  // A photon moving away from the ground (cos_zenith < 0 for n = (1, 0, 0), cx = 0.6) gets p = 0.
  const obdeect::propagation::CurvedPath up{kEarthRadius, kTelescope, 1.0, 0.0, 0.0};
  auto away = make_chunk(0, 1, false, 1);
  away[Column::cx][0] = 0.6F;
  away[Column::cy][0] = 0.0F;
  away[Column::wavelength][0] = 400.0F;
  std::vector<double> none(1, 1.0);
  obdeect::propagation::atmospheric_acceptance_curved<double>(away, table, depths, up, none);
  CHECK(none[0] == 0.0);

  // Invalid paths and sizes.
  std::vector<double> wrong(n + 1);
  CHECK(throws_invalid([&] {
    obdeect::propagation::atmospheric_acceptance_curved<double>(chunk, table, depths, path, wrong);
  }));
  const obdeect::propagation::CurvedPath not_unit{kEarthRadius, kTelescope, 1.0, 0.0, 0.1};
  CHECK(throws_invalid([&] { obdeect::propagation::validate(not_unit); }));
  CHECK(throws_invalid([&] {
    obdeect::propagation::atmospheric_acceptance_curved<double>(chunk, table, depths, not_unit, p);
  }));
  const obdeect::propagation::CurvedPath no_radius{0.0, kTelescope, 0.0, 0.0, 1.0};
  CHECK(throws_invalid([&] { obdeect::propagation::validate(no_radius); }));
  obdeect::propagation::validate(path);
  obdeect::propagation::validate(tilted);
}

// Settings for apply_extinction_chunk with the curved atmosphere, weighted output (p column).
void test_apply_extinction_curved() {
  const auto table = obdeect::propagation::TransmissionTable::from_file(kFixture, kSiteAltitude);
  const auto depths = table.telescope_depths(kTelescope, 28.0);
  const obdeect::propagation::Envelope envelope({200.0, 800.0}, {1.0, 1.0});
  const obdeect::propagation::CurvedPath path{kEarthRadius, kTelescope, 0.0, 0.0, 1.0};

  auto make_photons = [] {
    obdeect::stream::PhotonChunk chunk;
    chunk.resize(
        60, obdeect::stream::column_bit(Column::x) | obdeect::stream::column_bit(Column::y) |
                obdeect::stream::column_bit(Column::cx) | obdeect::stream::column_bit(Column::cy) |
                obdeect::stream::column_bit(Column::t) | obdeect::stream::column_bit(Column::w) |
                obdeect::stream::column_bit(Column::zem) |
                obdeect::stream::column_bit(Column::wavelength));
    for (std::size_t r = 0; r < chunk.rows(); ++r) {
      const double zenith = 1.4 * static_cast<double>(r) / 60.0;
      chunk[Column::cx][r] = static_cast<float>(std::sin(zenith));
      chunk[Column::cy][r] = 0.0F;
      chunk[Column::zem][r] = 25000.0F;
      chunk[Column::w][r] = 1.0F;
      chunk[Column::wavelength][r] = 0.0F; // assigned by the step
    }
    return chunk;
  };

  obdeect::propagation::ExtinctionSettings settings;
  settings.table = &table;
  settings.depths = &depths;
  settings.envelope = &envelope;
  settings.wavelengths = {250.0, 700.0};
  settings.seed = kSeed;
  settings.event = 1;
  settings.telescope = 2;
  settings.weighted = true;
  settings.curved = path;

  const auto curved_out = obdeect::propagation::apply_extinction_chunk(make_photons(), settings, 0);
  settings.curved.reset();
  const auto plane_out = obdeect::propagation::apply_extinction_chunk(make_photons(), settings, 0);
  CHECK(curved_out.size() == 1 && plane_out.size() == 1);
  if (curved_out.size() != 1 || plane_out.size() != 1) {
    return;
  }
  const auto &out = curved_out.front();
  const auto &plane = plane_out.front();
  CHECK(out.has(Column::p));
  double largest_difference = 0.0;
  for (std::size_t row = 0; row < out.rows(); ++row) {
    const double cos_t = std::sqrt(1.0 - std::pow(static_cast<double>(out[Column::cx][row]), 2));
    const auto iwl = static_cast<std::size_t>(out[Column::wavelength][row]);
    // p = transmission times the envelope (1 here).
    const double expected =
        table.transmission_curved<double>(out[Column::wavelength][row], out[Column::zem][row],
                                          cos_t, kTelescope, depths.od_focus[iwl], kEarthRadius);
    CHECK(std::abs(static_cast<double>(out[Column::p][row]) - expected) <= 1e-6);
    largest_difference =
        std::max(largest_difference, std::abs(static_cast<double>(out[Column::p][row]) -
                                              static_cast<double>(plane[Column::p][row])));
  }
  // Curvature changes the high-angle photons (the plane-parallel p is lower).
  CHECK(largest_difference > 1e-3);
}

} // namespace

int main() {
  return run_checks([] {
    test_flat_limit_fixture();
    test_flat_limit_and_brute_force_synthetic();
    test_curvature_at_large_angles();
    test_monotone_zero_and_validation();
    test_transmission_curved();
    test_acceptance_kernel();
    test_apply_extinction_curved();
  });
}
