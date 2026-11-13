#include "check.hpp"

#include "obdeect/propagation/scattering.hpp"
#include "obdeect/propagation/streams.hpp"
#include "obdeect/propagation/transmission.hpp"

#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using obdeect::stream::Column;

constexpr std::uint64_t seed = 20261107;
constexpr double lambda_nm = 355.0;

// Vertical beam, one telescope on the ground at horizontal distance D: the test geometry of
// the analytic integral and the step sum. Observation level obs_m, emission at h0_m.
constexpr double obs_m = 2000.0;
constexpr double h0_m = 50000.0;
constexpr double distance_m = 3000.0;
constexpr double radius_m = 12.5;

std::uint32_t mask_of(std::initializer_list<Column> columns) {
  std::uint32_t mask = 0;
  for (const Column c : columns) {
    mask |= obdeect::stream::column_bit(c);
  }
  return mask;
}

// Exponential molecular-only table as in the step-sum comparison: tau(h) = k exp(-h/scale),
// h in m a.s.l., nodes every 500 m from 1 to 60 km, wavelengths 354..356 nm.
std::string exponential_table_text(double k, double scale_m) {
  std::ostringstream text;
  text.precision(17);
  text << "# test table\n";
  for (int wl = 354; wl <= 356; ++wl) {
    for (int i = 0; i <= 118; ++i) {
      const double h_km = 1.0 + 0.5 * i;
      const double tau = k * std::exp(-h_km * 1000.0 / scale_m);
      text << wl << ' ' << h_km << ' ' << tau << " 0 1 0\n";
    }
  }
  return text.str();
}

// Site transmission table (obs 2000 m) shipped with the tests.
obdeect::propagation::TransmissionTable site_transmission() {
  return obdeect::propagation::TransmissionTable::from_file(
      std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/atm_trans_2000_1_10_0_0_2000.dat", obs_m);
}

// Analytic single-scattering integral (no sampling) for one downward photon of weight 1 from
// altitude h0 down to the observation level, into one telescope sphere on the ground:
//   Y = int dh exp(-(tau(h) - tau(h0)) [if attenuate]) beta_s(h) p_mix(theta) Omega(h) T(h)
// beta from the table (finite difference of the vertical optical depths), T and Omega as in
// scattering.hpp. Simpson's rule, 20000 intervals.
double analytic_yield(const obdeect::propagation::ScatteringTable &table,
                      const obdeect::propagation::TransmissionTable &tr, double lambda, double h0,
                      bool attenuate) {
  using obdeect::propagation::henyey_greenstein_phase;
  using obdeect::propagation::rayleigh_phase;
  const auto integrand = [&](double h) {
    const double dist = std::hypot(distance_m, h - obs_m);
    const double mu = (h - obs_m) / dist;
    const auto s = table.state<double>(lambda, h);
    const auto up = table.state<double>(lambda, h + 0.5);
    const auto down = table.state<double>(lambda, h - 0.5);
    const double beta_mol = down.tau_mol - up.tau_mol;
    const double beta_aer = down.tau_aer - up.tau_aer;
    const double phase =
        beta_mol * rayleigh_phase(mu) + s.ssa * beta_aer * henyey_greenstein_phase(mu, s.asymmetry);
    const double ratio = radius_m / dist;
    const double solid = 2.0 * std::numbers::pi * (1.0 - std::sqrt(1.0 - ratio * ratio));
    const double od = tr.optical_depth<double>(lambda, h) - tr.optical_depth<double>(lambda, obs_m);
    const double transmitted = std::exp(-od * dist / (h - obs_m));
    double attenuation = 1.0;
    if (attenuate) {
      const auto top = table.state<double>(lambda, h0);
      attenuation = std::exp(-((s.tau_mol + s.tau_aer) - (top.tau_mol + top.tau_aer)));
    }
    return phase * solid * transmitted * attenuation;
  };
  const int n = 20000;
  const double a = obs_m + 1e-3; // the integrand is 0/0 exactly at the observation level
  const double step = (h0_m - a) / n;
  double sum = integrand(a) + integrand(h0_m);
  for (int i = 1; i < n; ++i) {
    sum += (i % 2 == 1 ? 4.0 : 2.0) * integrand(a + i * step);
  }
  return sum * step / 3.0;
}

// Expected yield as a sum over beam steps (density hscale = exp(-z/7500 m), emission per step
// nphot_per_cm * hscale * steplength; Rayleigh angular distribution (3/16 pi)(1 + cos^2) of
// the same density; the sphere is hit with the solid angle Omega of the sphere). The sum
// neglects the attenuation of the beam before the scattering point. Per-step sum, 1 m steps,
// with q0 = w (K/H) exp(-obs/H) photons per metre at the observation level.
double step_sum_yield(double k, double scale_m, const obdeect::propagation::TransmissionTable &tr,
                      double lambda) {
  using obdeect::propagation::rayleigh_phase;
  const double q0 = (k / scale_m) * std::exp(-obs_m / scale_m);
  const double ds = 1.0;
  double sum = 0.0;
  for (double z = 0.0; z < h0_m - obs_m; z += ds) {
    const double z_mid = z + 0.5 * ds;
    const double h = obs_m + z_mid;
    const double hscale = std::exp(-z_mid / scale_m); // exp(-z/7500 m)
    const double dist = std::hypot(distance_m, z_mid);
    const double mu = z_mid / dist;
    const double ratio = radius_m / dist;
    const double solid = 2.0 * std::numbers::pi * (1.0 - std::sqrt(1.0 - ratio * ratio));
    const double od = tr.optical_depth<double>(lambda, h) - tr.optical_depth<double>(lambda, obs_m);
    const double transmitted = std::exp(-od * dist / z_mid);
    sum += q0 * hscale * ds * rayleigh_phase(mu) * solid * transmitted;
  }
  return sum;
}

// Per-primary scattered weight towards telescope 1: one vertical primary of weight 1 at
// altitude h0 for each first_index in [0, n).
std::vector<double> per_primary_totals(const obdeect::propagation::ScatteringTable &table,
                                       const obdeect::propagation::TransmissionTable &tr,
                                       double lambda, std::size_t n, std::uint32_t n_points) {
  obdeect::propagation::ScatteringConfig config;
  config.enabled = true;
  config.points_per_photon = n_points;
  const std::vector<obdeect::propagation::ScatterTelescope> sphere = {
      {1, distance_m, 0.0, 0.0, radius_m}};
  std::vector<double> totals(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    obdeect::stream::PhotonChunk chunk;
    chunk.frame = obdeect::stream::Frame::sky;
    chunk.first_index = i;
    chunk.resize(1, mask_of({Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::t,
                             Column::w, Column::wavelength}));
    chunk[Column::z][0] = static_cast<float>(h0_m - obs_m);
    chunk[Column::w][0] = 1.0F;
    chunk[Column::wavelength][0] = static_cast<float>(lambda);
    std::vector<std::vector<obdeect::stream::PhotonChunk>> lists;
    obdeect::propagation::scatter_chunk(config, table, tr, sphere, chunk, seed, 1, i, lists);
    for (const auto &piece : lists.at(0)) {
      for (const float w : piece[Column::w]) {
        totals[i] += w;
      }
    }
  }
  return totals;
}

void test_phase_functions() {
  // Normalisation over 4 pi: 2 pi int_{-1}^{1} p(mu) dmu = 1 (midpoint rule, 200000 steps).
  const int n = 200000;
  double norm_r = 0.0, norm_hg = 0.0, mean_hg = 0.0, mean_r = 0.0;
  const double g = 0.6;
  for (int i = 0; i < n; ++i) {
    const double mu = -1.0 + (i + 0.5) * 2.0 / n;
    const double dmu = 2.0 / n;
    norm_r += 2.0 * std::numbers::pi * obdeect::propagation::rayleigh_phase(mu) * dmu;
    mean_r += 2.0 * std::numbers::pi * mu * obdeect::propagation::rayleigh_phase(mu) * dmu;
    const double hg = obdeect::propagation::henyey_greenstein_phase(mu, g);
    norm_hg += 2.0 * std::numbers::pi * hg * dmu;
    mean_hg += 2.0 * std::numbers::pi * mu * hg * dmu;
  }
  CHECK_NEAR(norm_r, 1.0, 1e-8);
  CHECK_NEAR(norm_hg, 1.0, 1e-6); // peaked forward: midpoint rule converges more slowly
  CHECK_NEAR(mean_r, 0.0, 1e-8);
  CHECK_NEAR(mean_hg, g, 1e-5);
  // Mean cosine of the mixed phase with only aerosol scattering equals g.
  CHECK_NEAR(obdeect::propagation::mixed_phase(0.3, 0.0, 1.0, 1.0, g),
             obdeect::propagation::henyey_greenstein_phase(0.3, g), 1e-15);
}

void test_zero_optical_depth() {
  const obdeect::propagation::ScatteringTable table(exponential_table_text(0.0, 7500.0));
  const obdeect::propagation::TransmissionTable tr = site_transmission();
  obdeect::propagation::ScatteringConfig config;
  config.enabled = true;
  config.points_per_photon = 8;
  const std::vector<obdeect::propagation::ScatterTelescope> sphere = {{1, 3000.0, 0.0, 0.0, 12.5}};
  obdeect::stream::PhotonChunk chunk;
  chunk.frame = obdeect::stream::Frame::sky;
  chunk.resize(4, mask_of({Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::t,
                           Column::w, Column::wavelength}));
  for (std::size_t r = 0; r < 4; ++r) {
    chunk[Column::z][r] = 40000.0F;
    chunk[Column::w][r] = 1.0F;
    chunk[Column::wavelength][r] = 355.0F;
  }
  std::vector<std::vector<obdeect::stream::PhotonChunk>> lists;
  obdeect::propagation::scatter_chunk(config, table, tr, sphere, chunk, seed, 1, 0, lists);
  std::size_t rows = 0;
  for (const auto &piece : lists.at(0)) {
    rows += piece.rows();
  }
  CHECK(rows == 0);
  // Off by default: nothing is produced.
  obdeect::propagation::ScatteringConfig off;
  std::vector<std::vector<obdeect::stream::PhotonChunk>> none;
  obdeect::propagation::scatter_chunk(off, table, tr, sphere, chunk, seed, 1, 0, none);
  CHECK(none.empty());
}

// Analytic single-scattering integral against the Monte Carlo estimate (per-primary totals,
// standard error from their spread), for the exponential Rayleigh table, and against the
// step sum with the beam attenuation neglected.
void test_vertical_beam_against_integral_and_step_sum() {
  const double k = 0.02;
  const double scale = 7500.0;
  const obdeect::propagation::ScatteringTable table(exponential_table_text(k, scale));
  const obdeect::propagation::TransmissionTable tr = site_transmission();

  // Monte Carlo: 2000 primaries, 32 points each.
  const std::size_t n = 2000;
  const auto totals = per_primary_totals(table, tr, lambda_nm, n, 32);
  double mean = 0.0, var = 0.0;
  for (const double t : totals) {
    mean += t;
  }
  mean /= static_cast<double>(n);
  for (const double t : totals) {
    var += (t - mean) * (t - mean);
  }
  const double std_error = std::sqrt(var / static_cast<double>(n - 1) / static_cast<double>(n));

  const double with_attenuation = analytic_yield(table, tr, lambda_nm, h0_m, true);
  const double no_attenuation = analytic_yield(table, tr, lambda_nm, h0_m, false);
  const double step_sum = step_sum_yield(k, scale, tr, lambda_nm);
  std::fprintf(stderr,
               "vertical beam: MC %.6e +- %.2e, analytic %.6e, analytic w/o beam attenuation %.6e, "
               "step sum %.6e\n",
               mean, std_error, with_attenuation, no_attenuation, step_sum);
  CHECK(std::abs(mean - with_attenuation) <= 4.0 * std_error);
  // Step sum against the same integral without beam attenuation (1 m steps vs Simpson).
  CHECK(std::abs(step_sum - no_attenuation) <= 1e-3 * no_attenuation);
  // The step sum neglects the beam attenuation (here a few percent).
  CHECK(with_attenuation < no_attenuation);
}

// Example table (exponential atmosphere with molecular and aerosol scattering, g = 0.65). The
// segment discretisation error falls about as N_s^-2 (32 points: 3e-4 relative; 128: 1e-5).
void test_example_table_against_integral() {
  const obdeect::propagation::ScatteringTable table([] {
    std::ifstream file(std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) +
                       "/scattering_exponential.dat");
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
  }());
  const obdeect::propagation::TransmissionTable tr = site_transmission();
  const std::size_t n = 2000;
  const auto totals = per_primary_totals(table, tr, lambda_nm, n, 128);
  double mean = 0.0, var = 0.0;
  for (const double t : totals) {
    mean += t;
  }
  mean /= static_cast<double>(n);
  for (const double t : totals) {
    var += (t - mean) * (t - mean);
  }
  const double std_error = std::sqrt(var / static_cast<double>(n - 1) / static_cast<double>(n));
  const double analytic = analytic_yield(table, tr, lambda_nm, h0_m, true);
  std::fprintf(stderr, "example table: MC %.6e +- %.2e, analytic %.6e\n", mean, std_error,
               analytic);
  CHECK(std::abs(mean - analytic) <= 4.0 * std_error);
}

// Chunk splits and telescope order do not change the output: bit-identical.
void test_determinism_across_chunk_splits() {
  const obdeect::propagation::ScatteringTable table(exponential_table_text(0.05, 7500.0));
  const obdeect::propagation::TransmissionTable tr = site_transmission();
  obdeect::propagation::ScatteringConfig config;
  config.enabled = true;
  config.points_per_photon = 8;
  const std::vector<obdeect::propagation::ScatterTelescope> sphere = {
      {3, 3000.0, 0.0, 0.0, 12.5}, {9, -800.0, 1200.0, 0.0, 6.0}};
  const std::size_t n = 300;
  obdeect::stream::PhotonChunk all;
  all.frame = obdeect::stream::Frame::sky;
  all.resize(n, mask_of({Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::t,
                         Column::w, Column::wavelength}));
  for (std::size_t r = 0; r < n; ++r) {
    const double phi = 0.37 * static_cast<double>(r);
    all[Column::x][r] = static_cast<float>(200.0 * std::cos(phi));
    all[Column::y][r] = static_cast<float>(200.0 * std::sin(phi));
    all[Column::z][r] = static_cast<float>(5000.0 + 100.0 * static_cast<double>(r));
    all[Column::cx][r] = static_cast<float>(0.2 * std::cos(1.3 * phi));
    all[Column::cy][r] = static_cast<float>(0.2 * std::sin(1.3 * phi));
    all[Column::t][r] = static_cast<float>(r);
    all[Column::w][r] = 2.0F;
    all[Column::wavelength][r] = static_cast<float>(354 + static_cast<int>(r % 3));
  }

  const auto slice = [&](std::size_t begin, std::size_t end) {
    obdeect::stream::PhotonChunk part;
    part.frame = obdeect::stream::Frame::sky;
    part.first_index = begin;
    part.resize(end - begin, all.mask());
    for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c) {
      if (all.has(static_cast<Column>(c))) {
        for (std::size_t r = begin; r < end; ++r) {
          part.columns[c][r - begin] = all.columns[c][r];
        }
      }
    }
    return part;
  };
  const auto flatten = [](const std::vector<std::vector<obdeect::stream::PhotonChunk>> &lists) {
    std::vector<std::vector<float>> out(lists.size());
    for (std::size_t t = 0; t < lists.size(); ++t) {
      for (const auto &piece : lists[t]) {
        // Row-major, so the layout does not depend on how the rows are split into chunks.
        for (std::size_t row = 0; row < piece.rows(); ++row) {
          for (const auto column : obdeect::propagation::scatter_output_columns) {
            out[t].push_back(piece[column][row]);
          }
        }
      }
    }
    return out;
  };

  std::vector<std::vector<obdeect::stream::PhotonChunk>> whole;
  obdeect::propagation::scatter_chunk(config, table, tr, sphere, all, seed, 4, 0, whole);
  std::vector<std::vector<obdeect::stream::PhotonChunk>> split;
  obdeect::propagation::scatter_chunk(config, table, tr, sphere, slice(0, 70), seed, 4, 0, split);
  obdeect::propagation::scatter_chunk(config, table, tr, sphere, slice(70, 300), seed, 4, 70,
                                      split);

  const auto a = flatten(whole);
  const auto b = flatten(split);
  CHECK(a.size() == b.size());
  std::size_t total = 0;
  for (std::size_t t = 0; t < a.size() && t < b.size(); ++t) {
    CHECK(a[t] == b[t]);
    total += a[t].size();
  }
  CHECK(total > 0);
  // Two telescopes, same scatter points: one list each, telescope_id set.
  CHECK(whole.size() == 2);
  CHECK(whole[0].front().telescope_id == 3);
  CHECK(whole[1].front().telescope_id == 9);
}

// Invalid input is rejected: decreasing tau_mol with altitude is an error; so is a wavelength
// outside the table.
void test_rejects_invalid_input() {
  bool rejected = false;
  try {
    obdeect::propagation::ScatteringTable bad("340 1.0 0.1 0.0 1 0\n340 2.0 0.2 0.0 1 0\n");
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);

  const obdeect::propagation::ScatteringTable table(exponential_table_text(0.05, 7500.0));
  const obdeect::propagation::TransmissionTable tr = site_transmission();
  obdeect::propagation::ScatteringConfig config;
  config.enabled = true;
  const std::vector<obdeect::propagation::ScatterTelescope> sphere = {{1, 0.0, 0.0, 0.0, 12.5}};
  obdeect::stream::PhotonChunk chunk;
  chunk.frame = obdeect::stream::Frame::sky;
  chunk.resize(1, mask_of({Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::t,
                           Column::w, Column::wavelength}));
  chunk[Column::z][0] = 5000.0F;
  chunk[Column::w][0] = 1.0F;
  chunk[Column::wavelength][0] = 900.0F;
  std::vector<std::vector<obdeect::stream::PhotonChunk>> lists;
  rejected = false;
  try {
    obdeect::propagation::scatter_chunk(config, table, tr, sphere, chunk, seed, 1, 0, lists);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
}

// Scattered photons carry the time at their scattering point: t_emit + |emit - P| / c. The
// stage adds the rest of the path (checked in test_propagator.cpp).
void test_scattered_time_is_at_scattering_point() {
  const obdeect::propagation::ScatteringTable table(exponential_table_text(0.05, 7500.0));
  const obdeect::propagation::TransmissionTable tr = site_transmission();
  obdeect::propagation::ScatteringConfig config;
  config.enabled = true;
  config.points_per_photon = 16;
  const std::vector<obdeect::propagation::ScatterTelescope> sphere = {
      {1, distance_m, 0.0, 0.0, radius_m}};
  const double t_emit = 7.0;
  const double cx = 0.2;
  const double cy = 0.1;
  obdeect::stream::PhotonChunk chunk;
  chunk.frame = obdeect::stream::Frame::sky;
  chunk.resize(1, mask_of({Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::t,
                           Column::w, Column::wavelength}));
  chunk[Column::z][0] = static_cast<float>(h0_m - obs_m);
  chunk[Column::cx][0] = static_cast<float>(cx);
  chunk[Column::cy][0] = static_cast<float>(cy);
  chunk[Column::t][0] = static_cast<float>(t_emit);
  chunk[Column::w][0] = 1.0F;
  chunk[Column::wavelength][0] = static_cast<float>(lambda_nm);
  std::vector<std::vector<obdeect::stream::PhotonChunk>> lists;
  obdeect::propagation::scatter_chunk(config, table, tr, sphere, chunk, seed, 1, 0, lists);

  // Emission point: origin in x, y, at h0 (z above the observation level).
  const double emit_z = h0_m - obs_m;
  std::size_t rows = 0;
  double worst = 0.0;
  for (const auto &piece : lists.at(0)) {
    for (std::size_t r = 0; r < piece.rows(); ++r) {
      const double px = piece[Column::x][r];
      const double py = piece[Column::y][r];
      const double pz = piece[Column::z][r];
      const double path = std::sqrt(px * px + py * py + (emit_z - pz) * (emit_z - pz));
      const double expected = t_emit + path / obdeect::propagation::speed_of_light_m_per_ns;
      worst = std::max(worst, std::abs(static_cast<double>(piece[Column::t][r]) - expected));
      ++rows;
    }
  }
  CHECK(rows > 0);
  // float32 time near 1e5 ns resolves about 0.01 ns.
  CHECK(worst < 0.05);
}

// The state lookup rejects altitudes that are not finite and positive, and wavelengths that are
// not finite; far out-of-range wavelengths are clamped to the table edges.
void test_state_rejects_bad_arguments() {
  const obdeect::propagation::ScatteringTable table(exponential_table_text(0.05, 7500.0));
  const auto throws = [](auto &&f) {
    try {
      f();
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK(throws([&] { (void)table.state<double>(lambda_nm, -1000.0); }));
  CHECK(throws([&] { (void)table.state<double>(lambda_nm, 0.0); }));
  CHECK(throws([&] { (void)table.state<double>(lambda_nm, nan); }));
  CHECK(throws([&] { (void)table.state<double>(lambda_nm, inf); }));
  CHECK(throws([&] { (void)table.state<double>(nan, 5000.0); }));
  CHECK(!throws([&] { (void)table.state<double>(lambda_nm, 5000.0); }));

  // Out-of-range wavelengths (also huge ones, which overflow an integer cast) are clamped.
  const auto low = table.state<double>(-1e300, 5000.0);
  const auto first = table.state<double>(table.first_wavelength_nm(), 5000.0);
  const auto high = table.state<double>(1e300, 5000.0);
  const auto last = table.state<double>(table.last_wavelength_nm(), 5000.0);
  CHECK_NEAR(low.tau_mol, first.tau_mol, 0.0);
  CHECK_NEAR(high.tau_mol, last.tau_mol, 0.0);
}

// Missing documented columns and non-finite x, y, t are rejected.
void test_rejects_missing_and_non_finite_columns() {
  const obdeect::propagation::ScatteringTable table(exponential_table_text(0.05, 7500.0));
  const obdeect::propagation::TransmissionTable tr = site_transmission();
  obdeect::propagation::ScatteringConfig config;
  config.enabled = true;
  const std::vector<obdeect::propagation::ScatterTelescope> sphere = {
      {1, distance_m, 0.0, 0.0, radius_m}};
  const auto run = [&](std::initializer_list<Column> columns, Column poison, double value) {
    obdeect::stream::PhotonChunk chunk;
    chunk.frame = obdeect::stream::Frame::sky;
    chunk.resize(1, mask_of(columns));
    chunk[Column::z][0] = 5000.0F;
    chunk[Column::w][0] = 1.0F;
    chunk[Column::wavelength][0] = static_cast<float>(lambda_nm);
    if (chunk.has(poison)) {
      chunk[poison][0] = static_cast<float>(value);
    }
    std::vector<std::vector<obdeect::stream::PhotonChunk>> lists;
    try {
      obdeect::propagation::scatter_chunk(config, table, tr, sphere, chunk, seed, 1, 0, lists);
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  const std::initializer_list<Column> all = {Column::x,  Column::y, Column::z, Column::cx,
                                             Column::cy, Column::t, Column::w, Column::wavelength};
  CHECK(!run(all, Column::x, 0.0));
  for (const Column missing :
       {Column::x, Column::y, Column::cx, Column::cy, Column::t, Column::w}) {
    std::vector<Column> columns;
    for (const Column c : all) {
      if (c != missing) {
        columns.push_back(c);
      }
    }
    obdeect::stream::PhotonChunk chunk;
    chunk.frame = obdeect::stream::Frame::sky;
    std::uint32_t mask = 0;
    for (const Column c : columns) {
      mask |= obdeect::stream::column_bit(c);
    }
    chunk.resize(1, mask);
    chunk[Column::z][0] = 5000.0F;
    chunk[Column::wavelength][0] = static_cast<float>(lambda_nm);
    std::vector<std::vector<obdeect::stream::PhotonChunk>> lists;
    bool rejected = false;
    try {
      obdeect::propagation::scatter_chunk(config, table, tr, sphere, chunk, seed, 1, 0, lists);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    CHECK(rejected);
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK(run(all, Column::x, nan));
  CHECK(run(all, Column::y, inf));
  CHECK(run(all, Column::t, nan));
}

// Table wavelengths beyond int range are rejected before the integer conversion (1e300 is UB).
void test_rejects_huge_table_wavelength() {
  const auto throws_table = [](const std::string &text) {
    try {
      obdeect::propagation::ScatteringTable table(text);
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  // Two altitudes, tau decreasing with altitude; only the wavelength differs.
  const auto table_at = [](const std::string &wavelength) {
    return wavelength + " 1.0 0.1 0.0 1 0\n" + wavelength + " 2.0 0.05 0.0 1 0\n";
  };
  CHECK(throws_table(table_at("1e300")));
  CHECK(throws_table(table_at("2147483648")));
  CHECK(!throws_table(table_at("2147483647")));
}

// A telescope with an infinite radius is rejected; so is a NaN cz or zem column.
void test_rejects_infinite_radius_and_nan_columns() {
  const obdeect::propagation::ScatteringTable table(exponential_table_text(0.05, 7500.0));
  const obdeect::propagation::TransmissionTable tr = site_transmission();
  obdeect::propagation::ScatteringConfig config;
  config.enabled = true;
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::initializer_list<Column> base = {Column::x,  Column::y, Column::z, Column::cx,
                                              Column::cy, Column::t, Column::w, Column::wavelength};
  // One photon at z = 5 km with the given extra column set to value (if it is not zero-size).
  const auto make_chunk = [&](Column extra, double value) {
    obdeect::stream::PhotonChunk chunk;
    chunk.frame = obdeect::stream::Frame::sky;
    chunk.resize(1, mask_of(base) | obdeect::stream::column_bit(extra));
    chunk[Column::z][0] = 5000.0F;
    chunk[Column::w][0] = 1.0F;
    chunk[Column::wavelength][0] = static_cast<float>(lambda_nm);
    chunk[extra][0] = static_cast<float>(value);
    return chunk;
  };
  const auto throws_invalid = [](auto &&f) {
    try {
      f();
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  std::vector<std::vector<obdeect::stream::PhotonChunk>> lists;

  const std::vector<obdeect::propagation::ScatterTelescope> infinite = {
      {1, distance_m, 0.0, 0.0, inf}};
  const obdeect::stream::PhotonChunk plain = make_chunk(Column::cz, -1.0);
  CHECK(throws_invalid([&] {
    obdeect::propagation::scatter_chunk(config, table, tr, infinite, plain, seed, 1, 0, lists);
  }));

  const std::vector<obdeect::propagation::ScatterTelescope> sphere = {
      {1, distance_m, 0.0, 0.0, radius_m}};
  CHECK(throws_invalid([&] {
    obdeect::propagation::scatter_chunk(config, table, tr, sphere, make_chunk(Column::cz, nan),
                                        seed, 1, 0, lists);
  }));
  CHECK(throws_invalid([&] {
    obdeect::propagation::scatter_chunk(config, table, tr, sphere, make_chunk(Column::zem, nan),
                                        seed, 1, 0, lists);
  }));
}

} // namespace

int main() {
  return run_checks([] {
    test_phase_functions();
    test_zero_optical_depth();
    test_vertical_beam_against_integral_and_step_sum();
    test_example_table_against_integral();
    test_determinism_across_chunk_splits();
    test_rejects_invalid_input();
    test_scattered_time_is_at_scattering_point();
    test_state_rejects_bad_arguments();
    test_rejects_missing_and_non_finite_columns();
    test_rejects_huge_table_wavelength();
    test_rejects_infinite_radius_and_nan_columns();
  });
}
