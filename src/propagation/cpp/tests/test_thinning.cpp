#include "check.hpp"

#include "dlufd/pipeline.hpp"
#include "dlufd/streams.hpp"
#include "dlufd/thinning.hpp"
#include "dlufd/transmission.hpp"
#include "dlufd/wavelength.hpp"

#include "gstrickd/photons.hpp"
#include "gstrickd/random.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using gstrickd::Column;

constexpr std::uint64_t seed = 20261107;

std::uint32_t mask_of(std::initializer_list<Column> columns) {
  std::uint32_t mask = 0;
  for (const Column c : columns) {
    mask |= gstrickd::column_bit(c);
  }
  return mask;
}

template <class Function> bool throws_invalid(Function function) {
  try {
    function();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

void test_envelope() {
  const dlufd::Envelope envelope({300.0, 400.0, 500.0}, {0.0, 0.5, 1.0});
  CHECK_NEAR(envelope(300.0), 0.0, 0.0);
  CHECK_NEAR(envelope(350.0), 0.25, 1e-15);
  CHECK_NEAR(envelope(450.0), 0.75, 1e-15);
  CHECK_NEAR(envelope(500.0), 1.0, 0.0);
  CHECK(throws_invalid([] { dlufd::Envelope bad({300.0, 300.0}, {0.1, 0.2}); }));
  CHECK(throws_invalid([] { dlufd::Envelope bad({300.0, 400.0}, {0.1, std::nan("")}); }));
  CHECK(throws_invalid([] { dlufd::Envelope bad({300.0, std::nan("")}, {0.1, 0.2}); }));
  CHECK(throws_invalid([] { dlufd::Envelope bad({300.0, 400.0}, {0.1, 1.5}); }));
}

// The wavelength range must be finite, with 1 <= lo < hi <= 2^24 (exact integers in float32).
void test_wavelength_range_validation() {
  dlufd::WavelengthRange{250.0, 700.0}.validate();
  dlufd::WavelengthRange{1.0, dlufd::max_exact_wavelength_nm}.validate();
  constexpr double inf = std::numeric_limits<double>::infinity();
  CHECK(throws_invalid([] { (void)dlufd::WavelengthRange{250.0, inf}.validate(); }));
  CHECK(throws_invalid([] { (void)dlufd::WavelengthRange{std::nan(""), 700.0}.validate(); }));
  CHECK(throws_invalid([] { (void)dlufd::WavelengthRange{-inf, 700.0}.validate(); }));
  CHECK(throws_invalid([] { (void)dlufd::WavelengthRange{0.5, 700.0}.validate(); }));
  CHECK(throws_invalid([] { (void)dlufd::WavelengthRange{700.0, 700.0}.validate(); }));
  CHECK(throws_invalid([] { (void)dlufd::WavelengthRange{250.0, 1e10}.validate(); }));
  // The largest allowed range draws exact, non-negative integers below 2^24.
  gstrickd::PhotonChunk chunk;
  chunk.resize(4096, mask_of({Column::x, Column::y, Column::cx, Column::cy, Column::t, Column::w,
                              Column::wavelength}));
  dlufd::assign_wavelengths(chunk, {1.0, dlufd::max_exact_wavelength_nm}, seed, 0, 0);
  bool exact = true;
  for (const float lambda : chunk[Column::wavelength]) {
    exact = exact && lambda == std::floor(lambda) && lambda >= 1.0F &&
            lambda < static_cast<float>(dlufd::max_exact_wavelength_nm);
  }
  CHECK(exact);
}

// Non-finite or out-of-range probabilities and weights are rejected, never thinned to zero.
void test_non_finite_acceptance_rejected() {
  const double nan = std::nan("");
  gstrickd::PhotonChunk chunk;
  chunk.resize(3, mask_of({Column::x, Column::y, Column::cx, Column::cy, Column::t, Column::w,
                           Column::wavelength}));
  for (std::size_t r = 0; r < 3; ++r) {
    chunk[Column::w][r] = 2.0F;
  }
  CHECK(throws_invalid([&] {
    (void)dlufd::thin_chunk(chunk, std::vector<double>{0.5, nan, 0.5}, seed, 0, 0, 0, true);
  }));
  CHECK(throws_invalid(
      [&] { (void)dlufd::thin_chunk(chunk, std::vector<double>{0.5, 0.5}, seed, 0, 0, 0, true); }));
  CHECK(throws_invalid([&] {
    (void)dlufd::thin_chunk(chunk, std::vector<double>{0.5, 1.5, 0.5}, seed, 0, 0, 0, true);
  }));
  gstrickd::PhotonChunk weighted = chunk;
  CHECK(throws_invalid(
      [&] { (void)dlufd::weight_chunk(weighted, std::vector<double>{0.5, nan, 0.5}, true); }));
  CHECK(throws_invalid(
      [&] { (void)dlufd::weight_chunk(weighted, std::vector<double>{0.5, 0.5}, true); }));
  // Nothing changed by the rejected calls.
  CHECK(weighted[Column::w][0] == 2.0F && weighted[Column::w][1] == 2.0F);
  CHECK(!weighted.has(Column::p));

  gstrickd::PhotonChunk bad_weight = chunk;
  bad_weight[Column::w][1] = std::numeric_limits<float>::quiet_NaN();
  CHECK(throws_invalid(
      [&] { (void)dlufd::weight_chunk(bad_weight, std::vector<double>{0.5, 0.5, 0.5}, false); }));
  bad_weight[Column::w][1] = std::numeric_limits<float>::infinity();
  CHECK(throws_invalid(
      [&] { (void)dlufd::weight_chunk(bad_weight, std::vector<double>{0.5, 0.5, 0.5}, false); }));
  CHECK(throws_invalid([&] {
    (void)dlufd::thin_chunk(bad_weight, std::vector<double>{0.5, 0.5, 0.5}, seed, 0, 0, 0, true);
  }));
}

// Wavelengths follow density ~ 1/lambda^2 on [lo, hi]: P(lambda < x) analytic.
void test_wavelength_distribution() {
  const dlufd::WavelengthRange range{200.0, 700.0};
  const std::size_t n = 65536;
  gstrickd::PhotonChunk chunk;
  chunk.first_index = 0;
  chunk.resize(n, mask_of({Column::x, Column::y, Column::cx, Column::cy, Column::t, Column::w,
                           Column::wavelength}));
  dlufd::assign_wavelengths(chunk, range, seed, 1, 2);
  std::size_t below = 0;
  bool integral = true;
  bool inside = true;
  for (const float lambda : chunk[Column::wavelength]) {
    integral = integral && lambda == std::floor(lambda);
    inside = inside && lambda >= 200.0F && lambda <= 700.0F;
    below += lambda < 300.0F ? 1U : 0U;
  }
  CHECK(integral);
  CHECK(inside);
  const double p = (1.0 / 200.0 - 1.0 / 300.0) / (1.0 / 200.0 - 1.0 / 700.0);
  const double sigma = std::sqrt(p * (1.0 - p) / static_cast<double>(n));
  CHECK_NEAR(static_cast<double>(below) / static_cast<double>(n), p, 4.0 * sigma);

  // Same keys, same wavelengths; existing wavelengths are kept.
  gstrickd::PhotonChunk again = chunk;
  again[Column::wavelength].assign(n, 0.0F);
  dlufd::assign_wavelengths(again, range, seed, 1, 2);
  CHECK(again[Column::wavelength] == chunk[Column::wavelength]);
  again[Column::wavelength][0] = 350.0F;
  dlufd::assign_wavelengths(again, range, seed, 1, 2);
  CHECK(again[Column::wavelength][0] == 350.0F);
}

// Binomial(floor(w), p) + Bernoulli(frac(w) p): mean and variance from the draws.
void test_surviving_photons_moments() {
  const double w = 3.5;
  const double p = 0.3;
  const std::uint64_t n = 400000;
  double sum = 0.0;
  double sum_sq = 0.0;
  for (std::uint64_t i = 0; i < n; ++i) {
    const gstrickd::DrawKey key{seed, dlufd::stage_atmosphere, 7, 3, i, dlufd::stream_thinning};
    const auto k = static_cast<double>(dlufd::surviving_photons(w, p, key));
    sum += k;
    sum_sq += k * k;
  }
  const double mean = sum / static_cast<double>(n);
  const double variance = sum_sq / static_cast<double>(n) - mean * mean;
  const double expected_mean = w * p;
  const double expected_variance = 3.0 * p * (1.0 - p) + 0.5 * p * (1.0 - 0.5 * p);
  CHECK_NEAR(mean, expected_mean, 4.0 * std::sqrt(expected_variance / static_cast<double>(n)));
  CHECK_NEAR(variance, expected_variance, 0.01);
  // Edge cases: p = 0 and p = 1 are deterministic.
  const gstrickd::DrawKey key{seed, dlufd::stage_atmosphere, 0, 0, 0, dlufd::stream_thinning};
  CHECK(dlufd::surviving_photons(5.0, 0.0, key) == 0);
  CHECK(dlufd::surviving_photons(5.0, 1.0, key) == 5);
  CHECK(dlufd::surviving_photons(5.5, 1.0, key) == 6);
}

// Build a chunk of n photons with deterministic columns, starting at photon `offset`.
gstrickd::PhotonChunk make_input(std::size_t offset, std::size_t n) {
  gstrickd::PhotonChunk chunk;
  chunk.first_index = offset;
  chunk.resize(n, mask_of({Column::x, Column::y, Column::cx, Column::cy, Column::t, Column::w,
                           Column::zem, Column::wavelength}));
  for (std::size_t r = 0; r < n; ++r) {
    const std::size_t i = offset + r;
    chunk[Column::x][r] = static_cast<float>(i % 97);
    chunk[Column::y][r] = static_cast<float>(i % 89);
    chunk[Column::cx][r] = static_cast<float>(0.4 * std::sin(0.001 * static_cast<double>(i)));
    chunk[Column::cy][r] = static_cast<float>(0.4 * std::cos(0.0013 * static_cast<double>(i)));
    chunk[Column::t][r] = static_cast<float>(i % 13);
    chunk[Column::w][r] = static_cast<float>(1.0 + static_cast<double>(i % 4) * 2.5);
    chunk[Column::zem][r] = static_cast<float>(2500.0 + static_cast<double>(i % 7) * 1000.0);
    chunk[Column::wavelength][r] = 0.0F;
  }
  return chunk;
}

// Thinned mode on 10^6 rows: the number of output photons is Binomial-sum with analytic
// mean and variance. Weighted mode gives the same mean.
void test_thinning_statistics() {
  const std::size_t total = 1000000;
  const std::size_t chunk_rows = gstrickd::max_chunk_rows;
  double expected_mean = 0.0;
  double expected_variance = 0.0;
  std::uint64_t thinned_count = 0;
  double weighted_sum = 0.0;
  std::uint64_t out_index = 0;
  for (std::size_t offset = 0; offset < total; offset += chunk_rows) {
    const std::size_t rows = std::min(chunk_rows, total - offset);
    gstrickd::PhotonChunk in = make_input(offset, rows);
    std::vector<double> p(rows);
    for (std::size_t r = 0; r < rows; ++r) {
      const double x = std::fmod(static_cast<double>(offset + r) * 0.6180339887, 1.0);
      p[r] = 0.1 + 0.8 * x;
      const double w = in[Column::w][r];
      const double whole = std::floor(w);
      const double frac = w - whole;
      expected_mean += w * p[r];
      expected_variance += whole * p[r] * (1.0 - p[r]) + frac * p[r] * (1.0 - frac * p[r]);
    }
    gstrickd::PhotonChunk weighted = in;
    dlufd::weight_chunk(weighted, p, false);
    for (const float w : weighted[Column::w]) {
      weighted_sum += static_cast<double>(w);
    }
    for (const auto &chunk : dlufd::thin_chunk(in, p, seed, 3, 5, out_index, true)) {
      thinned_count += chunk.rows();
      out_index += chunk.rows();
      for (const float w : chunk[Column::w]) {
        CHECK(w == 1.0F);
      }
    }
  }
  const double sigma = std::sqrt(expected_variance);
  CHECK_NEAR(static_cast<double>(thinned_count), expected_mean, 4.0 * sigma);
  CHECK_NEAR(weighted_sum, expected_mean, 1e-3 * expected_mean);
  CHECK_NEAR(static_cast<double>(thinned_count), weighted_sum, 4.0 * sigma);
  std::fprintf(stderr, "thinned %llu, expected %.1f +- %.1f, weighted sum %.1f\n",
               static_cast<unsigned long long>(thinned_count), expected_mean, sigma, weighted_sum);
}

// The output does not depend on how the input is split into chunks.
void test_chunk_size_determinism() {
  const auto table = dlufd::TransmissionTable::from_file(
      std::string(DLUFD_TEST_DATA_DIR) + "/atm_trans_2000_1_10_0_0_2000.dat", 2000.0);
  const auto depths = table.telescope_depths(2030.0, 28.0);
  const dlufd::Envelope envelope({200.0, 800.0}, {0.9, 0.3});
  dlufd::ExtinctionSettings settings;
  settings.table = &table;
  settings.depths = &depths;
  settings.envelope = &envelope;
  settings.wavelengths = {250.0, 700.0};
  settings.seed = seed;
  settings.event = 11;
  settings.telescope = 4;
  settings.weighted = false;

  const std::size_t total = 200000;
  auto run = [&](std::size_t slice_rows) {
    std::vector<std::vector<float>> columns(gstrickd::n_columns);
    std::uint64_t out_index = 0;
    for (std::size_t offset = 0; offset < total; offset += slice_rows) {
      const std::size_t rows = std::min(slice_rows, total - offset);
      auto outputs = dlufd::apply_extinction_chunk(make_input(offset, rows), settings, out_index);
      for (const auto &chunk : outputs) {
        for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
          if (chunk.has(static_cast<Column>(c))) {
            columns[c].insert(columns[c].end(), chunk.columns[c].begin(), chunk.columns[c].end());
          }
        }
        out_index += chunk.rows();
      }
    }
    return columns;
  };
  const auto small = run(1000);
  const auto large = run(gstrickd::max_chunk_rows);
  CHECK(!small[static_cast<std::uint32_t>(Column::x)].empty());
  for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
    CHECK(small[c] == large[c]);
  }
}

} // namespace

int main() {
  return run_checks([] {
    test_envelope();
    test_wavelength_range_validation();
    test_non_finite_acceptance_rejected();
    test_wavelength_distribution();
    test_surviving_photons_moments();
    test_thinning_statistics();
    test_chunk_size_determinism();
  });
}
