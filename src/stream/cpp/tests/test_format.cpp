#include "check.hpp"

#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/random.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

using namespace obdeect::stream;

namespace {

// Pearson chi-square of values in [0, 1) over 64 equal bins: 63 degrees of freedom, whose
// upper 1e-4 quantile is 113.5 (scipy.stats.chi2.isf(1e-4, 63)).
constexpr double chi_square_limit_64_bins = 113.5;

double chi_square_64(const std::vector<double> &u) {
  std::array<double, 64> bins{};
  for (const double x : u)
    bins[std::min<std::size_t>(63, static_cast<std::size_t>(x * 64.0))] += 1.0;
  const double expected = static_cast<double>(u.size()) / 64.0;
  double chi = 0.0;
  for (const double b : bins)
    chi += (b - expected) * (b - expected) / expected;
  return chi;
}

// Pearson correlation of two equally long sequences.
double correlation(const std::vector<double> &x, const std::vector<double> &y) {
  const double n = static_cast<double>(x.size());
  double mx = 0.0, my = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    mx += x[i] / n;
    my += y[i] / n;
  }
  double sxy = 0.0, sxx = 0.0, syy = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) * (x[i] - mx);
    syy += (y[i] - my) * (y[i] - my);
  }
  return sxy / std::sqrt(sxx * syy);
}

// Lag-1 correlation is consistent with zero when |r| < 5 / sqrt(n) (5 sigma).
bool lag1_is_uncorrelated(const std::vector<double> &u) {
  const std::vector<double> head(u.begin(), u.end() - 1);
  const std::vector<double> tail(u.begin() + 1, u.end());
  return std::abs(correlation(head, tail)) < 5.0 / std::sqrt(static_cast<double>(head.size()));
}

void chunk_columns() {
  PhotonChunk chunk;
  chunk.resize(3, required_columns | column_bit(Column::wavelength));
  CHECK(chunk.rows() == 3);
  CHECK(chunk.has(Column::wavelength));
  CHECK(!chunk.has(Column::zem));
  CHECK(chunk.mask() == (required_columns | column_bit(Column::wavelength)));
  chunk.validate();
  chunk[Column::zem].assign(2, 1.0F);
  CHECK(throws_with([&] { chunk.validate(); }, "columns differ in length"));
}

void missing_required_column() {
  PhotonChunk chunk;
  chunk.resize(2, column_bit(Column::x));
  CHECK(throws_with([&] { chunk.validate(); }, "lacks a required column"));
}

void resize_rejects_bad_input() {
  PhotonChunk chunk;
  chunk.resize(2, required_columns);
  CHECK(throws_with([&] { chunk.resize(std::size_t{max_chunk_rows} + 1, required_columns); },
                    "exceeds 2^16 rows"));
  CHECK(chunk.rows() == 2); // unchanged
  CHECK(throws_with([&] { chunk.resize(2, 1U << 31); }, "column mask has undefined bits"));
  CHECK(chunk.mask() == required_columns);
  chunk.resize(0, required_columns);
  CHECK(chunk.rows() == 0 && chunk.mask() == 0);
  chunk.validate();
}

void random_is_keyed() {
  const DrawKey a{.seed = 1, .stage = 2, .event = 3, .telescope = 4, .photon = 5, .stream = 0};
  DrawKey b = a;
  CHECK(uniform4(a) == uniform4(b));
  b.stage = 3;
  CHECK(uniform4(a) != uniform4(b));
  b = a;
  b.seed = 2;
  CHECK(uniform4(a) != uniform4(b));
  b = a;
  b.stream = 1;
  CHECK(uniform4(a) != uniform4(b));
  b = a;
  b.photon = 6;
  CHECK(uniform4(a) != uniform4(b));
  CHECK(uniform4(a, 0) != uniform4(a, 1));
  for (const double u : uniform4(a)) {
    CHECK(u > 0.0 && u < 1.0);
  }
}

void random_fields_disjoint() {
  // Regression: event bits must not alias block bits (event 0x10000 vs block 1).
  const DrawKey a{
      .seed = 1, .stage = 1, .event = 0x10000U, .telescope = 0, .photon = 0, .stream = 0};
  DrawKey b = a;
  b.event = 0;
  CHECK(uniform4(a, 0) != uniform4(b, 1));
  DrawKey c = a;
  c.telescope = 1U << 16;
  CHECK(throws_with<std::out_of_range>([&] { (void)uniform4(c); },
                                       "random-number counter field out of range"));
}

// Direction helpers: cz_is_consistent accepts either sign of cz; the sign helpers are strict.
void direction_helpers() {
  const double nan = std::nan("");
  CHECK(cz_is_consistent(-0.8, 0.6, 0.0));      // downward
  CHECK(cz_is_consistent(0.8, 0.6, 0.0));       // upward
  CHECK(cz_is_consistent(0.0, 1.0, 0.0));       // horizontal, cx^2 + cy^2 = 1
  CHECK(cz_is_consistent(0.0, 0.0, -1.0));      // horizontal, cy = -1
  CHECK(!cz_is_consistent(-0.5, 0.6, 0.0));     // magnitude inconsistent
  CHECK(!cz_is_consistent(0.5, 0.6, 0.0));      // magnitude inconsistent, upward
  CHECK(!cz_is_consistent(nan, 0.6, 0.0));      // NaN cz
  CHECK(!cz_is_consistent(-0.8, nan, 0.0));     // NaN cx
  CHECK(!cz_is_consistent(INFINITY, 0.0, 0.0)); // infinite cz

  CHECK(is_downward(-0.8));
  CHECK(!is_downward(0.8));
  CHECK(!is_downward(0.0));
  CHECK(!is_downward(-0.0));
  CHECK(!is_downward(nan));
  CHECK(is_upward(0.8));
  CHECK(!is_upward(-0.8));
  CHECK(!is_upward(0.0));
  CHECK(!is_upward(nan));

  // The snap keeps |cz| for upward photons too; downward and in-disc directions as before.
  float cx = 1.0000001F;
  float cy = 0.0F;
  snap_direction_to_cz(cx, cy, 0.8);
  CHECK_NEAR(static_cast<double>(cx) * cx + static_cast<double>(cy) * cy, 1.0 - 0.64, 1e-6);
  CHECK(cx > 0.0F && cy == 0.0F);
  float cx_in = 0.6F;
  float cy_in = 0.0F;
  snap_direction_to_cz(cx_in, cy_in, 0.8);
  CHECK(cx_in == 0.6F && cy_in == 0.0F);
}

// Mean, chi-square uniformity and lag-1 correlation of the draws of many photons of one key.
void random_statistics() {
  constexpr std::uint64_t n = std::uint64_t{1} << 20;
  std::vector<double> u(n);
  for (std::uint64_t i = 0; i < n; ++i) {
    const DrawKey key{.seed = 7, .stage = 1, .event = 0, .telescope = 0, .photon = i, .stream = 0};
    u[i] = uniform4(key)[0];
  }
  CHECK_NEAR(std::accumulate(u.begin(), u.end(), 0.0) / static_cast<double>(n), 0.5, 0.001);
  CHECK(chi_square_64(u) < chi_square_limit_64_bins);
  CHECK(lag1_is_uncorrelated(u));
}

// The same statistics over all four components and across keys: seeds, stages, events,
// telescopes and streams, each pooled over 2^14 photons.
void random_across_keys() {
  constexpr std::uint64_t photons = std::uint64_t{1} << 14;
  std::vector<double> pooled;
  pooled.reserve(photons * 64 * 4);
  std::vector<double> first_stream, second_stream;
  for (std::uint64_t seed = 0; seed < 64; ++seed) {
    for (std::uint64_t i = 0; i < photons; ++i) {
      const DrawKey key{.seed = seed,
                        .stage = static_cast<std::uint32_t>(seed % 3),
                        .event = static_cast<std::uint32_t>(seed / 8),
                        .telescope = 0,
                        .photon = i,
                        .stream = 0};
      for (const double x : uniform4(key))
        pooled.push_back(x);
    }
  }
  CHECK(chi_square_64(pooled) < chi_square_limit_64_bins);

  // Streams of one photon are independent: lag-1 across photons within each stream, and the
  // correlation between two streams of the same photon.
  constexpr std::uint64_t n = std::uint64_t{1} << 20;
  first_stream.reserve(n);
  second_stream.reserve(n);
  for (std::uint64_t i = 0; i < n; ++i) {
    const DrawKey a{.seed = 3, .stage = 2, .event = 1, .telescope = 4, .photon = i, .stream = 0};
    DrawKey b = a;
    b.stream = 1;
    first_stream.push_back(uniform4(a)[1]);
    second_stream.push_back(uniform4(b)[1]);
  }
  CHECK(lag1_is_uncorrelated(first_stream));
  CHECK(lag1_is_uncorrelated(second_stream));
  CHECK(std::abs(correlation(first_stream, second_stream)) <
        5.0 / std::sqrt(static_cast<double>(n)));
}

// Known-answer tests of Philox-4x32-10 (Random123 kat_vectors): the generator itself is pinned.
void philox_known_answers() {
  CHECK((philox4x32({0, 0, 0, 0}, {0, 0}) ==
         std::array<std::uint32_t, 4>{0x6627e8d5U, 0xe169c58dU, 0xbc57ac4cU, 0x9b00dbd8U}));
  CHECK((philox4x32({0xffffffffU, 0xffffffffU, 0xffffffffU, 0xffffffffU},
                    {0xffffffffU, 0xffffffffU}) ==
         std::array<std::uint32_t, 4>{0x408f276dU, 0x41c83b0eU, 0xa20bc7c6U, 0x6d5451fdU}));
  CHECK((philox4x32({0x243f6a88U, 0x85a308d3U, 0x13198a2eU, 0x03707344U},
                    {0xa4093822U, 0x299f31d0U}) ==
         std::array<std::uint32_t, 4>{0xd16cfe09U, 0x94fdccebU, 0x5001e420U, 0x24126ea1U}));
}

} // namespace

int main() {
  return run_checks([] {
    chunk_columns();
    missing_required_column();
    resize_rejects_bad_input();
    random_is_keyed();
    random_fields_disjoint();
    philox_known_answers();
    random_statistics();
    random_across_keys();
    direction_helpers();
  });
}
