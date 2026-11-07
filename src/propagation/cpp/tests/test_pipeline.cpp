#include "check.hpp"

#include "dlufd/extinction.hpp"
#include "dlufd/pipeline.hpp"
#include "dlufd/streams.hpp"
#include "dlufd/thinning.hpp"
#include "dlufd/transmission.hpp"
#include "dlufd/wavelength.hpp"

#include "gstrickd/photons.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using gstrickd::Column;

constexpr std::uint64_t seed = 20261107;

template <class Function> bool throws_invalid(Function function) {
  try {
    function();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

// Constant optical depth 0.5 above 3 km, zero at the observation level (one row, 400 nm;
// other wavelengths are clamped to it).
dlufd::TransmissionTable constant_table() {
  return dlufd::TransmissionTable("# H2= 2.000, H1=    3.000     4.000\n400 0.5 0.5\n", 2000.0);
}

// n photons, straight down (cx = cy = 0, cz = -1), zem 3500 m, weight 1, all at `wavelength`.
gstrickd::PhotonChunk photons(std::size_t n, float wavelength) {
  gstrickd::PhotonChunk chunk;
  chunk.first_index = 0;
  chunk.resize(n, gstrickd::column_bit(Column::x) | gstrickd::column_bit(Column::y) |
                      gstrickd::column_bit(Column::cx) | gstrickd::column_bit(Column::cy) |
                      gstrickd::column_bit(Column::t) | gstrickd::column_bit(Column::w) |
                      gstrickd::column_bit(Column::zem) | gstrickd::column_bit(Column::wavelength));
  for (std::size_t r = 0; r < n; ++r) {
    chunk[Column::zem][r] = 3500.0F;
    chunk[Column::w][r] = 1.0F;
    chunk[Column::wavelength][r] = wavelength;
  }
  return chunk;
}

dlufd::ExtinctionSettings settings_for(const dlufd::TransmissionTable &table,
                                       const dlufd::TelescopeDepths &depths,
                                       const dlufd::Envelope &envelope, bool weighted) {
  dlufd::ExtinctionSettings settings;
  settings.table = &table;
  settings.depths = &depths;
  settings.envelope = &envelope;
  settings.wavelengths = {250.0, 700.0};
  settings.seed = seed;
  settings.event = 1;
  settings.telescope = 2;
  settings.weighted = weighted;
  return settings;
}

std::size_t total_rows(const std::vector<gstrickd::PhotonChunk> &chunks) {
  std::size_t rows = 0;
  for (const auto &chunk : chunks) {
    rows += chunk.rows();
  }
  return rows;
}

// cz is checked for finiteness as well as consistency with cx, cy.
void test_cz_column_must_be_finite() {
  const auto table = constant_table();
  const auto depths = table.telescope_depths(2030.0, 28.0);
  const dlufd::Envelope envelope({200.0, 800.0}, {0.9, 0.3});
  const auto settings = settings_for(table, depths, envelope, true);

  auto consistent = photons(4, 400.0F);
  consistent[Column::cz].assign(4, -1.0F);
  CHECK(!throws_invalid([&] { (void)dlufd::apply_extinction_chunk(consistent, settings, 0); }));

  for (const float bad :
       {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity(), 0.5F}) {
    auto chunk = photons(4, 400.0F);
    chunk[Column::cz].assign(4, -1.0F);
    chunk[Column::cz][2] = bad;
    CHECK(throws_invalid([&] { (void)dlufd::apply_extinction_chunk(chunk, settings, 0); }));
  }
}

// Wavelengths are finite, non-negative integers (0 means not assigned). 400.5, NaN, inf and
// negative values are rejected before any conversion to int.
void test_wavelength_must_be_integer() {
  const auto table = constant_table();
  const auto depths = table.telescope_depths(2030.0, 28.0);
  const dlufd::Envelope envelope({200.0, 800.0}, {0.9, 0.3});
  const auto settings = settings_for(table, depths, envelope, true);

  for (const float bad : {400.5F, std::numeric_limits<float>::quiet_NaN(),
                          std::numeric_limits<float>::infinity(), -1.0F}) {
    auto chunk = photons(3, 400.0F);
    chunk[Column::wavelength][1] = bad;
    CHECK(throws_invalid([&] { (void)dlufd::apply_extinction_chunk(chunk, settings, 0); }));
  }
  // Integer wavelengths and unassigned rows are accepted.
  auto valid = photons(3, 400.0F);
  valid[Column::wavelength][1] = 0.0F;
  CHECK(!throws_invalid([&] { (void)dlufd::apply_extinction_chunk(valid, settings, 0); }));
}

// Wavelengths at or above max_lambda_nm, including values beyond the int range, are dropped
// without any conversion to int: p = 0 and no thinned photon. A NaN or inf is rejected.
void test_wavelength_above_max_is_dropped() {
  const auto table = constant_table();
  const auto depths = table.telescope_depths(2030.0, 28.0);
  const dlufd::Envelope envelope({200.0, 2000.0}, {0.9, 0.9});
  for (const float lambda : {1000.0F, 1e10F, 3.0e38F}) {
    auto thinned = photons(50, lambda);
    CHECK(total_rows(dlufd::apply_extinction_chunk(
              thinned, settings_for(table, depths, envelope, false), 0)) == 0);

    auto weighted = photons(50, lambda);
    const auto out =
        dlufd::apply_extinction_chunk(weighted, settings_for(table, depths, envelope, true), 0);
    CHECK(out.size() == 1);
    for (const float w : out.front()[Column::w]) {
      CHECK(w == 0.0F);
    }
    for (const float p : out.front()[Column::p]) {
      CHECK(p == 0.0F);
    }
  }
  // 999 nm is still traced: the acceptance is not zero at the table edge.
  auto below = photons(50, 999.0F);
  CHECK(total_rows(dlufd::apply_extinction_chunk(
            below, settings_for(table, depths, envelope, false), 0)) > 0);
}

// The acceptance lookup applies the same rules directly.
void test_atmospheric_acceptance_wavelength() {
  const auto table = constant_table();
  const auto depths = table.telescope_depths(2030.0, 28.0);
  std::vector<double> p(2, -1.0);
  auto chunk = photons(2, 400.0F);
  chunk[Column::wavelength][1] = 1e10F;
  dlufd::atmospheric_acceptance<double>(chunk, table, depths, p);
  CHECK(p[0] > 0.0 && p[0] < 1.0);
  CHECK(p[1] == 0.0);

  chunk[Column::wavelength][1] = 400.5F;
  CHECK(throws_invalid(
      [&] { (void)dlufd::atmospheric_acceptance<double>(chunk, table, depths, p); }));
  chunk[Column::wavelength][1] = std::numeric_limits<float>::quiet_NaN();
  CHECK(throws_invalid(
      [&] { (void)dlufd::atmospheric_acceptance<double>(chunk, table, depths, p); }));
}

} // namespace

int main() {
  return run_checks([] {
    test_cz_column_must_be_finite();
    test_wavelength_must_be_integer();
    test_wavelength_above_max_is_dropped();
    test_atmospheric_acceptance_wavelength();
  });
}
