#include "check.hpp"

#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/random.hpp"

#include <cmath>
#include <stdexcept>

using namespace obdeect::stream;

namespace {

void chunk_columns() {
  PhotonChunk chunk;
  chunk.resize(3, required_columns | column_bit(Column::wavelength));
  CHECK(chunk.rows() == 3);
  CHECK(chunk.has(Column::wavelength));
  CHECK(!chunk.has(Column::zem));
  CHECK(chunk.mask() == (required_columns | column_bit(Column::wavelength)));
  chunk.validate();
  chunk[Column::zem].assign(2, 1.0F);
  bool threw = false;
  try {
    chunk.validate();
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
}

void missing_required_column() {
  PhotonChunk chunk;
  chunk.resize(2, column_bit(Column::x));
  bool threw = false;
  try {
    chunk.validate();
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
}

void resize_rejects_bad_input() {
  PhotonChunk chunk;
  chunk.resize(2, required_columns);
  bool threw = false;
  try {
    chunk.resize(std::size_t{max_chunk_rows} + 1, required_columns);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
  CHECK(chunk.rows() == 2); // unchanged
  threw = false;
  try {
    chunk.resize(2, 1U << 31);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
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
  bool threw = false;
  try {
    (void)uniform4(c);
  } catch (const std::out_of_range &) {
    threw = true;
  }
  CHECK(threw);
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

void random_mean() {
  double sum = 0.0;
  constexpr int n = 100000;
  for (int i = 0; i < n; ++i) {
    const DrawKey key{.seed = 7,
                      .stage = 1,
                      .event = 0,
                      .telescope = 0,
                      .photon = static_cast<std::uint64_t>(i),
                      .stream = 0};
    sum += uniform4(key)[0];
  }
  CHECK_NEAR(sum / n, 0.5, 0.005);
}

} // namespace

int main() {
  return run_checks([] {
    chunk_columns();
    missing_required_column();
    resize_rejects_bad_input();
    random_is_keyed();
    random_fields_disjoint();
    random_mean();
    direction_helpers();
  });
}
