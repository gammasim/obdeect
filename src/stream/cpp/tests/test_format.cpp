#include "check.hpp"

#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/random.hpp"

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
  });
}
