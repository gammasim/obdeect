#pragma once

// Wavelength assignment for photons without a wavelength (sim_telarray.c:2821):
// density proportional to 1/lambda^2 on [lo, hi] nm, drawn by inverse CDF. As in
// sim_telarray, the value is truncated to an integer number of nm (floor, as the value is
// positive).

#include "dlufd/streams.hpp"

#include "gstrickd/photons.hpp"
#include "gstrickd/random.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace dlufd {

// Integers are exact in float32 up to 2^24; the wavelength column is float32.
inline constexpr double max_exact_wavelength_nm = 16777216.0; // 2^24

struct WavelengthRange {
  double lo_nm{};
  double hi_nm{};

  // Requires finite 1 <= lo < hi <= 2^24 (nm): the lower bound keeps the integer wavelength
  // >= 1, the upper bound keeps the drawn integers exact in float32.
  void validate() const {
    if (!(std::isfinite(lo_nm) && std::isfinite(hi_nm) && lo_nm >= 1.0 && hi_nm > lo_nm &&
          hi_nm <= max_exact_wavelength_nm)) {
      throw std::invalid_argument("wavelength range needs finite 1 <= lo < hi <= 2^24 (nm)");
    }
  }
};

// Inverse CDF of density 1/lambda^2 on [lo, hi]; u in [0, 1).
[[nodiscard]] inline double wavelength_from_uniform(double u, const WavelengthRange &range) {
  return 1.0 / (1.0 / range.lo_nm - u * (1.0 / range.lo_nm - 1.0 / range.hi_nm));
}

// Set the wavelength of every row whose wavelength is 0. The draw of row r uses the key
// (seed, stage_atmosphere, event, telescope, first_index + r, stream_wavelength), so it
// does not depend on how the photons are split into chunks.
inline void assign_wavelengths(gstrickd::PhotonChunk &chunk, const WavelengthRange &range,
                               std::uint64_t seed, std::uint32_t event, std::uint32_t telescope) {
  auto &wavelength = chunk[gstrickd::Column::wavelength];
  for (std::size_t row = 0; row < chunk.rows(); ++row) {
    if (wavelength[row] != 0.0F) {
      continue;
    }
    const gstrickd::DrawKey key{seed,      stage_atmosphere,        event,
                                telescope, chunk.first_index + row, stream_wavelength};
    const double u = gstrickd::uniform4(key)[0];
    wavelength[row] = static_cast<float>(std::floor(wavelength_from_uniform(u, range)));
  }
}

} // namespace dlufd
