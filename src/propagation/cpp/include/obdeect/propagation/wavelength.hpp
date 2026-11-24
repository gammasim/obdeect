#pragma once

// Wavelength assignment for photons without a wavelength: density proportional to
// 1/lambda^2 on [lo, hi] nm, drawn by inverse CDF. The value is continuous. Transmission table
// lookups interpolate linearly between the integer-nm rows (default,
// WavelengthInterpolation::linear) or floor it to integer nm (WavelengthInterpolation::floor).

#include "obdeect/propagation/streams.hpp"

#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/random.hpp"
#include "obdeect/stream/wavelength.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace obdeect::propagation {

// Upper bound of the range; the wavelength column is float32.
inline constexpr double max_exact_wavelength_nm = 16777216.0; // 2^24

struct WavelengthRange {
  double lo_nm{};
  double hi_nm{};

  // Requires finite 1 <= lo <= hi <= 2^24 (nm): the lower bound keeps 1 / lo finite and
  // lookups at index >= 1, the upper bound is a float32-safe limit. lo == hi is a fixed
  // wavelength (a monochromatic light source), assigned as lo.
  void validate() const {
    if (!(std::isfinite(lo_nm) && std::isfinite(hi_nm) && lo_nm >= 1.0 && hi_nm >= lo_nm &&
          hi_nm <= max_exact_wavelength_nm)) {
      throw std::invalid_argument("wavelength range needs finite 1 <= lo <= hi <= 2^24 (nm)");
    }
  }
};

// Set the wavelength of every row whose wavelength is 0. The draw of row r uses the key
// (seed, stage_atmosphere, event, telescope, first_index + r, stream_wavelength), so it
// does not depend on how the photons are split into chunks.
// Rows [begin, end) only (default: all rows); rows are independent, so blocks may run in parallel.
inline void assign_wavelengths(obdeect::stream::PhotonChunk &chunk, const WavelengthRange &range,
                               std::uint64_t seed, std::uint32_t event, std::uint32_t telescope,
                               std::size_t begin = 0,
                               std::size_t end = std::numeric_limits<std::size_t>::max()) {
  auto &wavelength = chunk[obdeect::stream::Column::wavelength];
  end = std::min(end, chunk.rows());
  for (std::size_t row = begin; row < end; ++row) {
    if (wavelength[row] != 0.0F) {
      continue;
    }
    const obdeect::stream::DrawKey key{seed,      stage_atmosphere,        event,
                                       telescope, chunk.first_index + row, stream_wavelength};
    const double u = obdeect::stream::uniform4(key)[0];
    wavelength[row] = static_cast<float>(
        obdeect::stream::inverse_square_wavelength_nm(u, range.lo_nm, range.hi_nm));
  }
}

} // namespace obdeect::propagation
