#pragma once

#include "obdeect/stream/random.hpp"

#include <array>
#include <cstdint>

namespace obdeect::ray_tracing {

// Philox-4x32-10 (Salmon et al., SC'11), shared with the other components: a
// counter-based generator. A draw depends only on (seed, photon id, stream,
// block), never on the order in which photons are processed, so results are
// identical for any block size and thread count.
using obdeect::stream::philox4x32;

// Independent random streams of one photon.
// New values are appended; existing ones are never renumbered.
enum class Stream : std::uint32_t {
  source = 0,
  primary_scatter = 1,
  secondary_scatter = 2,
  thinning = 3,
  wavelength = 4
};

// Four uniform deviates in the open interval (0, 1).
[[nodiscard]] constexpr std::array<double, 4> uniform4(std::uint64_t seed, std::uint64_t photon_id,
                                                       Stream stream, std::uint32_t block = 0) {
  const auto r = philox4x32(
      {static_cast<std::uint32_t>(photon_id), static_cast<std::uint32_t>(photon_id >> 32),
       static_cast<std::uint32_t>(stream), block},
      {static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)});
  constexpr double scale = 1.0 / 4294967296.0;
  return {(r[0] + 0.5) * scale, (r[1] + 0.5) * scale, (r[2] + 0.5) * scale, (r[3] + 0.5) * scale};
}

} // namespace obdeect::ray_tracing
