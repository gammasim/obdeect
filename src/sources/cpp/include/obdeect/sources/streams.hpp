#pragma once

// Random-number streams and the photon-stream format version used by obdeect sources.

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/stages.hpp"

#include <cstdint>

namespace obdeect::sources {

// Philox stage of all photon, count and track draws (obdeect::stream::stage_sources).
inline constexpr std::uint32_t stage_sources = obdeect::stream::stage_sources;

// Major version of the photon-stream format that obdeect sources writes.
inline constexpr int photon_stream_format_major = 1;

// Streams of one photon: each is one uniform4 call and is never shared.
enum class Stream : std::uint32_t {
  direction = 1,  // emission angle and azimuth
  target = 2,     // point on the target or entrance disc
  wavelength = 3, // wavelength
  pulse = 4,      // emission time
  count = 5,      // photon count (photon index 0 only)
};

} // namespace obdeect::sources
