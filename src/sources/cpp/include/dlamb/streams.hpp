#pragma once

// Random-number streams and the photon-stream format version used by dlamb.

#include "gstrickd/format.hpp"
#include "gstrickd/stages.hpp"

#include <cstdint>

namespace dlamb {

// Philox stage of all photon, count and track draws (gstrickd::stage_sources).
inline constexpr std::uint32_t stage_sources = gstrickd::stage_sources;

// Major version of the photon-stream format that dlamb writes.
inline constexpr int photon_stream_format_major = 1;

// Streams of one photon: each is one uniform4 call and is never shared.
enum class Stream : std::uint32_t {
  direction = 1,  // emission angle and azimuth
  target = 2,     // point on the target or entrance disc
  wavelength = 3, // wavelength
  pulse = 4,      // emission time
  count = 5,      // photon count (photon index 0 only)
};

} // namespace dlamb
