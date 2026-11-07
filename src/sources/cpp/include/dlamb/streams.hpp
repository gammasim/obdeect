#pragma once

#include "gstrickd/stages.hpp"

#include <cstdint>

namespace dlamb {

// Stage id of dlamb draws (gstrickd::DrawKey::stage). Alias of the gstrickd constant.
inline constexpr std::uint32_t stage_sources = gstrickd::stage_sources;

// Photon-stream format major version that dlamb writes and reads. The C++ test
// gstrickd_contract checks it against gstrickd::format_major.
inline constexpr int photon_stream_format_major = 1;

// Streams of one photon (gstrickd::DrawKey::stream). Each stream is one uniform4 call,
// so the draws of a photon never depend on the other streams or on chunking.
enum class Stream : std::uint32_t {
  direction = 1, // Emission angle and azimuth (emitters without a target).
  target = 2,    // Point on a target disc (emitters) or on the entrance disc (stars).
  wavelength = 3,
  pulse = 4,
  count = 5 // Photon count of a source (photon index 0 only).
};

} // namespace dlamb
