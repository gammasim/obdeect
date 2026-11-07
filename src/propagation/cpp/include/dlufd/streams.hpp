#pragma once

#include <cstdint>

namespace dlufd {

// Philox stage id of the atmosphere module (gstrickd::DrawKey::stage). Every module
// uses its own stage id, so draws of different modules never coincide.
inline constexpr std::uint32_t stage_atmosphere = 0x444C5546U; // "DLUF"

// Streams of draws within stage_atmosphere (gstrickd::DrawKey::stream).
// Wavelength draws use block 0 only; thinning draws use blocks 0..512 (at most
// 2048 photons per bunch), so the two streams never share a counter.
inline constexpr std::uint32_t stream_wavelength = 1;
inline constexpr std::uint32_t stream_thinning = 2;

} // namespace dlufd
