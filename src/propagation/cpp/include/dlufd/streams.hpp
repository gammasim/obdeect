#pragma once

#include <cstdint>

namespace dlufd {

// Philox stage id of the atmosphere module (gstrickd::DrawKey::stage). Every module
// uses its own stage id, so draws of different modules never coincide.
inline constexpr std::uint32_t stage_atmosphere = 0x444C5546U; // "DLUF"

// Streams of draws within stage_atmosphere (gstrickd::DrawKey::stream).
// Wavelength draws use block 0 only; thinning draws use blocks 0..512 (at most
// 2048 photons per bunch), so the streams never share a counter.
inline constexpr std::uint32_t stream_wavelength = 1;
inline constexpr std::uint32_t stream_thinning = 2;
// Single scattering (scattering.hpp): deviates 0..N-1 place the scattering points, deviates
// N.. are the aim-point offsets; one key per primary photon.
inline constexpr std::uint32_t stream_scattering = 3;
// Thinning of scattered photons (envelope only; scattered photons carry their transmission in w).
// Their keys use the index of the scattered photon in its telescope's list, not the primary index.
inline constexpr std::uint32_t stream_scattered_thinning = 4;

} // namespace dlufd
