#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

namespace gstrickd {

// Philox stage ids (DrawKey::stage): one per module, so draws of different modules never
// coincide. The values are part of the results: never change a value in use. A new module
// takes the next free value and is added to stage_ids below.
inline constexpr std::uint32_t stage_sources = 0x444C0001U;    // dlamb: light sources
inline constexpr std::uint32_t stage_programs = 0x444C0002U;   // dlamb: calibration programs
inline constexpr std::uint32_t stage_obdeect = 0x444C0003U;    // obdeect: reserved (see docs)
inline constexpr std::uint32_t stage_atmosphere = 0x444C5546U; // dlufd: atmosphere ("DLUF")

namespace detail {

// True when no two entries are equal.
template <std::size_t N> constexpr bool all_distinct(const std::array<std::uint32_t, N> &ids) {
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = i + 1; j < N; ++j) {
      if (ids[i] == ids[j]) {
        return false;
      }
    }
  }
  return true;
}

} // namespace detail

inline constexpr std::array<std::uint32_t, 4> stage_ids{stage_sources, stage_programs,
                                                        stage_obdeect, stage_atmosphere};
static_assert(detail::all_distinct(stage_ids), "gstrickd: stage ids must be distinct");

} // namespace gstrickd
