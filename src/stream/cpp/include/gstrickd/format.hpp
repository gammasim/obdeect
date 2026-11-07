#pragma once

#include <cstdint>

namespace gstrickd {

// Photon-stream format constants (docs/format.md).
inline constexpr std::uint32_t record_magic = 0x52545347U; // "GSTR", little-endian.
inline constexpr int format_major = 1;
inline constexpr int format_minor = 0;
inline constexpr std::uint32_t max_chunk_rows = 1U << 16;

enum class RecordType : std::uint16_t {
  stream = 1,
  run = 2,
  event = 3,
  telescope = 4,
  photons = 5,
  telescope_end = 6,
  event_end = 7,
  run_end = 8,
  blob = 9,
  ledger = 10,
  stream_end = 15
};

enum class Frame : std::uint32_t { sky = 0, ground = 1, telescope = 2, camera = 3 };

enum class Column : std::uint32_t {
  x = 0,
  y = 1,
  z = 2,
  cx = 3,
  cy = 4,
  cz = 5,
  t = 6,
  zem = 7,
  w = 8,
  wavelength = 9,
  p = 10
};
inline constexpr std::uint32_t n_columns = 11;

[[nodiscard]] constexpr std::uint32_t column_bit(Column c) {
  return 1U << static_cast<std::uint32_t>(c);
}

inline constexpr std::uint32_t required_columns = column_bit(Column::x) | column_bit(Column::y) |
                                                  column_bit(Column::cx) | column_bit(Column::cy) |
                                                  column_bit(Column::t) | column_bit(Column::w);

[[nodiscard]] constexpr const char *column_name(Column c) {
  constexpr const char *names[n_columns] = {"x", "y",   "z", "cx",         "cy", "cz",
                                            "t", "zem", "w", "wavelength", "p"};
  return names[static_cast<std::uint32_t>(c)];
}

} // namespace gstrickd
