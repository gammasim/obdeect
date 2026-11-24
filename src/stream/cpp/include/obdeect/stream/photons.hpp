#pragma once

#include "obdeect/stream/format.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace obdeect::stream {

// Metadata value: int, float, string or float array (docs/stream/format.md, "Metadata").
using Value = std::variant<std::int64_t, double, std::string, std::vector<double>>;
using Metadata = std::map<std::string, Value>;

// Photons of one telescope in one frame, as float32 columns (SoA). Absent
// columns are empty; present columns all have rows() entries. A chunk without
// rows has no columns (column mask 0), as in the stream format.
struct PhotonChunk {
  std::uint32_t telescope_id{};
  Frame frame{Frame::ground};
  std::uint64_t first_index{};
  std::array<std::vector<float>, n_columns> columns{};

  [[nodiscard]] std::vector<float> &operator[](Column c) {
    return columns[static_cast<std::uint32_t>(c)];
  }
  [[nodiscard]] const std::vector<float> &operator[](Column c) const {
    return columns[static_cast<std::uint32_t>(c)];
  }
  [[nodiscard]] bool has(Column c) const { return !(*this)[c].empty(); }
  [[nodiscard]] std::size_t rows() const { return rows_; }

  // Column mask of the present (non-empty) columns.
  [[nodiscard]] std::uint32_t mask() const {
    std::uint32_t m = 0;
    for (std::uint32_t i = 0; i < n_columns; ++i) {
      if (!columns[i].empty()) {
        m |= 1U << i;
      }
    }
    return m;
  }

  // Size all columns in `column_mask` to n rows (zero-filled) and clear the
  // others; n = 0 clears all columns. Throws std::invalid_argument, before
  // changing the chunk, for n > 2^16 or bits of undefined columns.
  void resize(std::size_t n, std::uint32_t column_mask) {
    if (n > max_chunk_rows) {
      throw std::invalid_argument("photon chunk exceeds 2^16 rows");
    }
    if ((column_mask >> n_columns) != 0) {
      throw std::invalid_argument("photon chunk column mask has undefined bits");
    }
    rows_ = n;
    for (std::uint32_t i = 0; i < n_columns; ++i) {
      if (n > 0 && ((column_mask >> i) & 1U)) {
        columns[i].assign(n, 0.0F);
      } else {
        columns[i].clear();
      }
    }
  }

  // Clears all columns and sets rows() to n: the chunk of a photons record whose columns were not
  // decoded (Reader::keep_photons_of). validate() rejects such a chunk; it is not written.
  void set_rows_without_columns(std::size_t n) {
    if (n > max_chunk_rows) {
      throw std::invalid_argument("photon chunk exceeds 2^16 rows");
    }
    rows_ = n;
    for (auto &column : columns) {
      column.clear();
    }
  }

  // Keep the first n rows of every present column (n <= rows()).
  void truncate(std::size_t n) {
    if (n > rows_) {
      throw std::invalid_argument("photon chunk truncate beyond its rows");
    }
    rows_ = n;
    for (auto &column : columns) {
      if (!column.empty()) {
        column.resize(n);
      }
    }
  }

  // Throw unless all required columns are present with equal length.
  void validate() const {
    if ((mask() & required_columns) != required_columns && rows_ > 0) {
      throw std::invalid_argument("photon chunk lacks a required column (x y cx cy t w)");
    }
    for (const auto &column : columns) {
      if (!column.empty() && column.size() != rows_) {
        throw std::invalid_argument("photon chunk columns differ in length");
      }
    }
    if (rows_ > max_chunk_rows) {
      throw std::invalid_argument("photon chunk exceeds 2^16 rows");
    }
    if (has(Column::path)) {
      for (const float path : (*this)[Column::path]) {
        if (!std::isfinite(path) || path < 0.0F) {
          throw std::invalid_argument("photon chunk path must be finite and >= 0");
        }
      }
    }
  }

private:
  std::size_t rows_{};
};

// Largest allowed |cz^2 - (1 - cx^2 - cy^2)| for a stored direction. Float32 columns hold a unit
// vector to about 1e-7; the squared form keeps that rounding from failing nearly horizontal
// photons, for which an error in cz of the same size is amplified by 1 / (2 |cz|).
inline constexpr double direction_norm_tolerance = 1e-5;

// True if (cx, cy, cz) is a unit vector within the tolerance, for cz of either sign (cz = 0 is
// horizontal). False for any non-finite cz and for NaN.
[[nodiscard]] inline bool cz_is_consistent(double cz, double cx, double cy) {
  return std::isfinite(cz) &&
         std::abs(cz * cz - (1.0 - cx * cx - cy * cy)) <= direction_norm_tolerance;
}

// True if the photon moves downwards (cz < 0). False for upwards, horizontal and NaN.
[[nodiscard]] inline bool is_downward(double cz) { return cz < 0.0; }

// True if the photon moves upwards (cz > 0). False for downwards, horizontal and NaN.
[[nodiscard]] inline bool is_upward(double cz) { return cz > 0.0; }

// True if (cx, cy) lies in the unit disc within the tolerance. A stored direction can round to
// cx^2 + cy^2 slightly above 1 when cz is tiny; NaN is rejected. Use max(0, 1 - cx^2 - cy^2)
// for the magnitude of cz.
[[nodiscard]] inline bool direction_in_unit_disc(double cx, double cy) {
  return cx * cx + cy * cy <= 1.0 + direction_norm_tolerance;
}

// Reader checks of a decoded chunk beyond PhotonChunk::validate() (docs/stream/format.md): every
// column value is finite (x, y, t and w in particular); the direction is a unit vector, i.e. cz is
// consistent with cx, cy, or without a cz column (cx, cy) lies in the unit disc; and a photon with
// cz > 0 (upward) is only accepted if upward_allowed (format 1.1 and later). Throws
// std::invalid_argument. Writers do not call this: they may write what the reader rejects.
inline void validate_photon_values(const PhotonChunk &chunk, bool upward_allowed) {
  chunk.validate();
  for (const auto &column : chunk.columns) {
    for (const float value : column) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument("photon column holds a value that is not finite");
      }
    }
  }
  const bool has_cz = chunk.has(Column::cz);
  for (std::size_t i = 0; i < chunk.rows(); ++i) {
    const double cx = chunk[Column::cx][i];
    const double cy = chunk[Column::cy][i];
    if (!has_cz) {
      if (!direction_in_unit_disc(cx, cy)) {
        throw std::invalid_argument("photon direction (cx, cy) lies outside the unit disc");
      }
      continue;
    }
    const double cz = chunk[Column::cz][i];
    if (!cz_is_consistent(cz, cx, cy)) {
      throw std::invalid_argument("photon direction (cx, cy, cz) is not a unit vector");
    }
    if (!upward_allowed && is_upward(cz)) {
      throw std::invalid_argument("upward photon (cz > 0) in a stream of format 1.0");
    }
  }
}

// Rescale (cx, cy) of a photon with any finite cz, |cz| <= 1, whose cx^2 + cy^2 rounded above
// 1, so that cx^2 + cy^2 = 1 - cz^2 and the magnitude of the validated cz is kept (the sign of cz
// does not matter). Other directions are unchanged.
inline void snap_direction_to_cz(float &cx, float &cy, double cz) {
  const double r2 = static_cast<double>(cx) * cx + static_cast<double>(cy) * cy;
  if (std::isfinite(cz) && cz * cz <= 1.0 && r2 > 1.0) {
    const double scale = std::sqrt((1.0 - cz * cz) / r2);
    cx = static_cast<float>(cx * scale);
    cy = static_cast<float>(cy * scale);
  }
}

// Ground representation of a run (docs/stream/format.md, "Metadata"): `plane` (default) has rows at
// z = 0 with x, y relative to the telescope sphere centre; `sphere` has 3D rows on or inside the
// telescope sphere. Throws std::runtime_error for any other value or type.
enum class GroundRepresentation : std::uint8_t { plane, sphere };

[[nodiscard]] inline GroundRepresentation ground_representation(const Metadata &run) {
  const auto it = run.find("ground_representation");
  if (it == run.end()) {
    return GroundRepresentation::plane;
  }
  const auto *name = std::get_if<std::string>(&it->second);
  if (name == nullptr) {
    throw std::runtime_error("obdeect stream: ground_representation must be a string");
  }
  if (*name == "plane") {
    return GroundRepresentation::plane;
  }
  if (*name == "sphere") {
    return GroundRepresentation::sphere;
  }
  throw std::runtime_error("obdeect stream: ground_representation must be \"plane\" or "
                           "\"sphere\", not \"" +
                           *name + "\"");
}

// Tolerances for |r| against the telescope sphere radius R of sphere rows: a row is on the sphere
// if |r| >= R (1 - 1e-6) - 1e-4 m (rounding of float32 positions), and no row may lie outside
// R (1 + 1e-6) + 1e-4 m.
inline constexpr double sphere_relative_tolerance = 1e-6;
inline constexpr double sphere_absolute_tolerance_m = 1e-4;

// Throws std::invalid_argument unless the chunk is a valid set of sphere rows (G frame, see
// docs/stream/format.md) for a telescope sphere of radius radius_m: x y z cx cy cz t zem w present
// (z and cz included) and finite, cz consistent with (cx, cy), and |r| <= R (1 + 1e-6) + 1e-4 m for
// every row. Optional columns (path, wavelength) are checked by PhotonChunk::validate().
inline void validate_sphere_rows(const PhotonChunk &chunk, double radius_m) {
  if (chunk.frame != Frame::ground) {
    throw std::invalid_argument("sphere rows must be in the ground frame");
  }
  if (!(std::isfinite(radius_m) && radius_m > 0.0)) {
    throw std::invalid_argument("telescope sphere radius must be finite and positive");
  }
  constexpr std::uint32_t sphere_columns =
      column_bit(Column::x) | column_bit(Column::y) | column_bit(Column::z) |
      column_bit(Column::cx) | column_bit(Column::cy) | column_bit(Column::cz) |
      column_bit(Column::t) | column_bit(Column::zem) | column_bit(Column::w);
  if (chunk.rows() > 0 && (chunk.mask() & sphere_columns) != sphere_columns) {
    throw std::invalid_argument("sphere rows lack a required column (x y z cx cy cz t zem w)");
  }
  chunk.validate();
  const double limit = radius_m * (1.0 + sphere_relative_tolerance) + sphere_absolute_tolerance_m;
  for (std::size_t i = 0; i < chunk.rows(); ++i) {
    const double x = chunk[Column::x][i];
    const double y = chunk[Column::y][i];
    const double z = chunk[Column::z][i];
    const double cx = chunk[Column::cx][i];
    const double cy = chunk[Column::cy][i];
    const double cz = chunk[Column::cz][i];
    const double t = chunk[Column::t][i];
    const double zem = chunk[Column::zem][i];
    const double w = chunk[Column::w][i];
    if (!(std::isfinite(x) && std::isfinite(y) && std::isfinite(z) && std::isfinite(cx) &&
          std::isfinite(cy) && std::isfinite(t) && std::isfinite(zem) && std::isfinite(w))) {
      throw std::invalid_argument("sphere row has a non-finite value");
    }
    if (!cz_is_consistent(cz, cx, cy)) {
      throw std::invalid_argument("sphere row direction (cx, cy, cz) is not a unit vector");
    }
    if (std::sqrt(x * x + y * y + z * z) > limit) {
      throw std::invalid_argument("sphere row lies outside the telescope sphere");
    }
  }
}

// One decoded record. `metadata` is used by all metadata records, `chunk` by
// photons, `blob_*` by blob records.
struct Record {
  RecordType type{RecordType::stream};
  Metadata metadata{};
  PhotonChunk chunk{};
  std::uint32_t blob_eventio_type{};
  std::uint32_t blob_eventio_id{};
  std::vector<std::uint8_t> blob{};
};

} // namespace obdeect::stream
