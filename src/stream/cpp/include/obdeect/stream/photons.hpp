#pragma once

#include "obdeect/stream/format.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace obdeect::stream {

// Metadata value: int, float, string or float array (docs/format.md, "Metadata").
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
  }

private:
  std::size_t rows_{};
};

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
