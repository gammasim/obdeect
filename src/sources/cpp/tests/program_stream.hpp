#pragma once

// Reads back a stream written by dlamb::programs with gstrickd::Reader (test helper).

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/stream.hpp"

#include <cstddef>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct Row {
  float x{}, y{}, z{}, cx{}, cy{}, cz{}, t{}, w{}, wavelength{}, zem{};
  bool has_zem{};
  std::uint32_t frame{};
  std::uint32_t telescope{};
  std::uint64_t index{};
  std::int64_t event{};
};

struct StreamContents {
  std::vector<gstrickd::RecordType> types; // every record in order
  std::vector<std::pair<gstrickd::RecordType, gstrickd::Metadata>> metadata;
  std::vector<Row> rows; // all photons in order
};

inline StreamContents read_stream(const std::string &bytes) {
  std::istringstream in(bytes);
  gstrickd::Reader reader(in);
  gstrickd::Record record;
  StreamContents contents;
  std::int64_t event = 0;
  while (reader.next(record)) {
    contents.types.push_back(record.type);
    if (record.type == gstrickd::RecordType::photons) {
      const auto &chunk = record.chunk;
      const bool has_zem = chunk.has(gstrickd::Column::zem);
      for (std::size_t i = 0; i < chunk.rows(); ++i) {
        Row row;
        row.x = chunk[gstrickd::Column::x][i];
        row.y = chunk[gstrickd::Column::y][i];
        row.z = chunk.has(gstrickd::Column::z) ? chunk[gstrickd::Column::z][i] : 0.0F;
        row.cx = chunk[gstrickd::Column::cx][i];
        row.cy = chunk[gstrickd::Column::cy][i];
        row.cz = chunk[gstrickd::Column::cz][i];
        row.t = chunk[gstrickd::Column::t][i];
        row.w = chunk[gstrickd::Column::w][i];
        row.wavelength =
            chunk.has(gstrickd::Column::wavelength) ? chunk[gstrickd::Column::wavelength][i] : 0.0F;
        row.has_zem = has_zem;
        row.zem = has_zem ? chunk[gstrickd::Column::zem][i] : 0.0F;
        row.frame = static_cast<std::uint32_t>(chunk.frame);
        row.telescope = chunk.telescope_id;
        row.index = chunk.first_index + i;
        row.event = event;
        contents.rows.push_back(row);
      }
    } else if (gstrickd::is_metadata_record(record.type)) {
      if (record.type == gstrickd::RecordType::event) {
        event = std::get<std::int64_t>(record.metadata.at("event_number"));
      }
      contents.metadata.emplace_back(record.type, record.metadata);
    }
  }
  if (!reader.ended()) {
    throw std::runtime_error("stream without stream_end");
  }
  return contents;
}

// Number of records of one type.
inline std::size_t count_type(const StreamContents &c, gstrickd::RecordType type) {
  std::size_t n = 0;
  for (const auto t : c.types) {
    n += t == type;
  }
  return n;
}
