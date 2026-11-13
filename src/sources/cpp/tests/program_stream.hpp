#pragma once

// Reads back a stream written by obdeect::sources::programs with obdeect::stream::Reader (test
// helper).

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

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
  std::vector<obdeect::stream::RecordType> types; // every record in order
  std::vector<std::pair<obdeect::stream::RecordType, obdeect::stream::Metadata>> metadata;
  std::vector<Row> rows; // all photons in order
};

inline StreamContents read_stream(const std::string &bytes) {
  std::istringstream in(bytes);
  obdeect::stream::Reader reader(in);
  obdeect::stream::Record record;
  StreamContents contents;
  std::int64_t event = 0;
  while (reader.next(record)) {
    contents.types.push_back(record.type);
    if (record.type == obdeect::stream::RecordType::photons) {
      const auto &chunk = record.chunk;
      const bool has_zem = chunk.has(obdeect::stream::Column::zem);
      for (std::size_t i = 0; i < chunk.rows(); ++i) {
        Row row;
        row.x = chunk[obdeect::stream::Column::x][i];
        row.y = chunk[obdeect::stream::Column::y][i];
        row.z = chunk.has(obdeect::stream::Column::z) ? chunk[obdeect::stream::Column::z][i] : 0.0F;
        row.cx = chunk[obdeect::stream::Column::cx][i];
        row.cy = chunk[obdeect::stream::Column::cy][i];
        row.cz = chunk[obdeect::stream::Column::cz][i];
        row.t = chunk[obdeect::stream::Column::t][i];
        row.w = chunk[obdeect::stream::Column::w][i];
        row.wavelength = chunk.has(obdeect::stream::Column::wavelength)
                             ? chunk[obdeect::stream::Column::wavelength][i]
                             : 0.0F;
        row.has_zem = has_zem;
        row.zem = has_zem ? chunk[obdeect::stream::Column::zem][i] : 0.0F;
        row.frame = static_cast<std::uint32_t>(chunk.frame);
        row.telescope = chunk.telescope_id;
        row.index = chunk.first_index + i;
        row.event = event;
        contents.rows.push_back(row);
      }
    } else if (obdeect::stream::is_metadata_record(record.type)) {
      if (record.type == obdeect::stream::RecordType::event) {
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
inline std::size_t count_type(const StreamContents &c, obdeect::stream::RecordType type) {
  std::size_t n = 0;
  for (const auto t : c.types) {
    n += t == type;
  }
  return n;
}
