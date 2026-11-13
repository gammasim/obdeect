#pragma once

// Helpers for the sky-level program tests: write a stream to memory, read it back with
// obdeect::stream::Reader, and collect the record types, photons per telescope and chunks.

#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <cstdint>
#include <functional>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

struct StreamSummary {
  std::vector<std::uint16_t> types;
  std::map<std::uint32_t, std::uint64_t> photons; // rows per telescope id
  std::vector<obdeect::stream::PhotonChunk> chunks;
  std::uint64_t events{};
  obdeect::stream::Metadata run;
  bool ended{};
};

inline std::string write_to_string(const std::function<void(std::ostream &)> &write) {
  std::ostringstream out(std::ios::binary);
  write(out);
  return out.str();
}

inline StreamSummary read_stream(const std::string &bytes) {
  std::istringstream in(bytes, std::ios::binary);
  obdeect::stream::Reader reader(in);
  obdeect::stream::Record record;
  StreamSummary s;
  while (reader.next(record)) {
    s.types.push_back(static_cast<std::uint16_t>(record.type));
    if (record.type == obdeect::stream::RecordType::run)
      s.run = record.metadata;
    if (record.type == obdeect::stream::RecordType::event)
      ++s.events;
    if (record.type == obdeect::stream::RecordType::photons) {
      s.photons[record.chunk.telescope_id] += record.chunk.rows();
      s.chunks.push_back(record.chunk);
    }
  }
  s.ended = reader.ended();
  return s;
}
