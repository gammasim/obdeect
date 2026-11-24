#pragma once

// Shared fixtures of the propagation tests: reference-data paths, photon-stream writing and
// reading, the event and telescope metadata of a run, and the non-finite constants. Test-only;
// the run keys stay in each test because the geometry of every test differs.

#include "check.hpp"

#include "obdeect/stream/stream.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <functional>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace obdeect_test {

inline constexpr double kInf = std::numeric_limits<double>::infinity();
inline constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

// Path of a file in cpp/tests/data.
inline std::string data_file(const std::string &name) {
  return std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/" + name;
}

// A stream of the given records, written by body(writer) and finished.
inline std::string make_stream(const std::function<void(obdeect::stream::Writer &)> &body) {
  std::ostringstream out;
  {
    obdeect::stream::Writer writer(out, "test");
    body(writer);
    writer.finish();
  }
  return out.str();
}

// All records of a stream; the reader must reach the stream_end record.
inline std::vector<obdeect::stream::Record> records_of(const std::string &bytes) {
  std::istringstream in(bytes);
  obdeect::stream::Reader reader(in);
  std::vector<obdeect::stream::Record> records;
  obdeect::stream::Record record;
  while (reader.next(record)) {
    records.push_back(record);
  }
  CHECK(reader.ended());
  return records;
}

// Values of column c over all photon records of one telescope, in stream order.
inline std::vector<float> photon_column(const std::vector<obdeect::stream::Record> &records,
                                        std::uint32_t telescope, obdeect::stream::Column c) {
  std::vector<float> values;
  for (const auto &record : records) {
    if (record.type == obdeect::stream::RecordType::photons &&
        record.chunk.telescope_id == telescope && record.chunk.has(c)) {
      values.insert(values.end(), record.chunk[c].begin(), record.chunk[c].end());
    }
  }
  return values;
}

// One record as raw bytes (header and payload). The Writer refuses records out of order, so the
// streams that break the nesting are written with these.
inline std::string raw_record(obdeect::stream::RecordType type, const std::string &payload) {
  std::string out;
  const auto put = [&out](std::uint64_t value, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
      out.push_back(static_cast<char>(value >> (8 * i)));
    }
  };
  put(obdeect::stream::record_magic, 4);
  put(static_cast<std::uint16_t>(type), 2);
  put(0, 2);
  put(payload.size(), 8);
  return out + payload;
}

inline std::string raw_metadata(obdeect::stream::RecordType type,
                                const obdeect::stream::Metadata &metadata) {
  const auto payload = obdeect::stream::encode_metadata(metadata);
  return raw_record(type, std::string(payload.begin(), payload.end()));
}

inline obdeect::stream::Metadata event_keys() {
  obdeect::stream::Metadata m;
  m["event_number"] = std::int64_t{7};
  return m;
}

inline obdeect::stream::Metadata telescope_keys(std::uint32_t id, obdeect::stream::Frame frame) {
  obdeect::stream::Metadata m;
  m["telescope_id"] = static_cast<std::int64_t>(id);
  m["frame"] = static_cast<std::int64_t>(frame);
  return m;
}

// The telescope_end record of a telescope that has no further photons.
inline obdeect::stream::Metadata telescope_end_keys(std::uint32_t id) {
  obdeect::stream::Metadata m;
  m["telescope_id"] = static_cast<std::int64_t>(id);
  m["photons"] = 0.0;
  return m;
}

// True if action() runs without throwing; a thrown exception is reported and gives false.
template <class Action> bool accepts(Action &&action) {
  try {
    action();
    return true;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "accepts: unexpected exception: %s\n", error.what());
  }
  return false;
}

} // namespace obdeect_test
