#include "check.hpp"

#include "obdeect/stream/stream.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

using namespace obdeect::stream;

namespace {

template <class Fn> bool throws(Fn fn) {
  try {
    fn();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

template <class Fn> bool throws_invalid_argument(Fn fn) {
  try {
    fn();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

// True if fn throws an exception of type E (and nothing else).
template <class E, class Fn> bool throws_as(Fn fn) {
  try {
    fn();
  } catch (const E &) {
    return true;
  } catch (...) {
    return false;
  }
  return false;
}

// Little-endian bytes of the low n bytes of value.
std::string le_bytes(std::uint64_t value, std::size_t n) {
  std::string out;
  for (std::size_t i = 0; i < n; ++i) {
    out.push_back(static_cast<char>(value >> (8 * i)));
  }
  return out;
}

// Little-endian float32 bytes of the given value.
std::string float_bytes(float value) { return le_bytes(std::bit_cast<std::uint32_t>(value), 4); }

// One complete record with the given type and payload, little-endian header.
std::string raw_record(std::uint16_t type, const std::string &payload) {
  std::string out =
      le_bytes(record_magic, 4) + le_bytes(type, 2) + le_bytes(0, 2) + le_bytes(payload.size(), 8);
  return out + payload;
}

std::string as_string(const std::vector<std::uint8_t> &bytes) {
  return {bytes.begin(), bytes.end()};
}

// A metadata record of the given type (written raw, so that the writer's checks are not used).
std::string raw_metadata(RecordType type, const Metadata &metadata) {
  return raw_record(static_cast<std::uint16_t>(type), as_string(encode_metadata(metadata)));
}

// Stream header record of the given format_minor (major 1).
std::string stream_header(std::int64_t minor) {
  const Metadata header{{"format_major", std::int64_t{format_major}},
                        {"format_minor", minor},
                        {"producer", std::string("test")}};
  return raw_metadata(RecordType::stream, header);
}

// The stream end record.
std::string stream_end() {
  return raw_record(static_cast<std::uint16_t>(RecordType::stream_end), "");
}

// Header, records and stream end, as raw bytes.
std::string wrap(std::int64_t minor, const std::string &records) {
  return stream_header(minor) + records + stream_end();
}

// Photons payload of a chunk, written raw: the reader's checks are the subject of the test.
std::string chunk_payload(const PhotonChunk &chunk) {
  std::string out =
      le_bytes(chunk.telescope_id, 4) + le_bytes(static_cast<std::uint32_t>(chunk.frame), 4) +
      le_bytes(chunk.mask(), 4) + le_bytes(chunk.rows(), 4) + le_bytes(chunk.first_index, 8);
  for (std::uint32_t i = 0; i < n_columns; ++i) {
    for (const float value : chunk.columns[i]) {
      out += float_bytes(value);
    }
  }
  return out;
}

// A block of one telescope (id 1, ground frame) inside one event of one run, around `inner`.
std::string block(const std::string &inner) {
  return raw_metadata(RecordType::run, {}) + raw_metadata(RecordType::event, {}) +
         raw_metadata(RecordType::telescope,
                      {{"telescope_id", std::int64_t{1}}, {"frame", std::int64_t{1}}}) +
         inner +
         raw_metadata(RecordType::telescope_end,
                      {{"telescope_id", std::int64_t{1}}, {"photons", 0.0}}) +
         raw_metadata(RecordType::event_end, {}) + raw_metadata(RecordType::run_end, {});
}

// Photons record of a chunk (raw).
std::string photons_of(const PhotonChunk &chunk) {
  return raw_record(static_cast<std::uint16_t>(RecordType::photons), chunk_payload(chunk));
}

// Efficiency status record of (factor, state) pairs (raw).
std::string status_of(std::initializer_list<std::pair<std::string, std::string>> factors) {
  Metadata metadata;
  std::string names;
  for (const auto &[name, state] : factors) {
    names += (names.empty() ? "" : ",") + name;
    metadata[name + ".state"] = state;
  }
  metadata["factors"] = names;
  return raw_metadata(RecordType::efficiency_status, metadata);
}

// The metadata of an efficiency_status record (for the Writer).
Metadata status_of_metadata(std::initializer_list<std::pair<std::string, std::string>> factors) {
  Metadata metadata;
  std::string names;
  for (const auto &[name, state] : factors) {
    names += (names.empty() ? "" : ",") + name;
    metadata[name + ".state"] = state;
  }
  metadata["factors"] = names;
  return metadata;
}

PhotonChunk make_chunk(std::uint32_t telescope, Frame frame, std::uint64_t first, std::size_t rows,
                       std::uint32_t extra_columns) {
  PhotonChunk chunk;
  chunk.telescope_id = telescope;
  chunk.frame = frame;
  chunk.first_index = first;
  chunk.resize(rows, required_columns | extra_columns);
  for (std::uint32_t i = 0; i < n_columns; ++i) {
    if (chunk.columns[i].empty()) {
      continue;
    }
    for (std::size_t r = 0; r < rows; ++r) {
      chunk.columns[i][r] = static_cast<float>(i) * 1000.0F + static_cast<float>(r) * 0.25F;
    }
  }
  // Directions must be unit vectors (the reader checks them): small cx, cy and cz from them.
  for (std::size_t r = 0; r < rows; ++r) {
    chunk.columns[static_cast<std::size_t>(Column::cx)][r] = 0.001F * static_cast<float>(r + 1);
    chunk.columns[static_cast<std::size_t>(Column::cy)][r] = 0.002F * static_cast<float>(r);
    if (!chunk.columns[static_cast<std::size_t>(Column::cz)].empty()) {
      const double cx = chunk.columns[static_cast<std::size_t>(Column::cx)][r];
      const double cy = chunk.columns[static_cast<std::size_t>(Column::cy)][r];
      chunk.columns[static_cast<std::size_t>(Column::cz)][r] =
          static_cast<float>(-std::sqrt(1.0 - cx * cx - cy * cy));
    }
  }
  return chunk;
}

// Writes run, event, one telescope block per chunk (telescope id and frame from the chunk), then
// the closing records. The Writer checks the nesting but not the photon values, so this is the way
// to make a valid stream from chunks.
std::string blocks_stream(const Metadata &run, const std::vector<PhotonChunk> &chunks) {
  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.write_metadata(RecordType::run, run);
    writer.write_metadata(RecordType::event, Metadata{{"event_number", std::int64_t{1}}});
    for (const PhotonChunk &chunk : chunks) {
      const auto id = std::int64_t{chunk.telescope_id};
      writer.write_metadata(
          RecordType::telescope,
          Metadata{{"telescope_id", id}, {"frame", static_cast<std::int64_t>(chunk.frame)}});
      writer.write_chunk(chunk);
      writer.write_metadata(RecordType::telescope_end,
                            Metadata{{"telescope_id", id}, {"photons", 0.0}});
    }
    writer.write_metadata(RecordType::event_end, Metadata{});
    writer.write_metadata(RecordType::run_end, Metadata{});
    writer.finish();
  }
  return out.str();
}

// Reads all records; returns their types and keeps the chunks in `chunks`.
std::vector<RecordType> read_all(const std::string &bytes, std::vector<PhotonChunk> *chunks) {
  std::istringstream in(bytes);
  Reader reader(in);
  Record record;
  std::vector<RecordType> types;
  while (reader.next(record)) {
    types.push_back(record.type);
    if (record.type == RecordType::photons && chunks) {
      chunks->push_back(record.chunk);
    }
  }
  CHECK(reader.ended());
  return types;
}

void metadata_round_trip() {
  Metadata m;
  m["int_negative"] = std::int64_t{-42};
  m["int_big"] = std::int64_t{1} << 40;
  m["double"] = 3.5;
  m["double_tiny"] = 1e-300;
  m["text"] = std::string("hello, world");
  m["text_empty"] = std::string();
  m["array"] = std::vector<double>{1.0, 2.5, -3.0};
  m["array_empty"] = std::vector<double>{};
  const auto bytes = encode_metadata(m);
  CHECK(decode_metadata(bytes.data(), bytes.size()) == m);
}

void metadata_layout() {
  Metadata m;
  m["a"] = std::int64_t{1};
  const std::vector<std::uint8_t> expected = {1, 0, 0,   0,              // count
                                              1, 0, 'a',                 // key
                                              1,                         // int
                                              1, 0, 0,   0, 0, 0, 0, 0}; // value
  CHECK(encode_metadata(m) == expected);
}

void metadata_malformed() {
  Metadata m;
  m["a"] = std::int64_t{1};
  const auto bytes = encode_metadata(m);
  CHECK(throws([&] { static_cast<void>(decode_metadata(bytes.data(), bytes.size() - 1)); }));
  auto unknown_type = bytes;
  unknown_type[7] = 9;
  CHECK(throws(
      [&] { static_cast<void>(decode_metadata(unknown_type.data(), unknown_type.size())); }));
  auto trailing = bytes;
  trailing.push_back(0);
  CHECK(throws([&] { static_cast<void>(decode_metadata(trailing.data(), trailing.size())); }));
}

void stream_round_trip() {
  Metadata run{{"run_number", std::int64_t{7}}, {"telescope_x_m", std::vector<double>{1.0, 2.0}}};
  Metadata event{{"event_number", std::int64_t{1}}, {"energy_gev", 0.5}};
  Metadata telescope{{"telescope_id", std::int64_t{3}}, {"frame", std::int64_t{2}}};
  Metadata efficiency_status{{"factors", std::string("qe,mirror")},
                             {"qe.state", std::string("applied")},
                             {"mirror.state", std::string("pending")}};
  Metadata end{{"photons", 12.5}};

  const PhotonChunk first = make_chunk(3, Frame::ground, 0, 100, 0);
  const PhotonChunk second = make_chunk(3, Frame::telescope, 100, 50,
                                        column_bit(Column::zem) | column_bit(Column::wavelength));

  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.write_metadata(RecordType::run, run);
    writer.write_metadata(RecordType::efficiency_status, efficiency_status);
    writer.write_blob(1204, 0, {1, 2, 3, 4, 5, 6, 7, 8, 9});
    writer.write_metadata(RecordType::event, event);
    writer.write_metadata(RecordType::telescope, telescope);
    writer.write_chunk(first);
    writer.write_chunk(second);
    writer.write_metadata(RecordType::telescope_end, end);
    writer.write_metadata(RecordType::event_end, Metadata{});
    writer.write_metadata(RecordType::run_end, Metadata{});
    writer.finish();
    writer.finish();
  }

  std::vector<PhotonChunk> chunks;
  const auto types = read_all(out.str(), &chunks);
  const std::vector<RecordType> expected = {
      RecordType::stream,    RecordType::run,     RecordType::efficiency_status,
      RecordType::blob,      RecordType::event,   RecordType::telescope,
      RecordType::photons,   RecordType::photons, RecordType::telescope_end,
      RecordType::event_end, RecordType::run_end};
  CHECK(types == expected);
  CHECK(chunks.size() == 2);
  if (chunks.size() == 2) {
    CHECK(chunks[0].telescope_id == 3);
    CHECK(chunks[0].frame == Frame::ground);
    CHECK(chunks[0].first_index == 0);
    CHECK(chunks[0].rows() == 100);
    CHECK(chunks[0].columns == first.columns);
    CHECK(chunks[1].frame == Frame::telescope);
    CHECK(chunks[1].first_index == 100);
    CHECK(chunks[1].rows() == 50);
    CHECK(chunks[1].columns == second.columns);
    CHECK(chunks[1].has(Column::zem) && chunks[1].has(Column::wavelength));
    CHECK(!chunks[1].has(Column::p));
  }

  // Metadata and blob content survive the reader.
  std::istringstream in(out.str());
  Reader reader(in);
  Record record;
  CHECK(reader.next(record) && record.type == RecordType::stream);
  CHECK(record.metadata.at("format_major") == Value{std::int64_t{format_major}});
  CHECK(record.metadata.at("producer") == Value{std::string("test")});
  CHECK(reader.next(record) && record.metadata == run);
  CHECK(reader.next(record) && record.metadata == efficiency_status);
  CHECK(reader.next(record) && record.type == RecordType::blob);
  CHECK(record.blob_eventio_type == 1204 && record.blob_eventio_id == 0);
  CHECK(record.blob == (std::vector<std::uint8_t>{1, 2, 3, 4, 5, 6, 7, 8, 9}));
}

void unknown_record_skipped() {
  std::ostringstream out;
  Writer writer(out, "test");
  writer.write_metadata(RecordType::run, Metadata{{"run_number", std::int64_t{1}}});
  out << raw_record(12, "future record payload");
  writer.write_metadata(RecordType::run_end, Metadata{});
  writer.finish();
  const auto types = read_all(out.str(), nullptr);
  const std::vector<RecordType> expected = {RecordType::stream, RecordType::run,
                                            RecordType::run_end};
  CHECK(types == expected);
}

void empty_chunk() {
  const PhotonChunk empty = make_chunk(5, Frame::camera, 42, 0, 0);
  std::vector<PhotonChunk> chunks;
  CHECK(read_all(blocks_stream(Metadata{}, {empty}), &chunks).size() == 8);
  CHECK(chunks.size() == 1);
  if (chunks.size() == 1) {
    CHECK(chunks[0].rows() == 0);
    CHECK(chunks[0].mask() == 0);
    CHECK(chunks[0].telescope_id == 5);
    CHECK(chunks[0].frame == Frame::camera);
    CHECK(chunks[0].first_index == 42);
  }
}

void invalid_chunk_rejected() {
  PhotonChunk chunk;
  chunk.resize(2, column_bit(Column::x));
  std::ostringstream out;
  Writer writer(out, "test");
  CHECK(throws([&] { writer.write_chunk(chunk); }));
  CHECK(throws([&] { writer.write_metadata(RecordType::photons, Metadata{}); }));
}

// Photons payload with its fixed header and no column data.
std::string photons_header(std::uint32_t telescope, std::uint32_t frame, std::uint32_t mask,
                           std::uint32_t rows) {
  return le_bytes(telescope, 4) + le_bytes(frame, 4) + le_bytes(mask, 4) + le_bytes(rows, 4) +
         le_bytes(0, 8);
}

void filtered_photons() {
  const std::string bytes = blocks_stream(Metadata{}, {make_chunk(1, Frame::sky, 7, 10, 0),
                                                       make_chunk(2, Frame::telescope, 100, 5, 0),
                                                       make_chunk(3, Frame::camera, 555, 4, 0)});
  std::istringstream in(bytes);
  Reader reader(in);
  reader.keep_photons_of({2});
  Record record;
  std::vector<PhotonChunk> seen;
  CHECK(reader.next(record) && record.type == RecordType::stream);
  while (reader.next(record)) {
    if (record.type == RecordType::photons) {
      seen.push_back(record.chunk); // copies, so a skipped chunk after a kept one is checked too
    }
  }
  CHECK(reader.ended());
  CHECK(seen.size() == 3);
  if (seen.size() == 3) {
    // Skipped chunks report their header fields and the header row count, but no columns.
    CHECK(seen[0].telescope_id == 1 && seen[0].frame == Frame::sky);
    CHECK(seen[0].first_index == 7 && seen[0].rows() == 10 && seen[0].mask() == 0);
    CHECK(seen[1].telescope_id == 2 && seen[1].frame == Frame::telescope);
    CHECK(seen[1].first_index == 100 && seen[1].rows() == 5);
    CHECK(seen[1].columns == make_chunk(2, Frame::telescope, 100, 5, 0).columns);
    CHECK(seen[2].telescope_id == 3 && seen[2].frame == Frame::camera);
    CHECK(seen[2].first_index == 555 && seen[2].rows() == 4 && seen[2].mask() == 0);
  }
}

void filtered_photons_validated() {
  // Unselected chunks still have their header and payload size checked: unknown frame or wrong
  // payload size is rejected.
  const auto read_keeping_two = [](const std::string &bytes) {
    std::istringstream in(bytes);
    Reader reader(in);
    reader.keep_photons_of({2});
    Record record;
    while (reader.next(record)) {
    }
  };
  const std::uint32_t photons = static_cast<std::uint32_t>(RecordType::photons);
  for (const std::string &payload :
       {photons_header(9, 99, 0, 0), photons_header(9, 0, required_columns, 3)}) {
    CHECK(throws([&] {
      read_keeping_two(
          wrap(format_minor, block(raw_record(static_cast<std::uint16_t>(photons), payload))));
    }));
  }
}

void bad_magic() {
  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.finish();
  }
  std::string bytes = out.str();
  bytes[0] = 'X';
  CHECK(throws([&] { read_all(bytes, nullptr); }));
}

void truncated_record() {
  const std::string bytes = blocks_stream(Metadata{}, {make_chunk(1, Frame::sky, 0, 10, 0)});
  CHECK(throws([&] { read_all(bytes.substr(0, bytes.size() - 5), nullptr); }));
  CHECK(throws([&] { read_all(bytes.substr(0, bytes.size() - 20), nullptr); }));
  CHECK(throws([&] { read_all(bytes.substr(0, 10), nullptr); }));
}

// Reads a stream to its end; returns true if stream_end was read. Throws the reader's errors.
bool reads_complete(const std::string &bytes) {
  std::istringstream in(bytes);
  Reader reader(in);
  Record record;
  while (reader.next(record)) {
  }
  return reader.ended();
}

// Every proper prefix of a valid stream is an incomplete stream: the reader either throws its own
// error (a cut record) or reports that stream_end was not read. No prefix is accepted as a complete
// run, and none crashes the reader.
void every_prefix_is_incomplete() {
  const std::string valid = blocks_stream(Metadata{{"run_number", std::int64_t{1}}},
                                          {make_chunk(1, Frame::sky, 0, 10, 0)});
  CHECK(reads_complete(valid));
  for (std::size_t length = 0; length < valid.size(); ++length) {
    const std::string prefix = valid.substr(0, length);
    try {
      CHECK(!reads_complete(prefix));
    } catch (const std::runtime_error &error) {
      CHECK(std::string_view(error.what()).find("obdeect stream:") != std::string_view::npos);
    }
  }
}

// Flipping any byte of the header region (magic, stream record) either makes the reader throw
// its own error or still reads a well-formed stream; it never crashes or throws another type.
void header_byte_flips_are_contained() {
  const std::string valid = blocks_stream(Metadata{{"run_number", std::int64_t{1}}},
                                          {make_chunk(1, Frame::sky, 0, 10, 0)});
  constexpr std::size_t header_region = 64;
  CHECK(valid.size() > header_region);
  for (std::size_t i = 0; i < header_region; ++i) {
    for (const unsigned char mask : {0x01U, 0x80U, 0xFFU}) {
      std::string flipped = valid;
      flipped[i] = static_cast<char>(static_cast<unsigned char>(flipped[i]) ^ mask);
      try {
        read_all(flipped, nullptr);
      } catch (const std::runtime_error &error) {
        CHECK(std::string_view(error.what()).find("obdeect stream:") != std::string_view::npos);
      } catch (const std::exception &error) {
        std::fprintf(stderr, "byte %zu flip: unexpected exception: %s\n", i, error.what());
        CHECK(false);
      }
    }
  }
}

void wrong_format_major() {
  Metadata header{{"format_major", std::int64_t{format_major + 1}},
                  {"format_minor", std::int64_t{0}},
                  {"producer", std::string("future")}};
  const std::string bytes = raw_metadata(RecordType::stream, header);
  CHECK(throws([&] { read_all(bytes, nullptr); }));
}

// Readers accept every format_minor of their major version: 1.0 streams (no upward photons)
// and the current 1.1.
void format_minor_accepted() {
  for (const std::int64_t minor : {std::int64_t{0}, std::int64_t{format_minor}}) {
    CHECK(!throws([&] { read_all(wrap(minor, ""), nullptr); }));
  }
}

void header_required_once() {
  // Records before the stream header, and a second header, are malformed.
  CHECK(throws([&] { read_all(stream_end(), nullptr); }));
  CHECK(throws([&] { read_all(raw_metadata(RecordType::run, {}), nullptr); }));
  CHECK(throws([&] { read_all(wrap(format_minor, stream_header(format_minor)), nullptr); }));
}

void clean_end_required() {
  // Bytes after stream_end (here a whole extra record) are not a clean end.
  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.finish();
  }
  const std::string clean = out.str();
  for (const std::string &bytes : {clean, clean + raw_record(12, "trailing")}) {
    std::istringstream in(bytes);
    Reader reader(in);
    Record record;
    while (reader.next(record)) {
    }
    CHECK(reader.ended());
    CHECK(reader.at_clean_end() == (bytes == clean));
  }
}

// Photons record of telescope 1 in the ground frame: every column of the mask holds value + index,
// and the direction columns hold small values (a unit direction with cz absent).
std::string raw_photons(std::uint32_t mask, std::uint32_t rows, float value) {
  std::string payload = photons_header(1, static_cast<std::uint32_t>(Frame::ground), mask, rows);
  for (std::uint32_t i = 0; i < n_columns; ++i) {
    if (((mask >> i) & 1U) == 0) {
      continue;
    }
    const bool direction =
        i == static_cast<std::uint32_t>(Column::cx) || i == static_cast<std::uint32_t>(Column::cy);
    for (std::uint32_t r = 0; r < rows; ++r) {
      payload += float_bytes(direction ? 0.0F : value + static_cast<float>(i));
    }
  }
  return raw_record(static_cast<std::uint16_t>(RecordType::photons), payload);
}

// Path column: every frame round-trips bit-exactly, with and without the other optional columns.
void path_round_trip() {
  CHECK(column_name(Column::path) == std::string("path"));
  CHECK(column_bit(Column::path) == (1U << 11));
  for (const Frame frame : {Frame::sky, Frame::ground, Frame::telescope, Frame::camera}) {
    PhotonChunk chunk = make_chunk(4, frame, 9, 64, column_bit(Column::path));
    chunk[Column::path][0] = 0.0F;
    chunk[Column::path][1] = 1.0000001F;
    chunk[Column::path][2] = 3.4e38F;
    chunk[Column::path][3] = 1e-38F;
    std::vector<PhotonChunk> chunks;
    read_all(blocks_stream(Metadata{}, {chunk}), &chunks);
    CHECK(chunks.size() == 1);
    if (chunks.size() == 1) {
      CHECK(chunks[0].has(Column::path));
      CHECK(chunks[0].columns == chunk.columns);
      CHECK(chunks[0].mask() == (required_columns | column_bit(Column::path)));
    }
  }
  // Negative or non-finite path values are rejected by the writer and by the reader.
  for (const float bad : {-1.0F, std::nanf(""), INFINITY}) {
    PhotonChunk chunk = make_chunk(1, Frame::ground, 0, 2, column_bit(Column::path));
    chunk[Column::path][1] = bad;
    std::ostringstream out;
    Writer writer(out, "test");
    CHECK(throws([&] { writer.write_chunk(chunk); }));
    // A photons record of two rows whose path values are bad (the other columns are valid).
    const std::uint32_t mask = required_columns | column_bit(Column::path);
    std::string payload = photons_header(1, static_cast<std::uint32_t>(Frame::ground), mask, 2);
    for (std::uint32_t i = 0; i < n_columns; ++i) {
      for (std::uint32_t r = 0; r < 2 && ((mask >> i) & 1U); ++r) {
        const bool direction = i == static_cast<std::uint32_t>(Column::cx) ||
                               i == static_cast<std::uint32_t>(Column::cy);
        payload += float_bytes(i == static_cast<std::uint32_t>(Column::path) ? bad
                               : direction                                   ? 0.0F
                                                                             : 1.0F);
      }
    }
    const std::string bytes = wrap(
        format_minor, block(raw_record(static_cast<std::uint16_t>(RecordType::photons), payload)));
    CHECK(throws([&] { read_all(bytes, nullptr); }));
  }
}

// Streams of format 1.1 (no path column, no sphere runs) are read by the 1.2 reader; a path column
// or a sphere run in them is malformed. A later minor version is read (minor versions are
// compatible within a major version); a negative one is rejected.
void format_minor_columns() {
  const std::uint32_t plain =
      required_columns | column_bit(Column::zem) | column_bit(Column::wavelength);
  std::vector<PhotonChunk> chunks;
  CHECK(!throws([&] { read_all(wrap(1, block(raw_photons(plain, 3, 1.0F))), &chunks); }));
  CHECK(chunks.size() == 1);
  if (chunks.size() == 1) {
    CHECK(chunks[0].rows() == 3 && !chunks[0].has(Column::path));
    CHECK(chunks[0].has(Column::zem) && chunks[0].has(Column::wavelength));
    CHECK_NEAR(chunks[0][Column::x][2], 1.0, 0.0); // raw_photons: value 1.0 + column index 0
  }
  CHECK(throws([&] {
    read_all(wrap(1, block(raw_photons(plain | column_bit(Column::path), 3, 1.0F))), nullptr);
  }));
  const Metadata sphere{{"ground_representation", std::string("sphere")}};
  const std::string sphere_run =
      raw_metadata(RecordType::run, sphere) + raw_metadata(RecordType::run_end, {});
  CHECK(throws([&] { read_all(wrap(1, sphere_run), nullptr); }));
  CHECK(!throws([&] { read_all(wrap(2, sphere_run), nullptr); }));
  CHECK(!throws([&] { read_all(wrap(3, sphere_run), nullptr); }));
  CHECK(throws([&] { read_all(wrap(-1, ""), nullptr); }));
}

void ground_representation_helper() {
  CHECK(ground_representation(Metadata{}) == GroundRepresentation::plane);
  CHECK(ground_representation(Metadata{{"run_number", std::int64_t{1}}}) ==
        GroundRepresentation::plane);
  CHECK(ground_representation(Metadata{{"ground_representation", std::string("plane")}}) ==
        GroundRepresentation::plane);
  CHECK(ground_representation(Metadata{{"ground_representation", std::string("sphere")}}) ==
        GroundRepresentation::sphere);
  for (const Value &bad :
       {Value{std::string("Sphere")}, Value{std::string("")}, Value{std::int64_t{1}}, Value{2.5}}) {
    CHECK(throws([&] {
      static_cast<void>(ground_representation(Metadata{{"ground_representation", bad}}));
    }));
  }
  // The reader rejects an invalid value in a run record; the writer never writes one.
  CHECK(throws([&] {
    read_all(wrap(2, raw_metadata(RecordType::run,
                                  Metadata{{"ground_representation", std::string("flat")}})),
             nullptr);
  }));
  std::ostringstream out;
  Writer writer(out, "test");
  CHECK(throws([&] {
    writer.write_metadata(RecordType::run,
                          Metadata{{"ground_representation", std::string("flat")}});
  }));
  CHECK(!throws([&] {
    writer.write_metadata(RecordType::run,
                          Metadata{{"ground_representation", std::string("sphere")}});
  }));
}

// A ground-frame chunk with x y z cx cy cz t zem w (and the optional columns given in extra).
PhotonChunk sphere_rows(std::size_t rows, std::uint32_t extra = 0) {
  PhotonChunk chunk;
  chunk.telescope_id = 1;
  chunk.frame = Frame::ground;
  const std::uint32_t columns =
      column_bit(Column::x) | column_bit(Column::y) | column_bit(Column::z) |
      column_bit(Column::cx) | column_bit(Column::cy) | column_bit(Column::cz) |
      column_bit(Column::t) | column_bit(Column::zem) | column_bit(Column::w);
  chunk.resize(rows, columns | extra);
  for (std::size_t i = 0; i < rows; ++i) {
    // Row 0: on the sphere of radius 12 m above the centre (6 m); row 1: inside (3 m).
    const bool on = i == 0;
    chunk[Column::x][i] = on ? 10.392304F : 3.0F;
    chunk[Column::y][i] = 0.0F;
    chunk[Column::z][i] = on ? 6.0F : 0.0F;
    chunk[Column::cx][i] = 0.6F;
    chunk[Column::cy][i] = 0.0F;
    chunk[Column::cz][i] = -0.8F;
    chunk[Column::t][i] = 1.0F;
    chunk[Column::zem][i] = 2200.0F;
    chunk[Column::w][i] = 1.0F;
  }
  return chunk;
}

// validate_sphere_rows: sphere rows on and inside the sphere pass; outside rows, missing columns,
// missing z or cz, inconsistent cz, other frames and non-finite values are rejected.
void sphere_rows_validation() {
  constexpr double radius = 12.0;
  CHECK(!throws([&] { validate_sphere_rows(sphere_rows(1), radius); }));
  CHECK(!throws([&] { validate_sphere_rows(sphere_rows(2), radius); }));
  CHECK(!throws([&] { validate_sphere_rows(sphere_rows(2, column_bit(Column::path)), radius); }));
  // Empty chunks have nothing to check.
  CHECK(!throws([&] { validate_sphere_rows(PhotonChunk{}, radius); }));

  PhotonChunk outside = sphere_rows(1);
  outside[Column::x][0] = 0.0F;
  outside[Column::z][0] = 13.0F; // 13 m from the centre
  outside[Column::cz][0] = -1.0F;
  outside[Column::cx][0] = 0.0F;
  CHECK(throws([&] { validate_sphere_rows(outside, radius); }));
  CHECK(!throws([&] { validate_sphere_rows(outside, 13.5); }));

  for (const Column column : {Column::z, Column::cz, Column::zem, Column::t, Column::w}) {
    PhotonChunk missing = sphere_rows(2);
    missing[column].clear();
    CHECK(throws([&] { validate_sphere_rows(missing, radius); }));
  }
  PhotonChunk plane = sphere_rows(2);
  plane.frame = Frame::telescope;
  CHECK(throws([&] { validate_sphere_rows(plane, radius); }));
  plane.frame = Frame::sky;
  CHECK(throws([&] { validate_sphere_rows(plane, radius); }));

  PhotonChunk bad_cz = sphere_rows(2);
  bad_cz[Column::cz][1] = -0.5F; // cx^2 + cy^2 + cz^2 != 1
  CHECK(throws([&] { validate_sphere_rows(bad_cz, radius); }));
  PhotonChunk nan_x = sphere_rows(2);
  nan_x[Column::x][1] = std::nanf("");
  CHECK(throws([&] { validate_sphere_rows(nan_x, radius); }));
  CHECK(throws([&] { validate_sphere_rows(sphere_rows(1), 0.0); }));
  CHECK(throws([&] { validate_sphere_rows(sphere_rows(1), std::nan("")); }));
}

// Sphere runs in the reader: the run key and the sphere rows round-trip bit-exactly.
void sphere_stream_round_trip() {
  const std::string bytes =
      blocks_stream(Metadata{{"ground_representation", std::string("sphere")}},
                    {sphere_rows(2, column_bit(Column::path))});
  std::vector<PhotonChunk> chunks;
  CHECK(read_all(bytes, &chunks) ==
        (std::vector<RecordType>{RecordType::stream, RecordType::run, RecordType::event,
                                 RecordType::telescope, RecordType::photons,
                                 RecordType::telescope_end, RecordType::event_end,
                                 RecordType::run_end}));
  CHECK(chunks.size() == 1);
  if (chunks.size() == 1) {
    CHECK(chunks[0].columns == sphere_rows(2, column_bit(Column::path)).columns);
    CHECK(chunks[0].frame == Frame::ground);
  }
}

// Records that break the nesting of docs/stream/format.md are rejected by the reader; each stream
// below is otherwise valid.
void nesting_rejected() {
  const std::string run = raw_metadata(RecordType::run, {});
  const std::string event = raw_metadata(RecordType::event, {});
  const std::string event_end = raw_metadata(RecordType::event_end, {});
  const std::string run_end = raw_metadata(RecordType::run_end, {});
  const Metadata telescope_1{{"telescope_id", std::int64_t{1}}, {"frame", std::int64_t{1}}};
  const std::string telescope = raw_metadata(RecordType::telescope, telescope_1);
  const std::string telescope_end =
      raw_metadata(RecordType::telescope_end, {{"telescope_id", std::int64_t{1}}});
  const std::string photons = raw_photons(required_columns, 2, 1.0F);
  const std::string blob =
      raw_record(static_cast<std::uint16_t>(RecordType::blob), std::string(8, '\0') + "data");
  const std::string status = status_of({{"qe", "pending"}});

  // Controls: the nested forms are accepted.
  CHECK(!throws([&] {
    read_all(
        wrap(2, run + blob + event + telescope + photons + telescope_end + event_end + run_end),
        nullptr);
  }));
  CHECK(!throws([&] { read_all(wrap(2, run + run_end), nullptr); }));

  const std::vector<std::pair<const char *, std::string>> violations = {
      {"photons outside a telescope", run + event + photons + event_end + run_end},
      {"telescope_end without telescope", run + event + telescope_end + event_end + run_end},
      {"run_end without run", run + run_end + run_end},
      {"event inside an event", run + event + event + event_end + event_end + run_end},
      {"event outside a run", event + event_end},
      {"telescope outside an event", run + telescope + telescope_end + run_end},
      {"telescope inside a telescope",
       run + event + telescope + telescope + telescope_end + telescope_end + event_end + run_end},
      {"event_end inside a telescope",
       run + event + telescope + event_end + telescope_end + event_end + run_end},
      {"run inside a run", run + run + run_end},
      {"blob outside a run", blob},
      {"stream_end inside a run", run},
      {"run_end inside an event", run + event + run_end},
      {"efficiency status outside a run", status + run + run_end},
      {"run efficiency status after the first event", run + event + event_end + status + run_end},
      {"telescope efficiency status after its photons",
       run + event + telescope + photons + status + telescope_end + event_end + run_end},
  };
  for (const auto &violation : violations) {
    const std::string bytes = stream_header(format_minor) + violation.second + stream_end();
    CHECK(throws([&] { read_all(bytes, nullptr); }));
  }
}

// Efficiency status transitions (docs/stream/efficiency-status.md): pending -> envelope -> applied
// per factor and scope; a repeated or backward state is rejected.
void efficiency_transitions() {
  const std::string run = raw_metadata(RecordType::run, {});
  const std::string run_end = raw_metadata(RecordType::run_end, {});
  const std::string event = raw_metadata(RecordType::event, {});
  const std::string event_end = raw_metadata(RecordType::event_end, {});
  const std::string telescope = raw_metadata(
      RecordType::telescope, {{"telescope_id", std::int64_t{1}}, {"frame", std::int64_t{1}}});
  const std::string telescope_end =
      raw_metadata(RecordType::telescope_end, {{"telescope_id", std::int64_t{1}}});
  const std::string photons = raw_photons(required_columns, 2, 1.0F);

  const auto run_with = [&](std::initializer_list<std::string> statuses) {
    std::string records = run;
    for (const std::string &status : statuses) {
      records += status;
    }
    return wrap(format_minor, records + event + event_end + run_end);
  };
  CHECK(!throws([&] {
    read_all(run_with({status_of({{"atmosphere", "applied"}}), status_of({{"qe", "pending"}})}),
             nullptr);
  }));
  CHECK(!throws([&] {
    read_all(run_with({status_of({{"qe", "pending"}}), status_of({{"qe", "applied"}})}), nullptr);
  }));
  CHECK(!throws([&] {
    read_all(run_with({status_of({{"qe", "pending"}}), status_of({{"qe", "envelope"}}),
                       status_of({{"qe", "applied"}})}),
             nullptr);
  }));
  CHECK(throws([&] {
    read_all(run_with({status_of({{"qe", "pending"}}), status_of({{"qe", "pending"}})}), nullptr);
  }));
  CHECK(throws([&] {
    read_all(run_with({status_of({{"qe", "applied"}}), status_of({{"qe", "pending"}})}), nullptr);
  }));
  CHECK(throws([&] {
    read_all(run_with({status_of({{"qe", "envelope"}}), status_of({{"qe", "envelope"}})}), nullptr);
  }));
  CHECK(throws([&] {
    read_all(run_with({status_of({{"qe", "envelope"}}), status_of({{"qe", "pending"}})}), nullptr);
  }));

  // A run status after the first event is rejected (the run scope is fixed by then).
  CHECK(throws([&] {
    read_all(wrap(format_minor, run + event + event_end + status_of({{"qe", "pending"}}) + run_end),
             nullptr);
  }));

  // Telescope scope: checked separately from the run scope, and only before its photons.
  const std::string telescope_block = run + status_of({{"qe", "applied"}}) + event + telescope +
                                      status_of({{"qe", "pending"}}) + photons + telescope_end +
                                      event_end + run_end;
  CHECK(!throws([&] { read_all(wrap(format_minor, telescope_block), nullptr); }));
  CHECK(throws([&] {
    read_all(wrap(format_minor, run + event + telescope + status_of({{"qe", "applied"}}) +
                                    status_of({{"qe", "pending"}}) + telescope_end + event_end +
                                    run_end),
             nullptr);
  }));

  // Malformed status records: missing or unknown state, a factor listed twice.
  Metadata no_state{{"factors", std::string("qe")}};
  CHECK(throws([&] {
    read_all(run_with({raw_metadata(RecordType::efficiency_status, no_state)}), nullptr);
  }));
  Metadata unknown{{"factors", std::string("qe")}, {"qe.state", std::string("done")}};
  CHECK(throws([&] {
    read_all(run_with({raw_metadata(RecordType::efficiency_status, unknown)}), nullptr);
  }));
  Metadata twice{{"factors", std::string("qe,qe")}, {"qe.state", std::string("pending")}};
  CHECK(throws(
      [&] { read_all(run_with({raw_metadata(RecordType::efficiency_status, twice)}), nullptr); }));
}

// Values the format forbids: non-finite columns, a direction that is not a unit vector, and an
// upward photon in a format 1.0 stream (cz > 0); upward photons are valid in format 1.1.
void value_checks() {
  const auto pair = [](std::uint32_t extra) {
    PhotonChunk chunk = make_chunk(1, Frame::ground, 0, 2, extra);
    return chunk;
  };
  const std::uint32_t cz = column_bit(Column::cz);
  const std::uint32_t zem = column_bit(Column::zem);
  // Valid: downward with cz, and without cz (implied downward).
  CHECK(!throws([&] { read_all(blocks_stream(Metadata{}, {pair(cz)}), nullptr); }));
  CHECK(!throws([&] { read_all(blocks_stream(Metadata{}, {pair(0)}), nullptr); }));

  for (const Column column : {Column::x, Column::y, Column::t, Column::w, Column::zem}) {
    PhotonChunk chunk = pair(cz | zem);
    chunk[column][1] = std::nanf("");
    CHECK(throws([&] { read_all(blocks_stream(Metadata{}, {chunk}), nullptr); }));
  }
  PhotonChunk inf = pair(0);
  inf[Column::t][0] = INFINITY;
  CHECK(throws([&] { read_all(blocks_stream(Metadata{}, {inf}), nullptr); }));

  // cz inconsistent with (cx, cy).
  PhotonChunk nonunit = pair(cz);
  nonunit[Column::cz][0] = -0.5F;
  CHECK(throws([&] { read_all(blocks_stream(Metadata{}, {nonunit}), nullptr); }));
  // Without cz, (cx, cy) must lie in the unit disc.
  PhotonChunk outside = pair(0);
  outside[Column::cx][0] = 0.9F;
  outside[Column::cy][0] = 0.9F;
  CHECK(throws([&] { read_all(blocks_stream(Metadata{}, {outside}), nullptr); }));

  // Upward photons: rejected in format 1.0, accepted from 1.1 on.
  PhotonChunk upward = pair(cz);
  for (std::size_t r = 0; r < upward.rows(); ++r) {
    upward[Column::cz][r] = -upward[Column::cz][r];
  }
  CHECK(throws([&] { read_all(wrap(0, block(photons_of(upward))), nullptr); }));
  CHECK(!throws([&] { read_all(wrap(1, block(photons_of(upward))), nullptr); }));
  CHECK(!throws([&] { read_all(wrap(2, block(photons_of(upward))), nullptr); }));
  // A horizontal photon (cz = 0) with cx = cy = 0 is not upward and is accepted in 1.0.
  PhotonChunk horizontal = pair(cz);
  for (std::size_t r = 0; r < horizontal.rows(); ++r) {
    horizontal[Column::cx][r] = 1.0F;
    horizontal[Column::cy][r] = 0.0F;
    horizontal[Column::cz][r] = 0.0F;
  }
  CHECK(!throws([&] { read_all(wrap(0, block(photons_of(horizontal))), nullptr); }));
}

// After an error the reader has failed: every later call throws, as the C API does.
void reader_failed_state() {
  std::istringstream in(stream_header(format_minor) + raw_metadata(RecordType::run, {}) +
                        raw_metadata(RecordType::run, {}) + stream_end());
  Reader reader(in);
  Record record;
  CHECK(reader.next(record) && record.type == RecordType::stream);
  CHECK(reader.next(record) && record.type == RecordType::run);
  CHECK(throws([&] { static_cast<void>(reader.next(record)); }));
  CHECK(throws([&] { static_cast<void>(reader.next(record)); }));

  // A stream of an unsupported major version: the first call throws and so does every later one.
  const Metadata header{{"format_major", std::int64_t{format_major + 1}},
                        {"format_minor", std::int64_t{0}},
                        {"producer", std::string("future")}};
  std::istringstream future(raw_metadata(RecordType::stream, header) + stream_end());
  Reader future_reader(future);
  CHECK(throws([&] { static_cast<void>(future_reader.next(record)); }));
  CHECK(throws([&] { static_cast<void>(future_reader.next(record)); }));
  CHECK(!future_reader.ended());
}

// The Writer refuses the records the Reader would refuse (the same nesting rules). A refused record
// writes none of its bytes, and a refused efficiency status changes no state of the writer.
void writer_refuses_nesting() {
  const Metadata telescope_1{{"telescope_id", std::int64_t{1}}, {"frame", std::int64_t{1}}};
  const Metadata id_1{{"telescope_id", std::int64_t{1}}};
  const PhotonChunk chunk = make_chunk(1, Frame::sky, 0, 2, 0);
  using Setup = std::function<void(Writer &)>;
  const Setup none = [](Writer &) {};
  const Setup in_run = [](Writer &w) { w.write_metadata(RecordType::run, {}); };
  const Setup in_event = [](Writer &w) {
    w.write_metadata(RecordType::run, {});
    w.write_metadata(RecordType::event, {});
  };
  const Setup in_telescope = [&](Writer &w) {
    in_event(w);
    w.write_metadata(RecordType::telescope, telescope_1);
  };
  const std::vector<std::tuple<const char *, Setup, std::function<void(Writer &)>>> refused = {
      {"event outside a run", none, [](Writer &w) { w.write_metadata(RecordType::event, {}); }},
      {"run inside a run", in_run, [](Writer &w) { w.write_metadata(RecordType::run, {}); }},
      {"run_end without run", none, [](Writer &w) { w.write_metadata(RecordType::run_end, {}); }},
      {"event inside an event", in_event,
       [](Writer &w) { w.write_metadata(RecordType::event, {}); }},
      {"telescope outside an event", in_run,
       [&](Writer &w) { w.write_metadata(RecordType::telescope, telescope_1); }},
      {"telescope_end without telescope", in_event,
       [&](Writer &w) { w.write_metadata(RecordType::telescope_end, id_1); }},
      {"photons outside a telescope", in_event, [&](Writer &w) { w.write_chunk(chunk); }},
      {"event_end inside a telescope", in_telescope,
       [](Writer &w) { w.write_metadata(RecordType::event_end, {}); }},
      {"run_end inside an event", in_event,
       [](Writer &w) { w.write_metadata(RecordType::run_end, {}); }},
      {"blob outside a run", none, [](Writer &w) { w.write_blob(1204, 0, {1, 2}); }},
      {"stream_end inside a run", in_run, [](Writer &w) { w.finish(); }},
  };
  for (const auto &[name, setup, bad] : refused) {
    std::ostringstream out;
    Writer writer(out, "test");
    setup(writer);
    const std::string before = out.str();
    CHECK(throws_invalid_argument([&] { bad(writer); }));
    CHECK(out.str() == before);
  }

  // A rejected run end leaves the run open: the stream closes once the run does.
  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.write_metadata(RecordType::run, {});
    CHECK(throws_invalid_argument([&] { writer.finish(); }));
    writer.write_metadata(RecordType::run_end, {});
    writer.finish();
  }
  CHECK(read_all(out.str(), nullptr) ==
        (std::vector<RecordType>{RecordType::stream, RecordType::run, RecordType::run_end}));

  // A refused efficiency status (a factor goes backwards) changes no state: the factor that was
  // valid in the same record is still pending afterwards.
  std::ostringstream status_out;
  {
    Writer writer(status_out, "test");
    writer.write_metadata(RecordType::run, {});
    writer.write_metadata(RecordType::efficiency_status, status_of_metadata({{"a", "envelope"}}));
    CHECK(throws_invalid_argument([&] {
      writer.write_metadata(RecordType::efficiency_status,
                            status_of_metadata({{"b", "pending"}, {"a", "pending"}}));
    }));
    writer.write_metadata(RecordType::efficiency_status, status_of_metadata({{"b", "pending"}}));
    writer.write_metadata(RecordType::run_end, {});
    writer.finish();
  }
  CHECK(!throws([&] { read_all(status_out.str(), nullptr); }));
}

// The reader's errors have their own types (the C API maps them by type): a truncated stream is
// a truncated (and malformed) stream, a stream of another major version an unsupported version.
void typed_errors() {
  const std::string bytes = blocks_stream(Metadata{}, {make_chunk(1, Frame::sky, 0, 10, 0)});
  CHECK(throws_as<TruncatedStreamError>(
      [&] { read_all(bytes.substr(0, bytes.size() - 5), nullptr); }));
  CHECK(throws_as<MalformedStreamError>(
      [&] { read_all(bytes.substr(0, bytes.size() - 5), nullptr); }));
  const Metadata header{{"format_major", std::int64_t{format_major + 1}},
                        {"format_minor", std::int64_t{0}},
                        {"producer", std::string("future")}};
  CHECK(throws_as<UnsupportedVersionError>(
      [&] { read_all(raw_metadata(RecordType::stream, header), nullptr); }));
  CHECK(throws_as<MalformedStreamError>([&] {
    read_all(
        wrap(format_minor, raw_metadata(RecordType::run, {}) + raw_metadata(RecordType::run, {})),
        nullptr);
  }));
}

} // namespace

int main() {
  return run_checks([] {
    metadata_round_trip();
    metadata_layout();
    metadata_malformed();
    stream_round_trip();
    unknown_record_skipped();
    empty_chunk();
    filtered_photons();
    filtered_photons_validated();
    invalid_chunk_rejected();
    bad_magic();
    truncated_record();
    every_prefix_is_incomplete();
    header_byte_flips_are_contained();
    wrong_format_major();
    format_minor_accepted();
    header_required_once();
    clean_end_required();
    path_round_trip();
    format_minor_columns();
    ground_representation_helper();
    sphere_rows_validation();
    sphere_stream_round_trip();
    nesting_rejected();
    efficiency_transitions();
    value_checks();
    reader_failed_state();
    writer_refuses_nesting();
    typed_errors();
  });
}
