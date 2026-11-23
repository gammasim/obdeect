#include "check.hpp"

#include "obdeect/stream/stream.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
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

// One complete record with the given type and payload, little-endian header.
std::string raw_record(std::uint16_t type, const std::string &payload) {
  std::string out;
  const std::uint32_t magic = record_magic;
  for (int i = 0; i < 4; ++i) {
    out.push_back(static_cast<char>(magic >> (8 * i)));
  }
  out.push_back(static_cast<char>(type));
  out.push_back(static_cast<char>(type >> 8));
  out.append(2, '\0');
  const auto length = static_cast<std::uint64_t>(payload.size());
  for (int i = 0; i < 8; ++i) {
    out.push_back(static_cast<char>(length >> (8 * i)));
  }
  return out + payload;
}

std::string as_string(const std::vector<std::uint8_t> &bytes) {
  return {bytes.begin(), bytes.end()};
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
  return chunk;
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
                             {"qe.state", std::string("applied")}};
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
  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.write_chunk(empty);
    writer.finish();
  }
  std::vector<PhotonChunk> chunks;
  CHECK(read_all(out.str(), &chunks) ==
        (std::vector<RecordType>{RecordType::stream, RecordType::photons}));
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

// Little-endian bytes of the low n bytes of value.
std::string le_bytes(std::uint64_t value, std::size_t n) {
  std::string out;
  for (std::size_t i = 0; i < n; ++i) {
    out.push_back(static_cast<char>(value >> (8 * i)));
  }
  return out;
}

// Photons payload with its fixed header and no column data.
std::string photons_header(std::uint32_t telescope, std::uint32_t frame, std::uint32_t mask,
                           std::uint32_t rows) {
  return le_bytes(telescope, 4) + le_bytes(frame, 4) + le_bytes(mask, 4) + le_bytes(rows, 4) +
         le_bytes(0, 8);
}

void filtered_photons() {
  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.write_chunk(make_chunk(1, Frame::sky, 7, 10, 0));
    writer.write_chunk(make_chunk(2, Frame::telescope, 100, 5, 0));
    writer.write_chunk(make_chunk(3, Frame::camera, 555, 4, 0));
    writer.finish();
  }
  std::istringstream in(out.str());
  Reader reader(in);
  reader.keep_photons_of({2});
  Record record;
  CHECK(reader.next(record) && record.type == RecordType::stream);
  // Skipped chunks report their own header fields, also after a kept chunk reused the record.
  CHECK(reader.next(record) && record.type == RecordType::photons);
  CHECK(record.chunk.telescope_id == 1 && record.chunk.frame == Frame::sky);
  CHECK(record.chunk.first_index == 7 && record.chunk.rows() == 0 && record.chunk.mask() == 0);
  CHECK(reader.next(record) && record.type == RecordType::photons);
  CHECK(record.chunk.telescope_id == 2 && record.chunk.frame == Frame::telescope);
  CHECK(record.chunk.first_index == 100 && record.chunk.rows() == 5);
  CHECK(reader.next(record) && record.type == RecordType::photons);
  CHECK(record.chunk.telescope_id == 3 && record.chunk.frame == Frame::camera);
  CHECK(record.chunk.first_index == 555 && record.chunk.rows() == 0 && record.chunk.mask() == 0);
  CHECK(!reader.next(record) && reader.ended());
}

void filtered_photons_validated() {
  // Unselected chunks are still validated: unknown frame or wrong payload size is rejected.
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
    std::ostringstream bad;
    Writer writer(bad, "test");
    bad << raw_record(static_cast<std::uint16_t>(photons), payload);
    writer.finish();
    CHECK(throws([&] { read_keeping_two(bad.str()); }));
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
  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.write_chunk(make_chunk(1, Frame::sky, 0, 10, 0));
    writer.finish();
  }
  const std::string bytes = out.str();
  CHECK(throws([&] { read_all(bytes.substr(0, bytes.size() - 5), nullptr); }));
  CHECK(throws([&] { read_all(bytes.substr(0, bytes.size() - 20), nullptr); }));
  CHECK(throws([&] { read_all(bytes.substr(0, 10), nullptr); }));
}

void wrong_format_major() {
  Metadata header{{"format_major", std::int64_t{format_major + 1}},
                  {"format_minor", std::int64_t{0}},
                  {"producer", std::string("future")}};
  const std::string bytes = raw_record(1, as_string(encode_metadata(header)));
  CHECK(throws([&] { read_all(bytes, nullptr); }));
}

// Readers accept every format_minor of their major version: 1.0 streams (no upward photons)
// and the current 1.1.
void format_minor_accepted() {
  for (const std::int64_t minor : {std::int64_t{0}, std::int64_t{format_minor}}) {
    const Metadata header{{"format_major", std::int64_t{format_major}},
                          {"format_minor", minor},
                          {"producer", std::string("test")}};
    const std::string bytes =
        raw_record(1, as_string(encode_metadata(header))) + raw_record(15, "");
    CHECK(!throws([&] { read_all(bytes, nullptr); }));
  }
}

void header_required_once() {
  // Records before the stream header, and a second header, are malformed.
  CHECK(throws([&] { read_all(raw_record(15, ""), nullptr); }));
  CHECK(throws([&] { read_all(raw_record(2, as_string(encode_metadata({}))), nullptr); }));
  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.finish();
  }
  std::string twice = out.str();
  twice = twice.substr(0, twice.size() - 16) + twice; // header, header, stream_end
  CHECK(throws([&] { read_all(twice, nullptr); }));
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

// Little-endian float32 bytes of the given values.
std::string float_bytes(float value) { return le_bytes(std::bit_cast<std::uint32_t>(value), 4); }

// Stream header record of the given format_minor (major 1).
std::string stream_header(std::int64_t minor) {
  const Metadata header{{"format_major", std::int64_t{format_major}},
                        {"format_minor", minor},
                        {"producer", std::string("test")}};
  return raw_record(1, as_string(encode_metadata(header)));
}

// Photons record of telescope 1 in the ground frame: every column of the mask holds value + index.
std::string raw_photons(std::uint32_t mask, std::uint32_t rows, float value) {
  std::string payload = photons_header(1, static_cast<std::uint32_t>(Frame::ground), mask, rows);
  for (std::uint32_t i = 0; i < n_columns; ++i) {
    if (((mask >> i) & 1U) == 0) {
      continue;
    }
    for (std::uint32_t r = 0; r < rows; ++r) {
      payload += float_bytes(value + static_cast<float>(i));
    }
  }
  return raw_record(static_cast<std::uint16_t>(RecordType::photons), payload);
}

// Run record with the given metadata (written raw, so that the writer's checks are not used).
std::string raw_run(const Metadata &run) {
  return raw_record(static_cast<std::uint16_t>(RecordType::run), as_string(encode_metadata(run)));
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
    std::ostringstream out;
    {
      Writer writer(out, "test");
      writer.write_chunk(chunk);
      writer.finish();
    }
    std::vector<PhotonChunk> chunks;
    read_all(out.str(), &chunks);
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
    // A photons record of two rows whose path values are bad (the other columns are 1.0).
    const std::uint32_t mask = required_columns | column_bit(Column::path);
    std::string payload = photons_header(1, static_cast<std::uint32_t>(Frame::ground), mask, 2);
    for (std::uint32_t i = 0; i < n_columns; ++i) {
      for (std::uint32_t r = 0; r < 2 && ((mask >> i) & 1U); ++r) {
        payload += float_bytes(i == static_cast<std::uint32_t>(Column::path) ? bad : 1.0F);
      }
    }
    const std::string bytes = stream_header(format_minor) +
                              raw_record(static_cast<std::uint16_t>(RecordType::photons), payload) +
                              raw_record(15, "");
    CHECK(throws([&] { read_all(bytes, nullptr); }));
  }
}

// Streams of format 1.1 (no path column, no sphere runs) are read by the 1.2 reader; a path column
// or a sphere run in them is malformed. A later minor version is read (minor versions are
// compatible within a major version); a negative one is rejected.
void format_minor_columns() {
  const std::string end = raw_record(15, "");
  const std::uint32_t plain =
      required_columns | column_bit(Column::zem) | column_bit(Column::wavelength);
  std::vector<PhotonChunk> chunks;
  CHECK(!throws([&] { read_all(stream_header(1) + raw_photons(plain, 3, 1.0F) + end, &chunks); }));
  CHECK(chunks.size() == 1);
  if (chunks.size() == 1) {
    CHECK(chunks[0].rows() == 3 && !chunks[0].has(Column::path));
    CHECK(chunks[0].has(Column::zem) && chunks[0].has(Column::wavelength));
    CHECK_NEAR(chunks[0][Column::x][2], 1.0, 0.0); // raw_photons: value 1.0 + column index 0
  }
  CHECK(throws([&] {
    read_all(stream_header(1) + raw_photons(plain | column_bit(Column::path), 3, 1.0F) + end,
             nullptr);
  }));
  const Metadata sphere{{"ground_representation", std::string("sphere")}};
  CHECK(throws([&] { read_all(stream_header(1) + raw_run(sphere) + end, nullptr); }));
  CHECK(!throws([&] { read_all(stream_header(2) + raw_run(sphere) + end, nullptr); }));
  CHECK(!throws([&] { read_all(stream_header(3) + raw_run(sphere) + end, nullptr); }));
  CHECK(throws([&] { read_all(stream_header(-1) + end, nullptr); }));
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
    read_all(stream_header(2) + raw_run(Metadata{{"ground_representation", std::string("flat")}}) +
                 raw_record(15, ""),
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
  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.write_metadata(RecordType::run,
                          Metadata{{"ground_representation", std::string("sphere")}});
    writer.write_chunk(sphere_rows(2, column_bit(Column::path)));
    writer.finish();
  }
  std::vector<PhotonChunk> chunks;
  CHECK(read_all(out.str(), &chunks) ==
        (std::vector<RecordType>{RecordType::stream, RecordType::run, RecordType::photons}));
  CHECK(chunks.size() == 1);
  if (chunks.size() == 1) {
    CHECK(chunks[0].columns == sphere_rows(2, column_bit(Column::path)).columns);
    CHECK(chunks[0].frame == Frame::ground);
  }
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
    wrong_format_major();
    format_minor_accepted();
    header_required_once();
    clean_end_required();
    path_round_trip();
    format_minor_columns();
    ground_representation_helper();
    sphere_rows_validation();
    sphere_stream_round_trip();
  });
}
