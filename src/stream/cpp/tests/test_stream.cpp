#include "check.hpp"

#include "gstrickd/stream.hpp"

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace gstrickd;

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
  Metadata ledger{{"factors", std::string("qe,mirror")}, {"qe.state", std::string("applied")}};
  Metadata end{{"photons", 12.5}};

  const PhotonChunk first = make_chunk(3, Frame::ground, 0, 100, 0);
  const PhotonChunk second = make_chunk(3, Frame::telescope, 100, 50,
                                        column_bit(Column::zem) | column_bit(Column::wavelength));

  std::ostringstream out;
  {
    Writer writer(out, "test");
    writer.write_metadata(RecordType::run, run);
    writer.write_metadata(RecordType::ledger, ledger);
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
      RecordType::stream,        RecordType::run,       RecordType::ledger,  RecordType::blob,
      RecordType::event,         RecordType::telescope, RecordType::photons, RecordType::photons,
      RecordType::telescope_end, RecordType::event_end, RecordType::run_end};
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
  CHECK(reader.next(record) && record.metadata == ledger);
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

} // namespace

int main() {
  return run_checks([] {
    metadata_round_trip();
    metadata_layout();
    metadata_malformed();
    stream_round_trip();
    unknown_record_skipped();
    empty_chunk();
    invalid_chunk_rejected();
    bad_magic();
    truncated_record();
    wrong_format_major();
    header_required_once();
    clean_end_required();
  });
}
