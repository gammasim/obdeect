// Compares the C API with the C++ reader. Usage: test_capi_cpp <stream written by test_capi>
// <conformance dir>. Both readers must see the same records; the conformance streams must
// decode as described in cpp/tests/data/conformance/expected.json (values repeated here).

#include "check.hpp"

#include "obdeect/stream/obdeect_stream.h"
#include "obdeect/stream/stream.hpp"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

using namespace obdeect::stream;

namespace {

void compare_with_c_reader(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  CHECK(static_cast<bool>(in));
  Reader reader(in);
  Record record;
  obdeect_stream_reader *c = nullptr;
  CHECK(obdeect_stream_reader_open(path.c_str(), &c) == OBDEECT_STREAM_OK);
  int n_records = 0;
  while (reader.next(record)) {
    int type = 0;
    CHECK(obdeect_stream_reader_next(c, &type) == OBDEECT_STREAM_OK);
    CHECK(type == static_cast<int>(record.type));
    ++n_records;
    if (record.type == RecordType::photons) {
      std::uint32_t telescope = 0, frame = 0, rows = 0, mask = 0;
      std::uint64_t first = 0;
      CHECK(obdeect_stream_reader_photons(c, &telescope, &frame, &first, &rows, &mask) == OBDEECT_STREAM_OK);
      CHECK(telescope == record.chunk.telescope_id);
      CHECK(frame == static_cast<std::uint32_t>(record.chunk.frame));
      CHECK(first == record.chunk.first_index);
      CHECK(rows == record.chunk.rows());
      CHECK(mask == record.chunk.mask());
      for (std::uint32_t i = 0; i < n_columns; ++i) {
        const float *column = obdeect_stream_reader_column(c, static_cast<int>(i));
        const auto &expected = record.chunk.columns[i];
        CHECK((column != nullptr) == !expected.empty());
        if (column != nullptr) {
          CHECK(std::memcmp(column, expected.data(), expected.size() * sizeof(float)) == 0);
        }
      }
    } else if (record.type == RecordType::blob) {
      std::uint32_t t = 0, id = 0;
      const std::uint8_t *data = nullptr;
      std::size_t size = 0;
      CHECK(obdeect_stream_reader_blob(c, &t, &id, &data, &size) == OBDEECT_STREAM_OK);
      CHECK(t == record.blob_eventio_type && id == record.blob_eventio_id);
      CHECK(size == record.blob.size() && std::memcmp(data, record.blob.data(), size) == 0);
    } else {
      std::size_t n = 0;
      CHECK(obdeect_stream_reader_n_keys(c, &n) == OBDEECT_STREAM_OK);
      CHECK(n == record.metadata.size());
      std::size_t i = 0;
      for (const auto &[key, value] : record.metadata) {
        CHECK(key == obdeect_stream_reader_key(c, i++));
        int value_type = 0;
        CHECK(obdeect_stream_reader_meta_type(c, key.c_str(), &value_type) == OBDEECT_STREAM_OK);
        CHECK(value_type == static_cast<int>(value.index()) + 1);
        if (const auto *integer = std::get_if<std::int64_t>(&value)) {
          std::int64_t got = 0;
          CHECK(obdeect_stream_reader_meta_int(c, key.c_str(), &got) == OBDEECT_STREAM_OK && got == *integer);
        } else if (const auto *real = std::get_if<double>(&value)) {
          double got = 0.0;
          CHECK(obdeect_stream_reader_meta_float(c, key.c_str(), &got) == OBDEECT_STREAM_OK && got == *real);
        } else if (const auto *text = std::get_if<std::string>(&value)) {
          const char *got = nullptr;
          std::size_t length = 0;
          CHECK(obdeect_stream_reader_meta_string(c, key.c_str(), &got, &length) == OBDEECT_STREAM_OK);
          CHECK(*text == std::string(got, length));
        } else {
          const auto &array = std::get<std::vector<double>>(value);
          const double *got = nullptr;
          std::size_t length = 0;
          CHECK(obdeect_stream_reader_meta_array(c, key.c_str(), &got, &length) == OBDEECT_STREAM_OK);
          CHECK(length == array.size() &&
                std::memcmp(got, array.data(), length * sizeof(double)) == 0);
        }
      }
    }
  }
  int type = 0;
  CHECK(obdeect_stream_reader_next(c, &type) == OBDEECT_STREAM_END);
  CHECK(n_records == 9);
  obdeect_stream_reader_close(c);
}

std::vector<Record> read_all(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  Reader reader(in);
  std::vector<Record> records;
  Record record;
  while (reader.next(record)) {
    records.push_back(record);
  }
  if (!reader.ended()) {
    throw std::runtime_error("gstrickd: truncated stream (no stream_end record)");
  }
  return records;
}

void check_conformance(const std::string &directory) {
  const auto records = read_all(directory + "/valid.ps");
  const RecordType types[] = {
      RecordType::stream,    RecordType::run,       RecordType::efficiency_status,
      RecordType::event,     RecordType::telescope, RecordType::photons,
      RecordType::photons,   RecordType::blob,      RecordType::telescope_end,
      RecordType::event_end, RecordType::run_end};
  CHECK(records.size() == 11);
  for (std::size_t i = 0; i < records.size() && i < 11; ++i) {
    CHECK(records[i].type == types[i]);
  }
  if (records.size() != 11) {
    return;
  }
  CHECK(std::get<std::string>(records[0].metadata.at("producer")) == "conformance");
  CHECK(std::get<std::int64_t>(records[1].metadata.at("run_number")) == 7);
  CHECK(std::get<double>(records[1].metadata.at("energy_min_gev")) == 0.25);
  CHECK(std::get<std::string>(records[1].metadata.at("name")) == "site");
  CHECK((std::get<std::vector<double>>(records[1].metadata.at("telescope_x_m")) ==
         std::vector<double>{0.0, 10.5}));
  CHECK(std::get<std::string>(records[2].metadata.at("qe.state")) == "applied");
  const PhotonChunk &chunk = records[5].chunk;
  CHECK(chunk.telescope_id == 2 && chunk.frame == Frame::ground && chunk.first_index == 0);
  CHECK(chunk.rows() == 3);
  CHECK(chunk.mask() == (required_columns | column_bit(Column::zem)));
  CHECK((chunk[Column::x] == std::vector<float>{0.0F, 1.0F, 2.0F}));
  CHECK((chunk[Column::w] == std::vector<float>{1.0F, 1.0F, 0.5F}));
  CHECK((chunk[Column::zem] == std::vector<float>{8000.0F, 8100.0F, 8200.0F}));
  CHECK(records[6].chunk.rows() == 0 && records[6].chunk.first_index == 3);
  CHECK(records[7].blob_eventio_type == 1204 && records[7].blob_eventio_id == 9);
  CHECK((records[7].blob == std::vector<std::uint8_t>{1, 2, 3, 255}));
  CHECK(std::get<double>(records[8].metadata.at("photons")) == 2.5);

  const struct {
    const char *file;
    const char *message;
  } invalid[] = {{"bad_magic.ps", "bad magic"},
                 {"unknown_major.ps", "unsupported format major version"},
                 {"truncated_chunk.ps", "truncated"},
                 {"too_many_rows.ps", "2^16 rows"}};
  for (const auto &entry : invalid) {
    try {
      read_all(directory + "/" + entry.file);
      CHECK(false);
    } catch (const std::runtime_error &error) {
      CHECK(std::string(error.what()).find(entry.message) != std::string::npos);
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: test_capi_cpp <stream> <conformance dir>\n");
    return 2;
  }
  return run_checks([&] {
    compare_with_c_reader(argv[1]);
    check_conformance(argv[2]);
  });
}
