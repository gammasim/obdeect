// C API of gstrickd (cpp/include/obdeect/stream/obdeect_stream.h): a thin layer over the header-only C++
// reader and writer. Every exported function catches all exceptions.

#include "obdeect/stream/obdeect_stream.h"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stdin_buffer.hpp"
#include "obdeect/stream/stream.hpp"

#include <cerrno>
#include <cstring>
#include <exception>
#include <fstream>
#include <istream>
#include <limits>
#include <memory>
#include <new>
#include <ostream>
#include <streambuf>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <unistd.h>

#ifndef OBDEECT_STREAM_PACKAGE_VERSION
#define OBDEECT_STREAM_PACKAGE_VERSION "0.0.0"
#endif

static_assert(OBDEECT_STREAM_MAX_CHUNK_ROWS == obdeect::stream::max_chunk_rows);
static_assert(OBDEECT_STREAM_N_COLUMNS == obdeect::stream::n_columns);

namespace {

thread_local std::string last_error_message;

// Never throws: if the message cannot be stored it is empty.
int fail(int code, const char *message) noexcept {
  try {
    last_error_message = message;
  } catch (...) {
    last_error_message.clear();
  }
  return code;
}

int fail_key(int code, const char *prefix, const char *key, const char *suffix) noexcept {
  try {
    last_error_message = std::string(prefix) + key + suffix;
  } catch (...) {
    last_error_message.clear();
  }
  return code;
}

// Maps an escaping exception to a return code and records its message.
int from_exception() noexcept {
  try {
    throw;
  } catch (const std::bad_alloc &) {
    return fail(OBDEECT_STREAM_ERROR_MEMORY, "out of memory");
  } catch (const std::invalid_argument &error) {
    return fail(OBDEECT_STREAM_ERROR_ARGUMENT, error.what());
  } catch (const std::runtime_error &error) {
    const char *text = error.what();
    if (std::strstr(text, "unsupported format major version") != nullptr) {
      return fail(OBDEECT_STREAM_ERROR_VERSION, text);
    }
    if (std::strstr(text, "read error") != nullptr) {
      return fail(OBDEECT_STREAM_ERROR_IO, text);
    }
    return fail(OBDEECT_STREAM_ERROR_FORMAT, text);
  } catch (const std::exception &error) {
    return fail(OBDEECT_STREAM_ERROR_INTERNAL, error.what());
  } catch (...) {
    return fail(OBDEECT_STREAM_ERROR_INTERNAL, "unknown exception");
  }
}

int bad_argument(const char *what) { return fail(OBDEECT_STREAM_ERROR_ARGUMENT, what); }

// Unbuffered std::streambuf writing to a file descriptor.
class FdOutBuffer : public std::streambuf {
public:
  explicit FdOutBuffer(int fd) : fd_(fd) {}

protected:
  std::streamsize xsputn(const char *data, std::streamsize n) override {
    std::streamsize done = 0;
    while (done < n) {
      const auto got =
          static_cast<long long>(::write(fd_, data + done, static_cast<std::size_t>(n - done)));
      if (got < 0) {
        if (errno == EINTR) {
          continue;
        }
        break;
      }
      done += got;
    }
    return done;
  }
  int_type overflow(int_type c) override {
    if (traits_type::eq_int_type(c, traits_type::eof())) {
      return traits_type::not_eof(c);
    }
    const char ch = traits_type::to_char_type(c);
    return xsputn(&ch, 1) == 1 ? c : traits_type::eof();
  }

private:
  int fd_;
};

} // namespace

struct obdeect_stream_reader {
  std::unique_ptr<std::ifstream> file;
  obdeect::stream::StdinBuffer fd_buffer;
  std::istream fd_stream;
  std::istream *input{nullptr};
  std::unique_ptr<obdeect::stream::Reader> reader;
  obdeect::stream::Record record;
  std::vector<const std::string *> keys;
  bool has_record{false};
  bool failed{false};

  explicit obdeect_stream_reader(int fd)
      : fd_buffer(fd < 0 ? 1 : std::size_t{1} << 20, fd), fd_stream(&fd_buffer) {}

  bool has_metadata() const {
    return has_record && (record.type == obdeect::stream::RecordType::stream ||
                          obdeect::stream::is_metadata_record(record.type));
  }

  const obdeect::stream::Value *find(const char *key) const {
    const auto it = record.metadata.find(key);
    return it == record.metadata.end() ? nullptr : &it->second;
  }
};

struct obdeect_stream_writer {
  std::unique_ptr<std::ofstream> file;
  FdOutBuffer fd_buffer;
  std::ostream fd_stream;
  std::ostream *output{nullptr};
  std::unique_ptr<obdeect::stream::Writer> writer;
  obdeect::stream::Metadata pending;
  obdeect::stream::PhotonChunk chunk;

  explicit obdeect_stream_writer(int fd) : fd_buffer(fd), fd_stream(&fd_buffer) {}
};

namespace {

int check_metadata_access(const obdeect_stream_reader *reader, const char *key) {
  if (reader == nullptr || key == nullptr) {
    return bad_argument("NULL argument");
  }
  if (!reader->has_metadata()) {
    return fail(OBDEECT_STREAM_ERROR_STATE, "the current record has no metadata");
  }
  return OBDEECT_STREAM_OK;
}

template <class T> const T *typed_value(const obdeect_stream_reader *reader, const char *key, int &code) {
  code = check_metadata_access(reader, key);
  if (code != OBDEECT_STREAM_OK) {
    return nullptr;
  }
  const obdeect::stream::Value *value = reader->find(key);
  if (value == nullptr) {
    code = fail_key(OBDEECT_STREAM_ERROR_NOT_FOUND, "no metadata key \'", key, "\'");
    return nullptr;
  }
  const T *typed = std::get_if<T>(value);
  if (typed == nullptr) {
    code = fail_key(OBDEECT_STREAM_ERROR_TYPE, "metadata key '", key, "' has another type");
  }
  return typed;
}

int open_reader(std::unique_ptr<obdeect_stream_reader> handle, obdeect_stream_reader **out) {
  handle->reader = std::make_unique<obdeect::stream::Reader>(*handle->input);
  *out = handle.release();
  return OBDEECT_STREAM_OK;
}

int open_writer(std::unique_ptr<obdeect_stream_writer> handle, const char *producer,
                obdeect_stream_writer **out) {
  handle->writer = std::make_unique<obdeect::stream::Writer>(*handle->output, producer);
  if (handle->output->fail()) {
    return fail(OBDEECT_STREAM_ERROR_IO, "cannot write the stream record");
  }
  *out = handle.release();
  return OBDEECT_STREAM_OK;
}

// Writer calls that write to the stream: report a failed output stream as an I/O error.
int check_output(const obdeect_stream_writer *writer) {
  return writer->output->fail() ? fail(OBDEECT_STREAM_ERROR_IO, "write failed (e.g. disk full)")
                                : OBDEECT_STREAM_OK;
}

template <class T> int set_pending(obdeect_stream_writer *writer, const char *key, T &&value) {
  if (writer == nullptr || key == nullptr) {
    return bad_argument("NULL argument");
  }
  writer->pending[key] = std::forward<T>(value);
  return OBDEECT_STREAM_OK;
}

} // namespace

extern "C" {

const char *obdeect_stream_version(void) { return OBDEECT_STREAM_PACKAGE_VERSION; }
int obdeect_stream_abi_version(void) { return OBDEECT_STREAM_C_ABI_VERSION; }
int obdeect_stream_format_major(void) { return obdeect::stream::format_major; }
int obdeect_stream_format_minor(void) { return obdeect::stream::format_minor; }
const char *obdeect_stream_last_error(void) { return last_error_message.c_str(); }

const char *obdeect_stream_record_type_name(int record_type) {
  using obdeect::stream::RecordType;
  switch (static_cast<RecordType>(record_type)) {
  case RecordType::stream:
  case RecordType::run:
  case RecordType::event:
  case RecordType::telescope:
  case RecordType::photons:
  case RecordType::telescope_end:
  case RecordType::event_end:
  case RecordType::run_end:
  case RecordType::blob:
  case RecordType::efficiency_status:
  case RecordType::stream_end:
    return obdeect::stream::record_name(static_cast<RecordType>(record_type));
  }
  return nullptr;
}

const char *obdeect_stream_column_name(int column) {
  if (column < 0 || column >= static_cast<int>(obdeect::stream::n_columns)) {
    return nullptr;
  }
  return obdeect::stream::column_name(static_cast<obdeect::stream::Column>(column));
}

int obdeect_stream_reader_open(const char *path, obdeect_stream_reader **out) {
  try {
    if (path == nullptr || out == nullptr) {
      return bad_argument("NULL argument");
    }
    *out = nullptr;
    auto handle = std::make_unique<obdeect_stream_reader>(-1);
    handle->file = std::make_unique<std::ifstream>(path, std::ios::binary);
    if (!*handle->file) {
      return fail_key(OBDEECT_STREAM_ERROR_IO, "cannot open ", path, "");
    }
    handle->input = handle->file.get();
    return open_reader(std::move(handle), out);
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_reader_open_fd(int fd, obdeect_stream_reader **out) {
  try {
    if (fd < 0 || out == nullptr) {
      return bad_argument("invalid file descriptor or NULL argument");
    }
    *out = nullptr;
    auto handle = std::make_unique<obdeect_stream_reader>(fd);
    handle->input = &handle->fd_stream;
    return open_reader(std::move(handle), out);
  } catch (...) {
    return from_exception();
  }
}

void obdeect_stream_reader_close(obdeect_stream_reader *reader) {
  try {
    delete reader;
  } catch (...) {
  }
}

int obdeect_stream_reader_next(obdeect_stream_reader *reader, int *record_type) {
  try {
    if (reader == nullptr || record_type == nullptr) {
      return bad_argument("NULL argument");
    }
    *record_type = 0;
    reader->has_record = false;
    if (reader->failed) {
      return fail(OBDEECT_STREAM_ERROR_STATE, "the reader failed earlier");
    }
    if (!reader->reader->next(reader->record)) {
      if (!reader->reader->ended()) {
        reader->failed = true;
        return fail(OBDEECT_STREAM_ERROR_FORMAT, "truncated stream (no stream_end record)");
      }
      if (!reader->reader->at_clean_end()) {
        reader->failed = true;
        return fail(OBDEECT_STREAM_ERROR_FORMAT, "data after stream_end");
      }
      return OBDEECT_STREAM_END;
    }
    reader->has_record = true;
    reader->keys.clear();
    if (reader->has_metadata()) {
      for (const auto &entry : reader->record.metadata) {
        reader->keys.push_back(&entry.first);
      }
    }
    *record_type = static_cast<int>(reader->record.type);
    return OBDEECT_STREAM_OK;
  } catch (...) {
    reader->failed = true;
    return from_exception();
  }
}

int obdeect_stream_reader_n_keys(const obdeect_stream_reader *reader, size_t *n) {
  if (reader == nullptr || n == nullptr) {
    return bad_argument("NULL argument");
  }
  if (!reader->has_metadata()) {
    return fail(OBDEECT_STREAM_ERROR_STATE, "the current record has no metadata");
  }
  *n = reader->keys.size();
  return OBDEECT_STREAM_OK;
}

const char *obdeect_stream_reader_key(const obdeect_stream_reader *reader, size_t i) {
  if (reader == nullptr || !reader->has_metadata() || i >= reader->keys.size()) {
    return nullptr;
  }
  return reader->keys[i]->c_str();
}

int obdeect_stream_reader_meta_type(const obdeect_stream_reader *reader, const char *key, int *value_type) {
  try {
    const int code = check_metadata_access(reader, key);
    if (code != OBDEECT_STREAM_OK) {
      return code;
    }
    if (value_type == nullptr) {
      return bad_argument("NULL argument");
    }
    const obdeect::stream::Value *value = reader->find(key);
    if (value == nullptr) {
      return fail_key(OBDEECT_STREAM_ERROR_NOT_FOUND, "no metadata key \'", key, "\'");
    }
    *value_type = static_cast<int>(value->index()) + 1; // variant order = OBDEECT_STREAM_VALUE_*
    return OBDEECT_STREAM_OK;
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_reader_meta_int(const obdeect_stream_reader *reader, const char *key, int64_t *value) {
  try {
    int code = OBDEECT_STREAM_OK;
    const auto *typed = typed_value<std::int64_t>(reader, key, code);
    if (typed == nullptr) {
      return code;
    }
    if (value == nullptr) {
      return bad_argument("NULL argument");
    }
    *value = *typed;
    return OBDEECT_STREAM_OK;
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_reader_meta_float(const obdeect_stream_reader *reader, const char *key, double *value) {
  try {
    const int code = check_metadata_access(reader, key);
    if (code != OBDEECT_STREAM_OK) {
      return code;
    }
    if (value == nullptr) {
      return bad_argument("NULL argument");
    }
    const obdeect::stream::Value *found = reader->find(key);
    if (found == nullptr) {
      return fail_key(OBDEECT_STREAM_ERROR_NOT_FOUND, "no metadata key \'", key, "\'");
    }
    if (const auto *real = std::get_if<double>(found)) {
      *value = *real;
    } else if (const auto *integer = std::get_if<std::int64_t>(found)) {
      *value = static_cast<double>(*integer);
    } else {
      return fail_key(OBDEECT_STREAM_ERROR_TYPE, "metadata key '", key, "' is not a number");
    }
    return OBDEECT_STREAM_OK;
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_reader_meta_string(const obdeect_stream_reader *reader, const char *key, const char **value,
                                size_t *length) {
  try {
    int code = OBDEECT_STREAM_OK;
    const auto *typed = typed_value<std::string>(reader, key, code);
    if (typed == nullptr) {
      return code;
    }
    if (value == nullptr) {
      return bad_argument("NULL argument");
    }
    *value = typed->c_str();
    if (length != nullptr) {
      *length = typed->size();
    }
    return OBDEECT_STREAM_OK;
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_reader_meta_array(const obdeect_stream_reader *reader, const char *key,
                               const double **values, size_t *n) {
  try {
    int code = OBDEECT_STREAM_OK;
    const auto *typed = typed_value<std::vector<double>>(reader, key, code);
    if (typed == nullptr) {
      return code;
    }
    if (values == nullptr || n == nullptr) {
      return bad_argument("NULL argument");
    }
    *values = typed->data();
    *n = typed->size();
    return OBDEECT_STREAM_OK;
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_reader_photons(const obdeect_stream_reader *reader, uint32_t *telescope_id, uint32_t *frame,
                            uint64_t *first_index, uint32_t *rows, uint32_t *column_mask) {
  if (reader == nullptr) {
    return bad_argument("NULL argument");
  }
  if (!reader->has_record || reader->record.type != obdeect::stream::RecordType::photons) {
    return fail(OBDEECT_STREAM_ERROR_STATE, "the current record is not photons");
  }
  const obdeect::stream::PhotonChunk &chunk = reader->record.chunk;
  if (telescope_id != nullptr) {
    *telescope_id = chunk.telescope_id;
  }
  if (frame != nullptr) {
    *frame = static_cast<std::uint32_t>(chunk.frame);
  }
  if (first_index != nullptr) {
    *first_index = chunk.first_index;
  }
  if (rows != nullptr) {
    *rows = static_cast<std::uint32_t>(chunk.rows());
  }
  if (column_mask != nullptr) {
    *column_mask = chunk.mask();
  }
  return OBDEECT_STREAM_OK;
}

const float *obdeect_stream_reader_column(const obdeect_stream_reader *reader, int column) {
  if (reader == nullptr || !reader->has_record ||
      reader->record.type != obdeect::stream::RecordType::photons || column < 0 ||
      column >= static_cast<int>(obdeect::stream::n_columns)) {
    return nullptr;
  }
  const auto &values = reader->record.chunk.columns[static_cast<std::size_t>(column)];
  return values.empty() ? nullptr : values.data();
}

int obdeect_stream_reader_blob(const obdeect_stream_reader *reader, uint32_t *eventio_type,
                         uint32_t *eventio_id, const uint8_t **data, size_t *size) {
  if (reader == nullptr) {
    return bad_argument("NULL argument");
  }
  if (!reader->has_record || reader->record.type != obdeect::stream::RecordType::blob) {
    return fail(OBDEECT_STREAM_ERROR_STATE, "the current record is not a blob");
  }
  if (eventio_type != nullptr) {
    *eventio_type = reader->record.blob_eventio_type;
  }
  if (eventio_id != nullptr) {
    *eventio_id = reader->record.blob_eventio_id;
  }
  if (data != nullptr) {
    *data = reader->record.blob.data();
  }
  if (size != nullptr) {
    *size = reader->record.blob.size();
  }
  return OBDEECT_STREAM_OK;
}

int obdeect_stream_writer_open(const char *path, const char *producer, obdeect_stream_writer **out) {
  try {
    if (path == nullptr || producer == nullptr || out == nullptr) {
      return bad_argument("NULL argument");
    }
    *out = nullptr;
    auto handle = std::make_unique<obdeect_stream_writer>(-1);
    handle->file = std::make_unique<std::ofstream>(path, std::ios::binary);
    if (!*handle->file) {
      return fail_key(OBDEECT_STREAM_ERROR_IO, "cannot open ", path, " for writing");
    }
    handle->output = handle->file.get();
    return open_writer(std::move(handle), producer, out);
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_writer_open_fd(int fd, const char *producer, obdeect_stream_writer **out) {
  try {
    if (fd < 0 || producer == nullptr || out == nullptr) {
      return bad_argument("invalid file descriptor or NULL argument");
    }
    *out = nullptr;
    auto handle = std::make_unique<obdeect_stream_writer>(fd);
    handle->output = &handle->fd_stream;
    return open_writer(std::move(handle), producer, out);
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_writer_set_int(obdeect_stream_writer *writer, const char *key, int64_t value) {
  try {
    return set_pending(writer, key, static_cast<std::int64_t>(value));
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_writer_set_float(obdeect_stream_writer *writer, const char *key, double value) {
  try {
    return set_pending(writer, key, value);
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_writer_set_string(obdeect_stream_writer *writer, const char *key, const char *value) {
  try {
    if (value == nullptr) {
      return bad_argument("NULL argument");
    }
    return set_pending(writer, key, std::string(value));
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_writer_set_array(obdeect_stream_writer *writer, const char *key, const double *values,
                              size_t n) {
  try {
    if (values == nullptr && n > 0) {
      return bad_argument("NULL array");
    }
    return set_pending(writer, key,
                       (n == 0 ? std::vector<double>() : std::vector<double>(values, values + n)));
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_writer_write_record(obdeect_stream_writer *writer, int record_type) {
  try {
    if (writer == nullptr) {
      return bad_argument("NULL argument");
    }
    obdeect::stream::Metadata metadata;
    metadata.swap(writer->pending);
    const auto type = static_cast<obdeect::stream::RecordType>(record_type);
    if (record_type < 0 || record_type > 0xFFFF || !obdeect::stream::is_metadata_record(type)) {
      return bad_argument("not a metadata record type");
    }
    writer->writer->write_metadata(type, metadata);
    return check_output(writer);
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_writer_write_photons(obdeect_stream_writer *writer, uint32_t telescope_id, uint32_t frame,
                                  uint64_t first_index, uint32_t rows,
                                  const float *const *columns) {
  try {
    if (writer == nullptr || (rows > 0 && columns == nullptr)) {
      return bad_argument("NULL argument");
    }
    if (frame > static_cast<std::uint32_t>(obdeect::stream::Frame::camera)) {
      return bad_argument("frame must be 0 (sky), 1 (ground), 2 (telescope) or 3 (camera)");
    }
    obdeect::stream::PhotonChunk &chunk = writer->chunk;
    std::uint32_t mask = 0;
    for (std::uint32_t i = 0; rows > 0 && i < obdeect::stream::n_columns; ++i) {
      if (columns[i] != nullptr) {
        mask |= 1U << i;
      }
    }
    chunk.resize(rows, mask); // throws for rows > 2^16
    chunk.telescope_id = telescope_id;
    chunk.frame = static_cast<obdeect::stream::Frame>(frame);
    chunk.first_index = first_index;
    for (std::uint32_t i = 0; i < obdeect::stream::n_columns; ++i) {
      if ((mask >> i) & 1U) {
        std::memcpy(chunk.columns[i].data(), columns[i], std::size_t{rows} * sizeof(float));
      }
    }
    writer->writer->write_chunk(chunk);
    return check_output(writer);
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_writer_write_blob(obdeect_stream_writer *writer, uint32_t eventio_type, uint32_t eventio_id,
                               const uint8_t *data, size_t size) {
  try {
    if (writer == nullptr || (data == nullptr && size > 0)) {
      return bad_argument("NULL argument");
    }
    writer->writer->write_blob(eventio_type, eventio_id,
                               size == 0 ? std::vector<std::uint8_t>()
                                         : std::vector<std::uint8_t>(data, data + size));
    return check_output(writer);
  } catch (...) {
    return from_exception();
  }
}

int obdeect_stream_writer_close(obdeect_stream_writer *writer) {
  if (writer == nullptr) {
    return OBDEECT_STREAM_OK;
  }
  int code = OBDEECT_STREAM_OK;
  try {
    writer->writer->finish();
    writer->output->flush();
    if (writer->file && writer->file->is_open()) {
      writer->file->close();
    }
    code = check_output(writer);
  } catch (...) {
    code = from_exception();
  }
  try {
    delete writer;
  } catch (...) {
  }
  return code;
}

} // extern "C"
