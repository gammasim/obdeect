#pragma once

// Photon-stream reader and writer following docs/stream/format.md. Header-only,
// standard library only. All multi-byte fields are little-endian; hosts of
// other byte orders convert per value.

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <istream>
#include <limits>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace obdeect::stream {

inline constexpr bool little_endian_host = std::endian::native == std::endian::little;

// Name of a record type as used in docs/stream/format.md and the Python API.
[[nodiscard]] constexpr const char *record_name(RecordType type) {
  switch (type) {
  case RecordType::stream:
    return "stream";
  case RecordType::run:
    return "run";
  case RecordType::event:
    return "event";
  case RecordType::telescope:
    return "telescope";
  case RecordType::photons:
    return "photons";
  case RecordType::telescope_end:
    return "telescope_end";
  case RecordType::event_end:
    return "event_end";
  case RecordType::run_end:
    return "run_end";
  case RecordType::blob:
    return "blob";
  case RecordType::efficiency_status:
    return "efficiency_status";
  case RecordType::stream_end:
    return "stream_end";
  }
  return "unknown";
}

// True for the record types whose payload is metadata (written by write_metadata).
[[nodiscard]] constexpr bool is_metadata_record(RecordType type) {
  switch (type) {
  case RecordType::run:
  case RecordType::event:
  case RecordType::telescope:
  case RecordType::telescope_end:
  case RecordType::event_end:
  case RecordType::run_end:
  case RecordType::efficiency_status:
    return true;
  default:
    return false;
  }
}

// Typed errors of the reader and writer. All derive from std::runtime_error, so generic catches
// keep working; the C API maps them by type (not by message).
// A stream that breaks the format: malformed content, or a record out of order.
class MalformedStreamError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};
// A stream that ends inside a record or a payload.
class TruncatedStreamError : public MalformedStreamError {
public:
  using MalformedStreamError::MalformedStreamError;
};
// A stream whose format major version differs from format_major.
class UnsupportedVersionError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};
// A read error of the input stream (not the end of the input).
class StreamReadError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

namespace stream_detail {

inline constexpr std::size_t record_header_size = 16;
inline constexpr std::size_t chunk_header_size = 24;
inline constexpr std::size_t blob_header_size = 8;
// Payloads of unknown size are read in steps, so a corrupt length cannot
// trigger a huge allocation before the input runs out.
inline constexpr std::size_t read_step = std::size_t{1} << 20;

inline constexpr std::uint8_t value_int = 1;
inline constexpr std::uint8_t value_float = 2;
inline constexpr std::uint8_t value_string = 3;
inline constexpr std::uint8_t value_array = 4;

[[noreturn]] inline void malformed(const char *what) {
  throw MalformedStreamError(std::string("obdeect stream: malformed stream: ") + what);
}

[[noreturn]] inline void truncated(const char *what) {
  throw TruncatedStreamError(std::string("obdeect stream: malformed stream: truncated ") + what);
}

inline void put_le(std::vector<std::uint8_t> &out, std::uint64_t value, std::size_t nbytes) {
  for (std::size_t i = 0; i < nbytes; ++i) {
    out.push_back(static_cast<std::uint8_t>(value >> (8 * i)));
  }
}

inline std::uint64_t load_le(const std::uint8_t *p, std::size_t nbytes) {
  std::uint64_t value = 0;
  for (std::size_t i = 0; i < nbytes; ++i) {
    value |= std::uint64_t{p[i]} << (8 * i);
  }
  return value;
}

// Bounds-checked sequential reader over one payload.
class Cursor {
public:
  Cursor(const std::uint8_t *data, std::size_t size) : next_(data), end_(data + size) {}

  [[nodiscard]] std::size_t remaining() const { return static_cast<std::size_t>(end_ - next_); }

  const std::uint8_t *take(std::size_t n) {
    if (n > remaining()) {
      truncated("payload");
    }
    const std::uint8_t *at = next_;
    next_ += n;
    return at;
  }

  std::uint64_t le(std::size_t nbytes) { return load_le(take(nbytes), nbytes); }

  std::string text(std::size_t n) {
    const std::uint8_t *at = take(n);
    return {reinterpret_cast<const char *>(at), n};
  }

private:
  const std::uint8_t *next_;
  const std::uint8_t *end_;
};

inline void append_metadata(std::vector<std::uint8_t> &out, const Metadata &metadata) {
  if (metadata.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::invalid_argument("metadata has too many entries");
  }
  put_le(out, metadata.size(), 4);
  for (const auto &[key, value] : metadata) {
    if (key.size() > std::numeric_limits<std::uint16_t>::max()) {
      throw std::invalid_argument("metadata key longer than 65535 bytes");
    }
    put_le(out, key.size(), 2);
    out.insert(out.end(), key.begin(), key.end());
    switch (value.index()) {
    case 0:
      out.push_back(value_int);
      put_le(out, static_cast<std::uint64_t>(std::get<std::int64_t>(value)), 8);
      break;
    case 1:
      out.push_back(value_float);
      put_le(out, std::bit_cast<std::uint64_t>(std::get<double>(value)), 8);
      break;
    case 2: {
      const auto &text = std::get<std::string>(value);
      if (text.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument("metadata string too long");
      }
      out.push_back(value_string);
      put_le(out, text.size(), 4);
      out.insert(out.end(), text.begin(), text.end());
      break;
    }
    default: {
      const auto &array = std::get<std::vector<double>>(value);
      if (array.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument("metadata array too long");
      }
      out.push_back(value_array);
      put_le(out, array.size(), 4);
      for (const double v : array) {
        put_le(out, std::bit_cast<std::uint64_t>(v), 8);
      }
      break;
    }
    }
  }
}

inline void decode_metadata_into(Metadata &metadata, const std::uint8_t *data, std::size_t size) {
  Cursor cursor(data, size);
  metadata.clear();
  const std::uint64_t count = cursor.le(4);
  for (std::uint64_t i = 0; i < count; ++i) {
    const auto key_size = static_cast<std::size_t>(cursor.le(2));
    std::string key = cursor.text(key_size);
    const auto type = static_cast<std::uint8_t>(cursor.le(1));
    Value value;
    switch (type) {
    case value_int:
      value = static_cast<std::int64_t>(cursor.le(8));
      break;
    case value_float:
      value = std::bit_cast<double>(cursor.le(8));
      break;
    case value_string: {
      const auto length = static_cast<std::size_t>(cursor.le(4));
      value = cursor.text(length);
      break;
    }
    case value_array: {
      const std::uint64_t n = cursor.le(4);
      if (n > cursor.remaining() / 8) {
        malformed("metadata array longer than its payload");
      }
      std::vector<double> array;
      array.reserve(static_cast<std::size_t>(n));
      for (std::uint64_t k = 0; k < n; ++k) {
        array.push_back(std::bit_cast<double>(cursor.le(8)));
      }
      value = std::move(array);
      break;
    }
    default:
      malformed("unknown metadata value type");
    }
    if (!metadata.emplace(std::move(key), std::move(value)).second) {
      malformed("duplicate metadata key");
    }
  }
  if (cursor.remaining() != 0) {
    malformed("trailing bytes after metadata");
  }
}

inline void load_floats(const std::uint8_t *src, float *dst, std::size_t n) {
  if constexpr (little_endian_host) {
    if (n > 0) {
      std::memcpy(dst, src, n * sizeof(float));
    }
  } else {
    for (std::size_t i = 0; i < n; ++i) {
      dst[i] = std::bit_cast<float>(static_cast<std::uint32_t>(load_le(src + 4 * i, 4)));
    }
  }
}

// Decodes and validates a photons record. With load_columns false the header and payload size are
// still checked and the row count is set, but no column is allocated or copied, so the photon
// values are not checked (keep_photons_of). path_allowed is false for streams of format 1.0 and
// 1.1, which have no path column; upward_allowed is false for format 1.0, which has no upward
// photons.
inline void decode_chunk(PhotonChunk &chunk, const std::uint8_t *data, std::size_t size,
                         bool load_columns, bool path_allowed, bool upward_allowed) {
  Cursor cursor(data, size);
  chunk.telescope_id = static_cast<std::uint32_t>(cursor.le(4));
  const std::uint64_t frame = cursor.le(4);
  if (frame > static_cast<std::uint32_t>(Frame::camera)) {
    malformed("unknown frame");
  }
  chunk.frame = static_cast<Frame>(frame);
  const auto mask = static_cast<std::uint32_t>(cursor.le(4));
  const std::uint64_t rows = cursor.le(4);
  chunk.first_index = cursor.le(8);
  if ((mask >> n_columns) != 0) {
    malformed("unknown column bit");
  }
  if (!path_allowed && (mask & column_bit(Column::path)) != 0) {
    malformed("path column in a stream of format 1.0 or 1.1");
  }
  if (rows > max_chunk_rows) {
    malformed("chunk exceeds 2^16 rows");
  }
  if (rows == 0 && mask != 0) {
    malformed("empty chunk with columns");
  }
  if (rows > 0 && (mask & required_columns) != required_columns) {
    malformed("chunk lacks a required column");
  }
  const auto present = static_cast<std::size_t>(std::popcount(mask));
  const auto n = static_cast<std::size_t>(rows);
  if (cursor.remaining() != n * present * sizeof(float)) {
    malformed("chunk payload size does not match its header");
  }
  if (!load_columns) {
    chunk.set_rows_without_columns(n);
    return;
  }
  chunk.resize(n, mask);
  for (std::uint32_t i = 0; i < n_columns; ++i) {
    if ((mask >> i) & 1U) {
      load_floats(cursor.take(n * sizeof(float)), chunk.columns[i].data(), n);
    }
  }
  try {
    validate_photon_values(chunk, upward_allowed);
  } catch (const std::invalid_argument &error) {
    malformed(error.what());
  }
}

// Efficiency status of one record: factor -> state (docs/stream/efficiency-status.md).
using EfficiencyStates = std::map<std::string, std::string>;

// Rank of a state in the order pending -> envelope -> applied; a factor may only advance.
inline int efficiency_rank(const std::string &state) {
  if (state == "pending") {
    return 0;
  }
  return state == "envelope" ? 1 : 2;
}

// The factors and states of one efficiency_status record: every factor named in `factors` needs
// `<name>.state` with one of the three states.
inline EfficiencyStates efficiency_states(const Metadata &metadata) {
  const auto factors = metadata.find("factors");
  if (factors == metadata.end()) {
    malformed("efficiency status lacks 'factors'");
  }
  const auto *list = std::get_if<std::string>(&factors->second);
  if (list == nullptr) {
    malformed("efficiency status 'factors' is not a string");
  }
  EfficiencyStates states;
  std::size_t begin = 0;
  while (begin <= list->size()) {
    std::size_t end = list->find(',', begin);
    if (end == std::string::npos) {
      end = list->size();
    }
    std::string name = list->substr(begin, end - begin);
    begin = end + 1;
    const auto first = name.find_first_not_of(" \t");
    if (first == std::string::npos) {
      continue;
    }
    name = name.substr(first, name.find_last_not_of(" \t") - first + 1);
    const auto state_it = metadata.find(name + ".state");
    if (state_it == metadata.end()) {
      malformed(("efficiency status lacks '" + name + ".state'").c_str());
    }
    const auto *state = std::get_if<std::string>(&state_it->second);
    if (state == nullptr || (*state != "applied" && *state != "envelope" && *state != "pending")) {
      malformed(("efficiency status factor '" + name +
                 "' has an unknown state (applied, envelope or pending)")
                    .c_str());
    }
    if (!states.emplace(name, *state).second) {
      malformed(("efficiency status lists factor '" + name + "' twice").c_str());
    }
  }
  return states;
}

// Record nesting and efficiency-status transitions (docs/stream/format.md, "Record";
// docs/stream/efficiency-status.md). The Reader runs accept() on every record it returns, so
// `stream cat`, the Python reader and the converters all reject the same streams.
class Nesting {
public:
  // Throws std::runtime_error (malformed) if a record of this type may not follow the records
  // seen so far. Unknown record types are not passed here.
  void accept(RecordType type, const Metadata &metadata) {
    switch (type) {
    case RecordType::run:
      if (run_open_) {
        malformed("run inside a run (run_end missing)");
      }
      run_open_ = true;
      events_ = 0;
      run_states_.clear();
      return;
    case RecordType::event:
      if (!run_open_) {
        malformed("event outside a run");
      }
      if (event_open_) {
        malformed("event inside an event (event_end missing)");
      }
      event_open_ = true;
      ++events_;
      return;
    case RecordType::telescope:
      if (!event_open_) {
        malformed("telescope outside an event");
      }
      if (telescope_open_) {
        malformed("telescope inside a telescope (telescope_end missing)");
      }
      telescope_open_ = true;
      telescope_photons_ = false;
      telescope_states_.clear();
      return;
    case RecordType::photons:
      if (!telescope_open_) {
        malformed("photons outside a telescope record");
      }
      telescope_photons_ = true;
      return;
    case RecordType::telescope_end:
      if (!telescope_open_) {
        malformed("telescope_end without telescope");
      }
      telescope_open_ = false;
      return;
    case RecordType::event_end:
      if (!event_open_) {
        malformed("event_end without event");
      }
      if (telescope_open_) {
        malformed("event_end inside a telescope (telescope_end missing)");
      }
      event_open_ = false;
      return;
    case RecordType::run_end:
      if (!run_open_) {
        malformed("run_end without run");
      }
      if (event_open_) {
        malformed("run_end inside an event (event_end missing)");
      }
      run_open_ = false;
      return;
    case RecordType::blob:
      if (!run_open_) {
        malformed("blob outside a run");
      }
      return;
    case RecordType::efficiency_status:
      efficiency_status(metadata);
      return;
    case RecordType::stream_end:
      if (run_open_) {
        malformed("stream_end inside a run (run_end missing)");
      }
      return;
    default:
      return; // the stream header
    }
  }

private:
  // Run status records come before the first event; a telescope status record comes after its
  // telescope record and before its photons. The two scopes are checked separately: the telescope
  // status overrides the run status per factor (docs/stream/efficiency-status.md).
  void efficiency_status(const Metadata &metadata) {
    if (!run_open_) {
      malformed("efficiency status outside a run");
    }
    if (telescope_open_) {
      if (telescope_photons_) {
        malformed("efficiency status after the photons of its telescope");
      }
      advance(telescope_states_, metadata, "telescope status records");
      return;
    }
    if (events_ > 0) {
      malformed("run efficiency status after the first event");
    }
    advance(run_states_, metadata, "run status records");
  }

  // Checks every factor first and then updates the states, so a rejected record changes nothing.
  static void advance(EfficiencyStates &states, const Metadata &metadata, const char *scope) {
    const EfficiencyStates record = efficiency_states(metadata);
    for (const auto &[name, state] : record) {
      const auto previous = states.find(name);
      if (previous != states.end() && efficiency_rank(state) <= efficiency_rank(previous->second)) {
        malformed(("efficiency status factor '" + name + "' in " + scope +
                   " has an invalid transition from " + previous->second + " to " + state)
                      .c_str());
      }
    }
    for (const auto &[name, state] : record) {
      states[name] = state;
    }
  }

  bool run_open_{false};
  bool event_open_{false};
  bool telescope_open_{false};
  bool telescope_photons_{false};
  std::uint64_t events_{0};
  EfficiencyStates run_states_;
  EfficiencyStates telescope_states_;
};

inline std::int64_t int_field(const Metadata &metadata, const char *key) {
  const auto it = metadata.find(key);
  if (it == metadata.end() || !std::holds_alternative<std::int64_t>(it->second)) {
    throw std::runtime_error(std::string("obdeect stream: stream record lacks integer key '") +
                             key + "'");
  }
  return std::get<std::int64_t>(it->second);
}

} // namespace stream_detail

// Metadata payload: u32 count, then per entry u16 key length, key, u8 type, value.
[[nodiscard]] inline std::vector<std::uint8_t> encode_metadata(const Metadata &metadata) {
  std::vector<std::uint8_t> out;
  stream_detail::append_metadata(out, metadata);
  return out;
}

// Throws std::runtime_error on malformed input.
[[nodiscard]] inline Metadata decode_metadata(const std::uint8_t *data, std::size_t size) {
  Metadata metadata;
  stream_detail::decode_metadata_into(metadata, data, size);
  return metadata;
}

// Writes records to a stream. The constructor writes the `stream` record;
// finish() writes `stream_end` (idempotent). Nothing is written after finish().
// A record that breaks the nesting rules of the Reader (for example an event outside a run, or
// a run still open at finish()) is refused with std::invalid_argument before any of its bytes
// are written, so the Writer cannot produce a stream the Reader would reject.
class Writer {
public:
  explicit Writer(std::ostream &out, std::string producer) : out_(out) {
    Metadata header;
    header["format_major"] = static_cast<std::int64_t>(format_major);
    header["format_minor"] = static_cast<std::int64_t>(format_minor);
    header["producer"] = std::move(producer);
    begin(RecordType::stream);
    stream_detail::append_metadata(buf_, header);
    finish_record();
  }

  void write_metadata(RecordType type, const Metadata &metadata) {
    if (!is_metadata_record(type)) {
      throw std::invalid_argument("write_metadata: not a metadata record type");
    }
    if (type == RecordType::run) {
      static_cast<void>(ground_representation(metadata));
    }
    accept_record(type, metadata);
    begin(type);
    stream_detail::append_metadata(buf_, metadata);
    finish_record();
  }

  // Validates the chunk, then writes its header and each present column as
  // contiguous float32 values (no copy on little-endian hosts).
  void write_chunk(const PhotonChunk &chunk) {
    check_open();
    chunk.validate();
    if (static_cast<std::uint32_t>(chunk.frame) > static_cast<std::uint32_t>(Frame::camera)) {
      throw std::invalid_argument("photon chunk has an unknown frame");
    }
    accept_record(RecordType::photons, no_metadata);
    const std::uint32_t mask = chunk.mask();
    const auto rows = static_cast<std::uint32_t>(chunk.rows());
    const auto present = static_cast<std::uint64_t>(std::popcount(mask));
    begin(RecordType::photons);
    stream_detail::put_le(buf_, chunk.telescope_id, 4);
    stream_detail::put_le(buf_, static_cast<std::uint32_t>(chunk.frame), 4);
    stream_detail::put_le(buf_, mask, 4);
    stream_detail::put_le(buf_, rows, 4);
    stream_detail::put_le(buf_, chunk.first_index, 8);
    set_length(stream_detail::chunk_header_size + std::uint64_t{rows} * present * 4U);
    flush_buffer();
    for (std::uint32_t i = 0; i < n_columns; ++i) {
      if ((mask >> i) & 1U) {
        write_floats(chunk.columns[i]);
      }
    }
  }

  void write_blob(std::uint32_t eventio_type, std::uint32_t eventio_id,
                  const std::vector<std::uint8_t> &data) {
    accept_record(RecordType::blob, no_metadata);
    begin(RecordType::blob);
    stream_detail::put_le(buf_, eventio_type, 4);
    stream_detail::put_le(buf_, eventio_id, 4);
    buf_.insert(buf_.end(), data.begin(), data.end());
    finish_record();
  }

  void finish() {
    if (finished_) {
      return;
    }
    accept_record(RecordType::stream_end, no_metadata);
    begin(RecordType::stream_end);
    finish_record();
    out_.flush();
    finished_ = true;
  }

private:
  void check_open() const {
    if (finished_) {
      throw std::logic_error("obdeect::stream::Writer: stream already finished");
    }
  }

  // Refuses a record that the Reader would refuse (the same Nesting rules), before any byte of
  // it is written. Throws std::invalid_argument.
  void accept_record(RecordType type, const Metadata &metadata) {
    check_open();
    try {
      nesting_.accept(type, metadata);
    } catch (const MalformedStreamError &error) {
      throw std::invalid_argument(std::string("obdeect::stream::Writer: ") + error.what());
    }
  }

  // Starts a record in buf_; the length field is set by set_length().
  void begin(RecordType type) {
    check_open();
    buf_.clear();
    stream_detail::put_le(buf_, record_magic, 4);
    stream_detail::put_le(buf_, static_cast<std::uint16_t>(type), 2);
    stream_detail::put_le(buf_, 0, 2);
    stream_detail::put_le(buf_, 0, 8);
  }

  void set_length(std::uint64_t payload_size) {
    for (std::size_t i = 0; i < 8; ++i) {
      buf_[8 + i] = static_cast<std::uint8_t>(payload_size >> (8 * i));
    }
  }

  void flush_buffer() {
    out_.write(reinterpret_cast<const char *>(buf_.data()),
               static_cast<std::streamsize>(buf_.size()));
  }

  // Completes a record whose whole payload is in buf_ and writes it in one call.
  void finish_record() {
    set_length(buf_.size() - stream_detail::record_header_size);
    flush_buffer();
  }

  void write_floats(const std::vector<float> &values) {
    if constexpr (little_endian_host) {
      out_.write(reinterpret_cast<const char *>(values.data()),
                 static_cast<std::streamsize>(values.size() * sizeof(float)));
    } else {
      buf_.clear();
      for (const float v : values) {
        stream_detail::put_le(buf_, std::bit_cast<std::uint32_t>(v), 4);
      }
      flush_buffer();
    }
  }

  std::ostream &out_;
  std::vector<std::uint8_t> buf_;
  bool finished_{false};
  stream_detail::Nesting nesting_;
  inline static const Metadata no_metadata{};
};

// Reads the next complete record (header and payload) as raw bytes into
// `bytes`. Returns false at EOF on a record boundary. Throws on a bad magic,
// nonzero flags or a truncated record.
inline bool read_record_raw(std::istream &in, std::vector<std::uint8_t> &bytes) {
  using namespace stream_detail;
  bytes.resize(record_header_size);
  in.read(reinterpret_cast<char *>(bytes.data()), static_cast<std::streamsize>(record_header_size));
  const auto got = static_cast<std::size_t>(in.gcount());
  if (got == 0) {
    if (in.bad()) {
      throw StreamReadError("obdeect stream: read error on input stream");
    }
    bytes.clear();
    return false;
  }
  if (got != record_header_size) {
    truncated("record header");
  }
  if (load_le(bytes.data(), 4) != record_magic) {
    throw std::runtime_error("obdeect stream: bad magic, not a photon stream");
  }
  if (load_le(bytes.data() + 6, 2) != 0) {
    malformed("nonzero reserved flags");
  }
  std::uint64_t remaining = load_le(bytes.data() + 8, 8);
  while (remaining > 0) {
    const auto step = static_cast<std::size_t>(std::min<std::uint64_t>(remaining, read_step));
    const std::size_t old_size = bytes.size();
    bytes.resize(old_size + step);
    in.read(reinterpret_cast<char *>(bytes.data() + old_size), static_cast<std::streamsize>(step));
    if (static_cast<std::size_t>(in.gcount()) != step) {
      truncated("record payload");
    }
    remaining -= step;
  }
  return true;
}

// Reads records from a stream written by Writer. Reuse one Record for many calls.
class Reader {
public:
  explicit Reader(std::istream &in) : in_(in) {}

  // Decode photon chunk columns only for these telescope ids. Other chunks are read, and their
  // header and payload size are checked, but their columns are not decoded and their values are not
  // checked. They come back with telescope_id, frame, first_index and rows() (the header row count)
  // set and no columns. Default: decode and check all.
  void keep_photons_of(std::vector<std::uint32_t> telescope_ids) {
    kept_telescopes_ = std::move(telescope_ids);
    filter_photons_ = true;
  }

  // Fills `record` with the next known record. Returns false after stream_end
  // (which is consumed, not returned) or at EOF on a record boundary. Unknown
  // record types are skipped. Throws std::runtime_error on malformed input or
  // a stream whose format_major differs from obdeect::stream::format_major.
  // Throws std::runtime_error on malformed input, on a stream whose format_major differs from
  // obdeect::stream::format_major, and on any record that breaks the nesting or efficiency status
  // rules. After an exception the reader has failed: every later call throws too.
  bool next(Record &record) {
    if (failed_) {
      throw std::runtime_error("obdeect stream: reader failed earlier");
    }
    try {
      return read_next(record);
    } catch (...) {
      failed_ = true;
      throw;
    }
  }

  // True once stream_end has been read.
  [[nodiscard]] bool ended() const { return ended_; }

  // True if stream_end has been read and the input is at EOF, i.e. nothing
  // (not even a partial record) follows it. Throws on a read error.
  [[nodiscard]] bool at_clean_end() {
    if (!ended_) {
      return false;
    }
    const auto next = in_.peek();
    if (in_.bad()) {
      throw StreamReadError("obdeect stream: read error on input stream");
    }
    return next == std::istream::traits_type::eof();
  }

private:
  bool read_next(Record &record) {
    using namespace stream_detail;
    while (!ended_) {
      if (!read_record_raw(in_, raw_)) {
        return false;
      }
      const auto type = static_cast<std::uint16_t>(load_le(raw_.data() + 4, 2));
      if (static_cast<RecordType>(type) == RecordType::stream) {
        if (header_seen_) {
          malformed("second stream record");
        }
      } else if (!header_seen_) {
        malformed("stream does not start with a stream record");
      }
      const std::uint8_t *payload = raw_.data() + record_header_size;
      const std::size_t size = raw_.size() - record_header_size;
      const auto rtype = static_cast<RecordType>(type);
      switch (rtype) {
      case RecordType::stream: {
        decode_metadata_into(record.metadata, payload, size);
        const std::int64_t major = int_field(record.metadata, "format_major");
        if (major != format_major) {
          throw UnsupportedVersionError("obdeect stream: unsupported format major version " +
                                        std::to_string(major) + " (reader supports " +
                                        std::to_string(format_major) + ")");
        }
        // Minor versions are compatible within a major version: a later minor is read too (its
        // unknown columns fail the column-mask check). A missing minor is 0 (format 1.0).
        const std::int64_t minor = record.metadata.contains("format_minor")
                                       ? int_field(record.metadata, "format_minor")
                                       : 0;
        if (minor < 0) {
          throw std::runtime_error("obdeect stream: negative format minor version " +
                                   std::to_string(minor));
        }
        stream_minor_ = static_cast<int>(minor);
        header_seen_ = true;
        record.type = RecordType::stream;
        return true;
      }
      case RecordType::run:
        decode_metadata_into(record.metadata, payload, size);
        // Rejects unknown representations; sphere runs need format 1.2 (path column, sphere rows).
        if (ground_representation(record.metadata) == GroundRepresentation::sphere &&
            stream_minor_ < 2) {
          malformed("ground_representation sphere in a stream of format 1.0 or 1.1");
        }
        break;
      case RecordType::event:
      case RecordType::telescope:
      case RecordType::telescope_end:
      case RecordType::event_end:
      case RecordType::run_end:
      case RecordType::efficiency_status:
        decode_metadata_into(record.metadata, payload, size);
        break;
      case RecordType::photons:
        // Payloads shorter than the telescope id fail in decode_chunk.
        decode_chunk(record.chunk, payload, size,
                     !filter_photons_ || size < 4 || kept(chunk_telescope(payload)),
                     stream_minor_ >= 2, stream_minor_ >= 1);
        break;
      case RecordType::blob:
        if (size < blob_header_size) {
          malformed("blob record shorter than its header");
        }
        record.blob_eventio_type = static_cast<std::uint32_t>(load_le(payload, 4));
        record.blob_eventio_id = static_cast<std::uint32_t>(load_le(payload + 4, 4));
        record.blob.assign(payload + blob_header_size, payload + size);
        break;
      case RecordType::stream_end:
        if (size != 0) {
          malformed("stream_end with payload");
        }
        nesting_.accept(rtype, record.metadata);
        ended_ = true;
        return false;
      default:
        continue; // unknown record type: skipped
      }
      nesting_.accept(rtype, record.metadata);
      record.type = rtype;
      return true;
    }
    return false;
  }

  static std::uint32_t chunk_telescope(const std::uint8_t *payload) {
    return static_cast<std::uint32_t>(stream_detail::load_le(payload, 4));
  }
  [[nodiscard]] bool kept(std::uint32_t telescope_id) const {
    return std::find(kept_telescopes_.begin(), kept_telescopes_.end(), telescope_id) !=
           kept_telescopes_.end();
  }

  std::istream &in_;
  std::vector<std::uint8_t> raw_;
  std::vector<std::uint32_t> kept_telescopes_;
  bool filter_photons_{false};
  int stream_minor_{0}; // format_minor of the stream record read last
  bool header_seen_{false};
  bool ended_{false};
  bool failed_{false};
  stream_detail::Nesting nesting_;
};

// Writes one record as returned by Reader (stream and stream_end are not written).
inline void write_record(Writer &writer, const Record &record) {
  switch (record.type) {
  case RecordType::stream:
  case RecordType::stream_end:
    return;
  case RecordType::photons:
    writer.write_chunk(record.chunk);
    return;
  case RecordType::blob:
    writer.write_blob(record.blob_eventio_type, record.blob_eventio_id, record.blob);
    return;
  default:
    writer.write_metadata(record.type, record.metadata);
    return;
  }
}

} // namespace obdeect::stream
