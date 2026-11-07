#pragma once

// EventIO block container (Bernloehr's eventio.c) and the IACT payloads of
// CORSIKA output (types 1200-1217, mc_tel.h). Own little-endian code; only the
// formats needed by gstrickd. Byte layout and bit meanings follow eventio.c
// (block description at its top, item header code at lines 4120-4190 for the
// reader and 3968-3998 for the writer) and io_simtel.c (payloads).

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace gstrickd {

// Sync tag bytes 37 8A 1F D4 read as little-endian (eventio.c:4133). The
// swapped value 0x378A1FD4 marks big-endian data (eventio.c:4135).
inline constexpr std::uint32_t eventio_sync = 0xD41F8A37U;
inline constexpr std::uint32_t eventio_sync_swapped = 0x378A1FD4U;

// IACT item types (mc_tel.h:IO_TYPE_MC_*).
namespace iact {
inline constexpr std::uint32_t runh = 1200;            // run header, float array
inline constexpr std::uint32_t telpos = 1201;          // telescope positions
inline constexpr std::uint32_t evth = 1202;            // event header, float array
inline constexpr std::uint32_t teloff = 1203;          // array offsets (and weights)
inline constexpr std::uint32_t telarray = 1204;        // container of 1205 blocks
inline constexpr std::uint32_t photons = 1205;         // photon bunches of one telescope
inline constexpr std::uint32_t evte = 1209;            // event end, float array
inline constexpr std::uint32_t rune = 1210;            // run end, float array
inline constexpr std::uint32_t inputcfg = 1212;        // CORSIKA input lines
inline constexpr std::uint32_t telarray_head = 1213;   // split mode: array start marker
inline constexpr std::uint32_t telarray_end = 1214;    // split mode: array end marker
inline constexpr std::uint32_t atmprof = 1216;         // atmospheric profile
inline constexpr std::uint32_t photons3d = 1217;       // 3D photon bunches (40 bytes each)
inline constexpr std::uint32_t compact_version = 1000; // photon version for compact bunches
} // namespace iact

struct EventioBlock {
  std::uint32_t type{};
  std::uint32_t version{};
  std::uint32_t ident{};
  bool user_flag{};
  bool only_subobjects{}; // length bit 30: the item consists of sub-items only
  // Bytes as read: sync tag (top-level blocks only) and item header, incl.
  // extended length when present. Empty for blocks built by make_block until
  // they are encoded.
  std::vector<std::uint8_t> header;
  std::vector<std::uint8_t> body; // data, including padding to 4 bytes

  // Raw bytes as stored in a file: header followed by body.
  [[nodiscard]] std::vector<std::uint8_t> raw() const {
    std::vector<std::uint8_t> out(header);
    out.insert(out.end(), body.begin(), body.end());
    return out;
  }
};

// Little-endian reader and writer for payloads.
class ByteReader {
public:
  ByteReader(const std::uint8_t *data, std::size_t size) : data_{data}, size_{size} {}
  explicit ByteReader(const std::vector<std::uint8_t> &bytes)
      : ByteReader(bytes.data(), bytes.size()) {}

  [[nodiscard]] std::uint64_t take(std::size_t n) {
    if (n > size_ - pos_) {
      throw std::runtime_error("EventIO: payload truncated");
    }
    std::uint64_t value = 0;
    for (std::size_t i = 0; i < n; ++i) {
      value |= std::uint64_t{data_[pos_ + i]} << (8 * i);
    }
    pos_ += n;
    return value;
  }
  [[nodiscard]] std::uint32_t u32() { return static_cast<std::uint32_t>(take(4)); }
  [[nodiscard]] std::int32_t i32() { return static_cast<std::int32_t>(take(4)); }
  [[nodiscard]] std::int16_t i16() { return static_cast<std::int16_t>(take(2)); }
  [[nodiscard]] float f32() { return std::bit_cast<float>(u32()); }
  [[nodiscard]] std::size_t remaining() const { return size_ - pos_; }

private:
  const std::uint8_t *data_;
  std::size_t size_;
  std::size_t pos_{};
};

class ByteWriter {
public:
  explicit ByteWriter(std::vector<std::uint8_t> &out) : out_{out} {}

  void put(std::uint64_t value, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
      out_.push_back(static_cast<std::uint8_t>(value >> (8 * i)));
    }
  }
  void u32(std::uint32_t value) { put(value, 4); }
  void i32(std::int32_t value) { put(static_cast<std::uint32_t>(value), 4); }
  void i16(std::int16_t value) { put(static_cast<std::uint16_t>(value), 2); }
  void f32(float value) { u32(std::bit_cast<std::uint32_t>(value)); }

private:
  std::vector<std::uint8_t> &out_;
};

namespace detail {

struct ItemHeader {
  std::uint32_t type{};
  std::uint32_t version{};
  std::uint32_t ident{};
  bool user_flag{};
  bool only_subobjects{};
  std::uint64_t length{};
  std::size_t size{}; // header bytes after the sync tag: 12, or 16 if extended
};

// Parse an item header (without sync tag) at p. Throws on truncation or on
// fields that eventio.c rejects.
inline ItemHeader parse_item_header(const std::uint8_t *p, std::size_t avail) {
  if (avail < 12) {
    throw std::runtime_error("EventIO: truncated item header");
  }
  ByteReader in{p, avail};
  const std::uint32_t type_version = in.u32();
  if ((type_version & 0x000C0000U) != 0) {
    throw std::runtime_error("EventIO: reserved bits set in type/version word");
  }
  ItemHeader h;
  h.type = type_version & 0xFFFFU;
  h.user_flag = (type_version & 0x00010000U) != 0;
  const bool extension_in_type = (type_version & 0x00020000U) != 0;
  h.version = (type_version >> 20) & 0xFFFU;
  if ((h.version & 0x800U) != 0) {
    throw std::runtime_error("EventIO: invalid version number (corrupt data)");
  }
  h.ident = in.u32();
  const std::uint32_t length_word = in.u32();
  h.only_subobjects = (length_word & 0x40000000U) != 0;
  // eventio.c:4170-4187 decodes the extension from length bit 31 (and writes
  // both flags, eventio.c:3994 and 3832). Files where the two flags disagree are
  // rejected.
  const bool extended = (length_word & 0x80000000U) != 0;
  if (extension_in_type != extended) {
    throw std::runtime_error("EventIO: inconsistent extension flags");
  }
  h.length = length_word & 0x3FFFFFFFU;
  h.size = 12;
  if (extended) {
    if (avail < 16) {
      throw std::runtime_error("EventIO: truncated item header");
    }
    const std::uint32_t extension = in.u32();
    if ((extension & ~0x0FFFU) != 0) {
      throw std::runtime_error("EventIO: reserved bits set in length extension");
    }
    h.length |= std::uint64_t{extension} << 30; // eventio.c:4184
    h.size = 16;
  }
  return h;
}

inline void copy_fields(EventioBlock &block, const ItemHeader &h) {
  block.type = h.type;
  block.version = h.version;
  block.ident = h.ident;
  block.user_flag = h.user_flag;
  block.only_subobjects = h.only_subobjects;
}

// Header bytes for a block, as eventio.c writes them. top_level adds the sync tag.
inline std::vector<std::uint8_t> encode_header(const EventioBlock &block, bool top_level) {
  if (block.type > 0xFFFFU) {
    throw std::invalid_argument("EventIO: type exceeds 16 bits");
  }
  if (block.version > 0x7FFU) {
    throw std::invalid_argument("EventIO: version exceeds 11 bits");
  }
  if (block.body.size() % 4 != 0) {
    throw std::invalid_argument("EventIO: body is not padded to 4 bytes");
  }
  const std::uint64_t length = block.body.size();
  if ((length >> 42) != 0) {
    throw std::invalid_argument("EventIO: body too long");
  }
  const bool extended = length >= (std::uint64_t{1} << 30);
  std::vector<std::uint8_t> out;
  ByteWriter w{out};
  if (top_level) {
    w.u32(eventio_sync);
  }
  std::uint32_t type_version = (block.type & 0xFFFFU) | (block.version << 20);
  if (block.user_flag) {
    type_version |= 0x00010000U;
  }
  if (extended) {
    type_version |= 0x00020000U;
  }
  w.u32(type_version);
  w.u32(block.ident);
  std::uint32_t length_word = static_cast<std::uint32_t>(length & 0x3FFFFFFFU);
  if (block.only_subobjects) {
    length_word |= 0x40000000U;
  }
  if (extended) {
    length_word |= 0x80000000U;
  }
  w.u32(length_word);
  if (extended) {
    w.u32(static_cast<std::uint32_t>(length >> 30));
  }
  return out;
}

} // namespace detail

// Sub-item form of a block (no sync tag), as stored inside a container.
inline std::vector<std::uint8_t> encode_item(const EventioBlock &block) {
  std::vector<std::uint8_t> out = detail::encode_header(block, false);
  out.insert(out.end(), block.body.begin(), block.body.end());
  return out;
}

// Reads top-level blocks. Little-endian data only.
class EventioReader {
public:
  explicit EventioReader(std::istream &in) : in_{in} {}

  // Read the next top-level block. Returns false at a clean end of stream.
  bool next(EventioBlock &block) {
    std::uint8_t sync[4]{};
    in_.read(reinterpret_cast<char *>(sync), 4);
    const std::streamsize got = in_.gcount();
    if (got == 0) {
      return false;
    }
    if (got != 4) {
      throw std::runtime_error("EventIO: truncated sync tag");
    }
    const std::uint32_t tag = ByteReader{sync, 4}.u32();
    if (tag == eventio_sync_swapped) {
      throw std::runtime_error("EventIO: big-endian data is not supported");
    }
    if (tag != eventio_sync) {
      throw std::runtime_error("EventIO: bad sync tag");
    }

    std::vector<std::uint8_t> item(12);
    read_exact(item.data(), 12);
    const bool extended = (ByteReader{item.data() + 8, 4}.u32() & 0x80000000U) != 0;
    if (extended) {
      item.resize(16);
      read_exact(item.data() + 12, 4);
    }
    const detail::ItemHeader h = detail::parse_item_header(item.data(), item.size());

    block = EventioBlock{};
    detail::copy_fields(block, h);
    block.header.assign(sync, sync + 4);
    block.header.insert(block.header.end(), item.begin(),
                        item.begin() + static_cast<std::ptrdiff_t>(h.size));
    // Read the body in chunks: a corrupt length must not allocate its size.
    block.body.clear();
    for (std::uint64_t left = h.length; left > 0;) {
      const auto chunk = static_cast<std::size_t>(std::min<std::uint64_t>(left, 1U << 20));
      const std::size_t old_size = block.body.size();
      block.body.resize(old_size + chunk);
      read_exact(block.body.data() + old_size, chunk);
      left -= chunk;
    }
    return true;
  }

private:
  void read_exact(std::uint8_t *dst, std::size_t n) {
    if (n == 0) {
      return;
    }
    in_.read(reinterpret_cast<char *>(dst), static_cast<std::streamsize>(n));
    if (static_cast<std::size_t>(in_.gcount()) != n) {
      throw std::runtime_error("EventIO: truncated block");
    }
  }

  std::istream &in_;
};

// Sub-blocks of a container block (for example 1205 inside 1204). The
// sub-blocks have no sync tag; their header holds 12 or 16 bytes.
inline std::vector<EventioBlock> subblocks(const EventioBlock &parent) {
  std::vector<EventioBlock> out;
  const std::vector<std::uint8_t> &data = parent.body;
  std::size_t pos = 0;
  while (pos < data.size()) {
    const detail::ItemHeader h = detail::parse_item_header(data.data() + pos, data.size() - pos);
    if (h.length > data.size() - pos - h.size) {
      throw std::runtime_error("EventIO: sub-block exceeds its container");
    }
    const auto start = static_cast<std::ptrdiff_t>(pos);
    const auto body_start = start + static_cast<std::ptrdiff_t>(h.size);
    const auto end = body_start + static_cast<std::ptrdiff_t>(h.length);
    EventioBlock block;
    detail::copy_fields(block, h);
    block.header.assign(data.begin() + start, data.begin() + body_start);
    block.body.assign(data.begin() + body_start, data.begin() + end);
    out.push_back(std::move(block));
    pos += h.size + static_cast<std::size_t>(h.length);
  }
  return out;
}

// Build a block. The body is padded with zero bytes to 4 bytes. The header is
// computed from the fields (length word, extension when the body is 1 GiB or
// more).
inline EventioBlock make_block(std::uint32_t type, std::uint32_t version, std::uint32_t ident,
                               std::vector<std::uint8_t> body, bool only_subobjects = false) {
  while (body.size() % 4 != 0) {
    body.push_back(0);
  }
  EventioBlock block;
  block.type = type;
  block.version = version;
  block.ident = ident;
  block.only_subobjects = only_subobjects;
  block.body = std::move(body);
  block.header = detail::encode_header(block, true);
  return block;
}

// Container block holding the given sub-blocks (length bit 30 set).
inline EventioBlock make_container(std::uint32_t type, std::uint32_t version, std::uint32_t ident,
                                   const std::vector<EventioBlock> &children) {
  std::vector<std::uint8_t> body;
  for (const EventioBlock &child : children) {
    const std::vector<std::uint8_t> item = encode_item(child);
    body.insert(body.end(), item.begin(), item.end());
  }
  return make_block(type, version, ident, std::move(body), true);
}

// Write a top-level block. The header is recomputed from the fields.
inline void write_block(std::ostream &out, const EventioBlock &block) {
  const std::vector<std::uint8_t> header = detail::encode_header(block, true);
  out.write(reinterpret_cast<const char *>(header.data()),
            static_cast<std::streamsize>(header.size()));
  out.write(reinterpret_cast<const char *>(block.body.data()),
            static_cast<std::streamsize>(block.body.size()));
  if (!out) {
    throw std::runtime_error("EventIO: write failed");
  }
}

// ---------------------------------------------------------------------------
// Float arrays: RUNH 1200, EVTH 1202, EVTE 1209, RUNE 1210 (write_tel_block,
// io_simtel.c:72). Body: int32 count, then count float32. The first value is
// stored as the bit pattern of a float (it holds the block name).

inline std::vector<float> decode_float_block(const EventioBlock &block) {
  ByteReader in{block.body};
  const std::int32_t count = in.i32();
  if (count < 0 || static_cast<std::size_t>(count) * 4 > in.remaining()) {
    throw std::runtime_error("EventIO: float block count does not match its size");
  }
  std::vector<float> values(static_cast<std::size_t>(count));
  for (float &v : values) {
    v = in.f32();
  }
  return values;
}

inline EventioBlock make_float_block(std::uint32_t type, std::uint32_t ident,
                                     const std::vector<float> &values) {
  std::vector<std::uint8_t> body;
  ByteWriter out{body};
  out.i32(static_cast<std::int32_t>(values.size()));
  for (const float v : values) {
    out.f32(v);
  }
  return make_block(type, 0, ident, std::move(body));
}

// ---------------------------------------------------------------------------
// Telescope positions, 1201 (write_tel_pos, io_simtel.c:484). Values in cm.

struct TelPos {
  std::vector<double> x, y, z, r;
};

inline TelPos decode_telpos(const EventioBlock &block) {
  if (block.version != 0) {
    throw std::runtime_error("EventIO: unsupported TELPOS version");
  }
  ByteReader in{block.body};
  const std::int32_t ntel = in.i32();
  if (ntel < 0 || static_cast<std::size_t>(ntel) * 16 > in.remaining()) {
    throw std::runtime_error("EventIO: TELPOS count does not match its size");
  }
  TelPos pos;
  for (std::vector<double> *column : {&pos.x, &pos.y, &pos.z, &pos.r}) {
    column->resize(static_cast<std::size_t>(ntel));
    for (double &v : *column) {
      v = in.f32();
    }
  }
  return pos;
}

inline EventioBlock make_telpos_block(const TelPos &pos) {
  const std::size_t n = pos.x.size();
  if (pos.y.size() != n || pos.z.size() != n || pos.r.size() != n) {
    throw std::invalid_argument("TELPOS columns differ in length");
  }
  std::vector<std::uint8_t> body;
  ByteWriter out{body};
  out.i32(static_cast<std::int32_t>(n));
  for (const std::vector<double> *column : {&pos.x, &pos.y, &pos.z, &pos.r}) {
    for (const double v : *column) {
      out.f32(static_cast<float>(v));
    }
  }
  return make_block(iact::telpos, 0, 0, std::move(body));
}

// ---------------------------------------------------------------------------
// Array offsets, 1203 (write_tel_offset_w, io_simtel.c:650). Version 1 adds
// one weight per array. Core offsets in cm, toff in ns.

struct TelOffset {
  double toff{};
  std::vector<double> xoff, yoff, weight;
  bool has_weight{};
};

inline TelOffset decode_teloff(const EventioBlock &block) {
  if (block.version > 1) {
    throw std::runtime_error("EventIO: unsupported TELOFF version");
  }
  ByteReader in{block.body};
  const std::int32_t narray = in.i32();
  const std::size_t columns = block.version == 1 ? 3 : 2;
  if (narray < 0 || static_cast<std::size_t>(narray) * 4 * columns + 4 > block.body.size()) {
    throw std::runtime_error("EventIO: TELOFF count does not match its size");
  }
  TelOffset off;
  off.has_weight = block.version == 1;
  off.toff = in.f32();
  const auto n = static_cast<std::size_t>(narray);
  off.xoff.resize(n);
  off.yoff.resize(n);
  for (double &v : off.xoff) {
    v = in.f32();
  }
  for (double &v : off.yoff) {
    v = in.f32();
  }
  if (off.has_weight) {
    off.weight.resize(n);
    for (double &v : off.weight) {
      v = in.f32();
    }
  }
  return off;
}

inline EventioBlock make_teloff_block(const TelOffset &off) {
  const std::size_t n = off.xoff.size();
  if (off.yoff.size() != n || (off.has_weight && off.weight.size() != n)) {
    throw std::invalid_argument("TELOFF columns differ in length");
  }
  std::vector<std::uint8_t> body;
  ByteWriter out{body};
  out.i32(static_cast<std::int32_t>(n));
  out.f32(static_cast<float>(off.toff));
  for (const double v : off.xoff) {
    out.f32(static_cast<float>(v));
  }
  for (const double v : off.yoff) {
    out.f32(static_cast<float>(v));
  }
  if (off.has_weight) {
    for (const double v : off.weight) {
      out.f32(static_cast<float>(v));
    }
  }
  return make_block(iact::teloff, off.has_weight ? 1 : 0, 0, std::move(body));
}

// ---------------------------------------------------------------------------
// Photon bunches, 1205 (write_tel_photons, io_simtel.c:1035; compact variant
// write_tel_compact_photons, io_simtel.c:1297). Header: int16 array, int16
// telescope, float32 photon sum, int32 bunch count. Units are CORSIKA's: cm,
// ns, cm a.s.l., nm.

struct Bunch {
  float x{}, y{}, cx{}, cy{}, ctime{}, zem{}, photons{}, lambda{};
};

struct TelPhotons {
  int array{};
  int tel{};       // 0-based index into TELPOS, as written in the block ident (1000*array+tel)
  float photons{}; // photon sum as stored in the header
  bool compact{};  // compact format (version 1000), decoded to the same units
  std::vector<Bunch> bunches;
};

inline TelPhotons decode_photons(const EventioBlock &block) {
  if (block.version % iact::compact_version != 0 || block.version / iact::compact_version > 1) {
    throw std::runtime_error("EventIO: unsupported photon block version");
  }
  ByteReader in{block.body};
  TelPhotons p;
  p.array = in.i16();
  p.tel = in.i16();
  p.photons = in.f32();
  p.compact = block.version == iact::compact_version;
  const std::int32_t n = in.i32();
  const std::size_t bytes_per_bunch = p.compact ? 16 : 32;
  if (n < 0 || static_cast<std::size_t>(n) * bytes_per_bunch > in.remaining()) {
    throw std::runtime_error("EventIO: photon count does not match block size");
  }
  p.bunches.resize(static_cast<std::size_t>(n));
  for (Bunch &b : p.bunches) {
    if (!p.compact) {
      b.x = in.f32();
      b.y = in.f32();
      b.cx = in.f32();
      b.cy = in.f32();
      b.ctime = in.f32();
      b.zem = in.f32();
      b.photons = in.f32();
      b.lambda = in.f32();
      continue;
    }
    // Compact decoding as io_simtel.c:1457-1479 (read_tel_photons).
    b.x = static_cast<float>(0.1 * in.i16());
    b.y = static_cast<float>(0.1 * in.i16());
    float cx = static_cast<float>(in.i16() / 30000.0);
    float cy = static_cast<float>(in.i16() / 30000.0);
    cx = cx > 1.0F ? 1.0F : (cx < -1.0F ? -1.0F : cx);
    cy = cy > 1.0F ? 1.0F : (cy < -1.0F ? -1.0F : cy);
    b.cx = cx;
    b.cy = cy;
    b.ctime = static_cast<float>(0.1 * in.i16());
    b.zem = static_cast<float>(std::pow(10.0, 0.001 * in.i16()));
    b.photons = static_cast<float>(0.01 * in.i16());
    b.lambda = static_cast<float>(in.i16());
  }
  return p;
}

// Encode photon bunches as 1205. Long format: version 0. Compact format:
// version 1000 with the quantisation of io_simtel.c; values that do not fit
// throw std::out_of_range.
inline EventioBlock make_photons_block(const TelPhotons &p, bool compact = false) {
  std::vector<std::uint8_t> body;
  ByteWriter out{body};
  out.i16(static_cast<std::int16_t>(p.array));
  out.i16(static_cast<std::int16_t>(p.tel));
  out.f32(p.photons);
  out.i32(static_cast<std::int32_t>(p.bunches.size()));
  const auto quantise = [](double value, double scale) {
    const double q = std::round(value * scale);
    if (!(q >= -32768.0 && q <= 32767.0)) {
      throw std::out_of_range("EventIO: compact photon value out of range");
    }
    return static_cast<std::int16_t>(q);
  };
  for (const Bunch &b : p.bunches) {
    if (!compact) {
      out.f32(b.x);
      out.f32(b.y);
      out.f32(b.cx);
      out.f32(b.cy);
      out.f32(b.ctime);
      out.f32(b.zem);
      out.f32(b.photons);
      out.f32(b.lambda);
      continue;
    }
    if (b.zem <= 0.0F) {
      throw std::out_of_range("EventIO: compact photon emission height must be positive");
    }
    out.i16(quantise(b.x, 10.0));
    out.i16(quantise(b.y, 10.0));
    out.i16(quantise(b.cx, 30000.0));
    out.i16(quantise(b.cy, 30000.0));
    out.i16(quantise(b.ctime, 10.0));
    out.i16(quantise(std::log10(b.zem), 1000.0));
    out.i16(quantise(b.photons, 100.0));
    out.i16(quantise(b.lambda, 1.0));
  }
  const std::uint32_t version = compact ? iact::compact_version : 0;
  const auto ident = static_cast<std::uint32_t>(1000 * p.array + p.tel);
  return make_block(iact::photons, version, ident, std::move(body));
}

// 3D photon bunches, 1217 (write_tel_photons3d, io_simtel.c:1164). Decoding only.
struct Bunch3d {
  float x{}, y{}, z{}, cx{}, cy{}, cz{}, ctime{}, dist{}, photons{}, lambda{};
};

struct TelPhotons3d {
  int array{};
  int tel{};
  float photons{};
  std::vector<Bunch3d> bunches;
};

inline TelPhotons3d decode_photons3d(const EventioBlock &block) {
  if (block.version != 0) {
    throw std::runtime_error("EventIO: unsupported 3D photon block version");
  }
  ByteReader in{block.body};
  TelPhotons3d p;
  p.array = in.i16();
  p.tel = in.i16();
  p.photons = in.f32();
  const std::int32_t n = in.i32();
  if (n < 0 || static_cast<std::size_t>(n) * 40 > in.remaining()) {
    throw std::runtime_error("EventIO: 3D photon count does not match block size");
  }
  p.bunches.resize(static_cast<std::size_t>(n));
  for (Bunch3d &b : p.bunches) {
    b.x = in.f32();
    b.y = in.f32();
    b.z = in.f32();
    b.cx = in.f32();
    b.cy = in.f32();
    b.cz = in.f32();
    b.ctime = in.f32();
    b.dist = in.f32();
    b.photons = in.f32();
    b.lambda = in.f32();
  }
  return p;
}

} // namespace gstrickd
