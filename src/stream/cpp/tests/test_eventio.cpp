#include "check.hpp"

#include "gstrickd/eventio.hpp"

#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace gstrickd;

namespace {

using Bytes = std::vector<std::uint8_t>;

Bytes to_bytes(const std::string &text) { return Bytes(text.begin(), text.end()); }

// Little-endian word, for hand-built headers.
void put_u32(Bytes &out, std::uint32_t value) {
  for (int i = 0; i < 4; ++i) {
    out.push_back(static_cast<std::uint8_t>(value >> (8 * i)));
  }
}

template <class Action> bool throws(Action action) {
  try {
    action();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

void header_round_trip() {
  const EventioBlock block = make_block(1204, 0, 7, Bytes(8, 0x5A), true);
  CHECK(block.header.size() == 16);
  CHECK(block.header[0] == 0x37 && block.header[1] == 0x8A && block.header[2] == 0x1F &&
        block.header[3] == 0xD4); // sync tag, little-endian
  std::ostringstream out;
  write_block(out, block);
  const std::string text = out.str();
  CHECK(text.size() == 24);

  std::istringstream in{text};
  EventioReader reader{in};
  EventioBlock read;
  CHECK(reader.next(read));
  CHECK(read.type == 1204);
  CHECK(read.version == 0);
  CHECK(read.ident == 7);
  CHECK(read.only_subobjects);
  CHECK(!read.user_flag);
  CHECK(read.body == block.body);
  CHECK(read.raw() == to_bytes(text));
  EventioBlock none;
  CHECK(!reader.next(none)); // clean end of stream
}

void user_flag_and_version_round_trip() {
  EventioBlock block = make_block(1210, 3, 0xFFFFFFFFU, Bytes(4, 1));
  block.user_flag = true;
  block.header = detail::encode_header(block, true);
  std::ostringstream out;
  write_block(out, block);
  std::istringstream in{out.str()};
  EventioReader reader{in};
  EventioBlock read;
  CHECK(reader.next(read));
  CHECK(read.type == 1210 && read.version == 3 && read.ident == 0xFFFFFFFFU);
  CHECK(read.user_flag && !read.only_subobjects);
  CHECK(read.raw() == to_bytes(out.str()));
}

// Extension: type bit 17 and length bit 31 both set, then an extension word
// (length bits 30-41). eventio.c:4170-4184.
void extended_length_is_read() {
  Bytes raw;
  put_u32(raw, eventio_sync);
  put_u32(raw, 1205U | 0x00020000U); // type 1205, version 0, extension flag
  put_u32(raw, 3);                   // ident
  put_u32(raw, 8U | 0x80000000U);    // length 8, extension flag
  put_u32(raw, 0);                   // extension word: length bits 30-41 = 0
  for (std::uint8_t i = 0; i < 8; ++i) {
    raw.push_back(i);
  }
  std::istringstream in{std::string(raw.begin(), raw.end())};
  EventioReader reader{in};
  EventioBlock read;
  CHECK(reader.next(read));
  CHECK(read.type == 1205 && read.ident == 3);
  CHECK(read.header.size() == 20);
  CHECK(read.body.size() == 8 && read.body[7] == 7);
  CHECK(read.raw() == raw); // header kept as read
}

void extension_length_bits_are_decoded() {
  // Extension word 1 means 2^30 bytes of body: a container of 8 bytes cannot
  // hold it, so both the reader and sub-block parsing must reject it.
  Bytes sub;
  put_u32(sub, 1205U | 0x00020000U);
  put_u32(sub, 0);
  put_u32(sub, 8U | 0x80000000U);
  put_u32(sub, 1U);
  for (int i = 0; i < 8; ++i) {
    sub.push_back(0);
  }
  EventioBlock container;
  container.type = 1204;
  container.body = sub;
  CHECK(throws([&] { (void)subblocks(container); }));

  Bytes top;
  put_u32(top, eventio_sync);
  put_u32(top, 1205U | 0x00020000U);
  put_u32(top, 0);
  put_u32(top, 8U | 0x80000000U);
  put_u32(top, 1U);
  for (int i = 0; i < 8; ++i) {
    top.push_back(0);
  }
  std::istringstream in{std::string(top.begin(), top.end())};
  EventioReader reader{in};
  EventioBlock read;
  CHECK(throws([&] { (void)reader.next(read); }));
}

void rejects_bad_input() {
  // Big-endian sync tag (bytes D4 1F 8A 37).
  std::istringstream big{std::string("\xD4\x1F\x8A\x37", 4)};
  EventioBlock read;
  CHECK(throws([&] {
    EventioReader r{big};
    (void)r.next(read);
  }));
  std::istringstream bad{std::string(4, '\0')};
  CHECK(throws([&] {
    EventioReader r{bad};
    (void)r.next(read);
  }));
  std::istringstream sync_only{std::string("\x37\x8A\x1F\xD4", 4)};
  CHECK(throws([&] {
    EventioReader r{sync_only};
    (void)r.next(read);
  }));
  // Body shorter than the length word claims.
  Bytes short_body;
  put_u32(short_body, eventio_sync);
  put_u32(short_body, 1205);
  put_u32(short_body, 0);
  put_u32(short_body, 8);
  short_body.push_back(0);
  std::istringstream truncated{std::string(short_body.begin(), short_body.end())};
  CHECK(throws([&] {
    EventioReader r{truncated};
    (void)r.next(read);
  }));
  // Version with bit 11 set is invalid (eventio.c:4156).
  Bytes bad_version;
  put_u32(bad_version, eventio_sync);
  put_u32(bad_version, 1205U | (0x800U << 20));
  put_u32(bad_version, 0);
  put_u32(bad_version, 0);
  std::istringstream bv{std::string(bad_version.begin(), bad_version.end())};
  CHECK(throws([&] {
    EventioReader r{bv};
    (void)r.next(read);
  }));
}

void subblocks_of_container() {
  const EventioBlock a = make_block(1205, 0, 1, Bytes(5, 0xA1)); // padded to 8
  const EventioBlock b = make_block(1205, 0, 2, Bytes(4, 0xB2));
  const EventioBlock container = make_container(1204, 0, 0, {a, b});
  CHECK(container.only_subobjects);
  CHECK(container.body.size() == (12 + 8) + (12 + 4)); // sub-item: header + padded body
  const std::vector<EventioBlock> subs = subblocks(container);
  CHECK(subs.size() == 2);
  if (subs.size() == 2) {
    CHECK(subs[0].type == 1205 && subs[0].ident == 1);
    CHECK(subs[0].body == a.body);
    CHECK(subs[0].header.size() == 12);
    CHECK(subs[1].ident == 2 && subs[1].body == b.body);
    CHECK(subs[1].raw() == encode_item(b));
  }
}

void float_block_round_trip() {
  const std::vector<float> values{1.0F, 2.5F, -3.0F, 1e7F};
  const EventioBlock block = make_float_block(1202, 42, values);
  CHECK(block.type == 1202 && block.ident == 42);
  CHECK(block.body.size() == 4 + 4 * values.size());
  CHECK(decode_float_block(block) == values);
  EventioBlock broken = block;
  broken.body[0] = 9; // count larger than the data
  CHECK(throws([&] { (void)decode_float_block(broken); }));
}

void telpos_round_trip() {
  TelPos pos;
  pos.x = {0.0, 1200.0, -500.0};
  pos.y = {0.0, -300.0, 800.0};
  pos.z = {1800.0, 1900.0, 2000.0};
  pos.r = {1200.0, 1200.0, 1200.0};
  const TelPos back = decode_telpos(make_telpos_block(pos));
  CHECK(back.x == pos.x && back.y == pos.y && back.z == pos.z && back.r == pos.r);
}

void teloff_round_trip() {
  TelOffset off;
  off.toff = 12.5;
  off.xoff = {0.0, 1500.0};
  off.yoff = {-300.0, 0.0};
  off.weight = {0.5, 0.5};
  off.has_weight = true;
  const EventioBlock block = make_teloff_block(off);
  CHECK(block.version == 1);
  const TelOffset back = decode_teloff(block);
  CHECK(back.has_weight);
  CHECK(back.toff == 12.5);
  CHECK(back.xoff == off.xoff && back.yoff == off.yoff && back.weight == off.weight);

  off.has_weight = false;
  off.weight.clear();
  const EventioBlock plain = make_teloff_block(off);
  CHECK(plain.version == 0);
  CHECK(!decode_teloff(plain).has_weight);
}

void long_photons_layout_and_round_trip() {
  TelPhotons p;
  p.array = 2;
  p.tel = 5;
  p.photons = 3.0F;
  p.bunches = {{100.0F, -200.0F, 0.1F, 0.2F, 35.5F, 2.2e6F, 1.5F, 0.0F},
               {1200.0F, 300.0F, -0.1F, 0.05F, 36.25F, 2.2e6F, 1.5F, 350.0F}};
  const EventioBlock block = make_photons_block(p);
  CHECK(block.type == 1205 && block.version == 0);
  CHECK(block.ident == 2005U);
  // short array, short tel, float photons, long n, 8 floats per bunch.
  CHECK(block.body.size() == 2 + 2 + 4 + 4 + 2 * 32);
  CHECK(block.body[0] == 2 && block.body[1] == 0 && block.body[2] == 5 && block.body[3] == 0);
  const TelPhotons back = decode_photons(block);
  CHECK(back.array == 2 && back.tel == 5 && back.photons == 3.0F && !back.compact);
  CHECK(back.bunches.size() == 2);
  if (back.bunches.size() == 2) {
    CHECK(back.bunches[1].x == 1200.0F && back.bunches[1].lambda == 350.0F);
    CHECK(back.bunches[0].zem == 2.2e6F && back.bunches[0].ctime == 35.5F);
  }
}

// Compact bunches: decoded with the formulas of io_simtel.c:1457-1479. The raw
// values are written by hand, not by the encoder.
void compact_decoding_matches_io_simtel() {
  Bytes body;
  ByteWriter out{body};
  out.i16(1);
  out.i16(2);
  out.f32(3.0F);
  out.i32(2);
  // bunch A
  out.i16(123);    // x: 0.1 * 123 cm
  out.i16(-45);    // y
  out.i16(15000);  // cx = 15000 / 30000
  out.i16(-30000); // cy = -1
  out.i16(-7);     // ctime: 0.1 * -7 ns
  out.i16(6342);   // zem = 10^(0.001 * 6342) cm
  out.i16(150);    // photons: 0.01 * 150
  out.i16(0);      // lambda
  // bunch B: cx and cy beyond +-1 are clamped
  out.i16(0);
  out.i16(0);
  out.i16(32767);
  out.i16(-32768);
  out.i16(0);
  out.i16(0);
  out.i16(1);
  out.i16(9999);
  const EventioBlock block = make_block(1205, 1000, 1002, body);
  const TelPhotons p = decode_photons(block);
  CHECK(p.compact && p.array == 1 && p.tel == 2 && p.photons == 3.0F);
  CHECK(p.bunches.size() == 2);
  if (p.bunches.size() == 2) {
    const Bunch &a = p.bunches[0];
    CHECK(a.x == static_cast<float>(0.1 * 123));
    CHECK(a.y == static_cast<float>(0.1 * -45));
    CHECK(a.cx == static_cast<float>(15000 / 30000.0));
    CHECK(a.cy == -1.0F);
    CHECK(a.ctime == static_cast<float>(0.1 * -7));
    CHECK(a.zem == static_cast<float>(std::pow(10.0, 0.001 * 6342)));
    CHECK(a.photons == static_cast<float>(0.01 * 150));
    CHECK(a.lambda == 0.0F);
    const Bunch &b = p.bunches[1];
    CHECK(b.cx == 1.0F);
    CHECK(b.cy == -1.0F);
    CHECK(b.photons == static_cast<float>(0.01 * 1));
    CHECK(b.lambda == 9999.0F);
  }
}

void compact_encoding_quantises() {
  TelPhotons p;
  p.array = 0;
  p.tel = 1;
  p.photons = 2.0F;
  p.bunches = {{12.34F, -4.56F, 0.25F, -0.5F, 40.0F, 2.2e6F, 1.0F, 0.0F}};
  const TelPhotons back = decode_photons(make_photons_block(p, true));
  CHECK(back.compact);
  CHECK(back.bunches.size() == 1);
  if (back.bunches.size() == 1) {
    const Bunch &b = back.bunches[0];
    CHECK_NEAR(b.x, 12.34, 0.05);
    CHECK_NEAR(b.y, -4.56, 0.05);
    CHECK_NEAR(b.cx, 0.25, 1e-4);
    CHECK_NEAR(b.zem, 2.2e6, 2.2e6 * 0.0024);
    CHECK_NEAR(b.ctime, 40.0, 0.05);
  }
  // A bunch with x = 10000 cm needs 100000 mm: out of the int16 range.
  p.bunches = {{10000.0F, 0.0F, 0.0F, 0.0F, 0.0F, 2.2e6F, 1.0F, 0.0F}};
  CHECK(throws([&] { (void)make_photons_block(p, true); }));
}

void photons3d_decoding() {
  Bytes body;
  ByteWriter out{body};
  out.i16(0);
  out.i16(3);
  out.f32(1.0F);
  out.i32(1);
  const float values[10] = {1.0F, 2.0F, 3.0F, 0.0F, 0.0F, -1.0F, 35.0F, 100.0F, 1.0F, 500.0F};
  for (const float v : values) {
    out.f32(v);
  }
  const TelPhotons3d p = decode_photons3d(make_block(1217, 0, 3, body));
  CHECK(p.tel == 3 && p.bunches.size() == 1);
  if (p.bunches.size() == 1) {
    CHECK(p.bunches[0].z == 3.0F && p.bunches[0].cz == -1.0F);
    CHECK(p.bunches[0].dist == 100.0F && p.bunches[0].lambda == 500.0F);
  }
}

void photon_count_is_checked() {
  Bytes body;
  ByteWriter out{body};
  out.i16(0);
  out.i16(0);
  out.f32(1.0F);
  out.i32(5); // five bunches claimed, none present
  CHECK(throws([&] { (void)decode_photons(make_block(1205, 0, 0, body)); }));
}

// Array and telescope indices are int16 on disk: only 0..32767 is encoded.
void photon_block_indices_are_range_checked() {
  TelPhotons p;
  p.photons = 1.0F;
  p.array = -1;
  CHECK(throws([&] { (void)make_photons_block(p); }));
  p.array = 32768;
  CHECK(throws([&] { (void)make_photons_block(p); }));
  p.array = 0;
  p.tel = -1;
  CHECK(throws([&] { (void)make_photons_block(p); }));
  p.tel = 32768;
  CHECK(throws([&] { (void)make_photons_block(p); }));
  p.array = 32767;
  p.tel = 32767;
  const TelPhotons back = decode_photons(make_photons_block(p));
  CHECK(back.array == 32767 && back.tel == 32767);
}

} // namespace

// Variable-length counts: eventio.c byte layout at the branch boundaries, and round trips.
void count_encoding() {
  Bytes out;
  ByteWriter w{out};
  put_count(w, 128);
  CHECK(out == Bytes({0x80, 0x80}));
  out.clear();
  put_count(w, 16384);
  CHECK(out == Bytes({0xC0, 0x40, 0x00}));
  for (const std::uint64_t n :
       {std::uint64_t{0}, std::uint64_t{127}, std::uint64_t{128}, std::uint64_t{16383},
        std::uint64_t{16384}, (std::uint64_t{1} << 21) - 1, std::uint64_t{1} << 21,
        (std::uint64_t{1} << 28) - 1, std::uint64_t{1} << 28, (std::uint64_t{1} << 32) - 1}) {
    out.clear();
    put_count(w, n);
    ByteReader r{out};
    CHECK(take_count(r) == n);
    CHECK(r.remaining() == 0);
  }
  out.clear();
  CHECK(throws([&] { put_count(w, std::uint64_t{1} << 32); }));
}

// ATMPROF (1216) block: layout, file round trip through EventioReader, and a rejected 5-layer
// block (the writer never produces one).
void atmprof_round_trip() {
  AtmosphereTable table;
  table.profile_id = 99;
  table.file_name = "atmprof1.dat";
  table.observation_level_cm = 220000.0;
  table.alt_km = {0.0, 1.0, 5.0};
  table.rho = {1.2, 1.0, 0.5};
  table.thick = {1000.0, 900.0, 500.0};
  table.refidx_m1 = {2.9e-4, 2.7e-4, 2.2e-4};
  const EventioBlock block = make_atmprof_block(table);
  CHECK(block.type == iact::atmprof);
  CHECK(block.version == 1);
  CHECK(block.ident == 99);
  CHECK(block.body.size() == 120); // 13 (name) + 8 (level) + 1 (count) + 96 (levels) + 1 + 1 pad
  CHECK(block.body[0] == 12);      // name length, then the name
  CHECK(block.body[21] == 3);      // level count after the level

  std::ostringstream out;
  write_block(out, block);
  std::istringstream in(out.str());
  EventioReader reader{in};
  EventioBlock read;
  CHECK(reader.next(read));
  const AtmosphereTable back = decode_atmprof(read);
  CHECK(back.profile_id == 99);
  CHECK(back.file_name == table.file_name);
  CHECK(back.observation_level_cm == table.observation_level_cm);
  CHECK(back.alt_km == table.alt_km);
  CHECK(back.rho == table.rho);
  CHECK(back.thick == table.thick);
  CHECK(back.refidx_m1 == table.refidx_m1);

  AtmosphereTable bad = table;
  bad.rho.pop_back();
  CHECK(throws([&] { make_atmprof_block(bad); }));

  Bytes layers;
  ByteWriter w{layers};
  put_count(w, 0);
  w.f64(1.0);
  put_count(w, 0);
  put_count(w, 5); // 5-layer parameters: not read
  CHECK(throws([&] { decode_atmprof(make_block(iact::atmprof, 1, 99, layers)); }));
  CHECK(throws([&] { decode_atmprof(make_block(iact::atmprof, 0, 99, layers)); }));
}

int main() {
  return run_checks([] {
    header_round_trip();
    user_flag_and_version_round_trip();
    extended_length_is_read();
    extension_length_bits_are_decoded();
    rejects_bad_input();
    subblocks_of_container();
    float_block_round_trip();
    telpos_round_trip();
    teloff_round_trip();
    long_photons_layout_and_round_trip();
    compact_decoding_matches_io_simtel();
    compact_encoding_quantises();
    photons3d_decoding();
    photon_count_is_checked();
    photon_block_indices_are_range_checked();
    count_encoding();
    atmprof_round_trip();
  });
}
