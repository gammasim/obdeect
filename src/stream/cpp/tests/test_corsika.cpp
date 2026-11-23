#include "check.hpp"

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/corsika.hpp"
#include "obdeect/stream/simtel.hpp"
#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace obdeect::stream;

namespace {

using Bytes = std::vector<std::uint8_t>;

template <class Action> bool throws_error(Action action) {
  try {
    action();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

// A small CORSIKA-like IACT file built with make_block. Stream A: one event in
// normal mode, two long-format telescopes. Stream B adds a compact telescope to
// event 1 and a second event in split mode (1213, 1205, 1214).
EventioBlock runh_block() {
  std::vector<float> w(273, 0.0F);
  w[1] = 12345.0F;   // run number
  w[3] = 7.41F;      // version
  w[4] = 1.0F;       // number of observation levels
  w[5] = 2200000.0F; // observation level 1 [cm]
  w[15] = -2.0F;     // spectral slope
  w[16] = 10.0F;     // energy min [GeV]
  w[17] = 100000.0F; // energy max [GeV]
  return make_float_block(iact::runh, 12345, w);
}

EventioBlock evth_block(float event, float curved_flag = 0.0F) {
  std::vector<float> w(273, 0.0F);
  w[1] = event;
  w[2] = 14.0F;        // gamma
  w[3] = 1000.0F;      // energy [GeV]
  w[6] = 2500000.0F;   // first interaction at 25 km [cm]
  w[10] = 0.3F;        // zenith [rad]
  w[11] = 1.2F;        // azimuth [rad]
  w[76] = 5136.0F;     // atmosphere 5 (bits 11-20), refraction (bit 5)
  w[78] = curved_flag; // word 79: CURVED flag (0 = standard, 1 = curved, 2 = curved with upgoing)
  w[95] = 250.0F;      // Cherenkov wavelength range [nm]
  w[96] = 800.0F;
  return make_float_block(iact::evth, static_cast<std::uint32_t>(event), w);
}

EventioBlock evte_block(float event) {
  std::vector<float> w(273, 0.0F);
  w[1] = event;
  return make_float_block(iact::evte, static_cast<std::uint32_t>(event), w);
}

EventioBlock rune_block() {
  std::vector<float> w(3, 0.0F);
  w[1] = 12345.0F;
  w[2] = 2.0F;
  return make_float_block(iact::rune, 12345, w);
}

EventioBlock telpos_block() {
  TelPos pos;
  pos.x = {0.0, 1200.0, -500.0};
  pos.y = {0.0, -300.0, 800.0};
  pos.z = {1800.0, 1900.0, 2000.0};
  pos.r = {1200.0, 1200.0, 1200.0};
  return make_telpos_block(pos);
}

EventioBlock teloff_block() {
  TelOffset off;
  off.toff = 12.5;
  off.xoff = {0.0, 1500.0};
  off.yoff = {-300.0, 0.0};
  off.weight = {0.5, 0.5};
  off.has_weight = true;
  return make_teloff_block(off);
}

EventioBlock long_block(int tel, float photons, std::vector<Bunch> bunches) {
  TelPhotons p;
  p.array = 0;
  p.tel = tel;
  p.photons = photons;
  p.bunches = std::move(bunches);
  return make_photons_block(p);
}

// Bunch positions are multiples of 100 cm, so that cm -> m -> cm is exact.
const std::vector<Bunch> tel0_bunches{{100.0F, -200.0F, 0.1F, 0.2F, 35.5F, 2.2e6F, 1.5F, 0.0F},
                                      {1200.0F, 300.0F, -0.1F, 0.05F, 36.25F, 2.2e6F, 1.5F, 0.0F}};
const std::vector<Bunch> tel1_bunches{{-300.0F, 0.0F, 0.0F, 0.0F, 40.0F, 2.2e6F, 2.0F, 350.0F}};

std::string serialize(const std::vector<EventioBlock> &blocks) {
  std::ostringstream out;
  for (const EventioBlock &block : blocks) {
    write_block(out, block);
  }
  return out.str();
}

std::string build_stream_a(float curved_flag = 0.0F) {
  return serialize(
      {runh_block(), make_block(iact::inputcfg, 0, 0, Bytes{'R', 'U', 'N', '\n'}), telpos_block(),
       evth_block(1, curved_flag), teloff_block(),
       make_container(iact::telarray, 0, 0,
                      {long_block(0, 3.0F, tel0_bunches), long_block(1, 2.0F, tel1_bunches)}),
       evte_block(1), rune_block()});
}

std::string build_stream_b() {
  // Compact telescope 2 in event 1: raw values chosen so that the quantised
  // values are exact enough to check.
  TelPhotons compact;
  compact.array = 0;
  compact.tel = 2;
  compact.photons = 3.0F;
  compact.bunches = {{12.3F, -4.5F, 0.5F, 0.25F, 40.0F, 2.2e6F, 1.5F, 0.0F},
                     {-1200.0F, 0.0F, 0.0F, -0.5F, 41.0F, 1.0e5F, 1.5F, 0.0F}};
  const EventioBlock compact_block = make_photons_block(compact, true);
  TelPhotons split;
  split.array = 1;
  split.tel = 0;
  split.photons = 1.0F;
  split.bunches = {{5.0F, 6.0F, 0.0F, 0.0F, 10.0F, 2.2e6F, 1.0F, 0.0F}};
  const EventioBlock split_photons = make_photons_block(split);
  return serialize({runh_block(), make_block(iact::inputcfg, 0, 0, Bytes{'R', 'U', 'N', '\n'}),
                    telpos_block(), evth_block(1), teloff_block(),
                    make_container(iact::telarray, 0, 0,
                                   {long_block(0, 3.0F, tel0_bunches),
                                    long_block(1, 2.0F, tel1_bunches), compact_block}),
                    evte_block(1), evth_block(2), teloff_block(),
                    make_block(iact::telarray_head, 0, 1, {}), split_photons,
                    make_block(iact::telarray_end, 0, 1, {}), evte_block(2), rune_block()});
}

std::vector<Record> convert(const std::string &bytes) {
  std::istringstream in{bytes};
  std::vector<Record> records;
  corsika_to_records(in, [&](const Record &r) { records.push_back(r); });
  return records;
}

std::vector<RecordType> types_of(const std::vector<Record> &records) {
  std::vector<RecordType> types;
  for (const Record &r : records) {
    types.push_back(r.type);
  }
  return types;
}

template <class T> T meta(const Metadata &m, const std::string &key) {
  return std::get<T>(m.at(key));
}

void units_are_double_then_float() {
  // x_m = float(double(x_cm) / 100), as documented in corsika.hpp.
  const float x_cm = 12345.678F;
  CHECK(corsika_cm_to_m(x_cm) == static_cast<float>(static_cast<double>(x_cm) / 100.0));
  CHECK(corsika_m_to_cm(1.0F) == 100.0F);
}

// Round trip cm -> m -> cm: never more than 1 ulp off, not always exact.
void float_round_trip_is_within_one_ulp() {
  std::mt19937_64 rng{12345};
  std::uniform_real_distribution<double> uniform{-1.0e5, 1.0e5};
  std::uint64_t exact = 0;
  std::uint64_t off_by_one = 0;
  std::uint64_t worse = 0;
  constexpr std::uint64_t n = 1000000;
  for (std::uint64_t i = 0; i < n; ++i) {
    const auto x = static_cast<float>(uniform(rng));
    const float back = corsika_m_to_cm(corsika_cm_to_m(x));
    if (back == x) {
      ++exact;
    } else if (std::abs(back - x) <= std::nextafter(std::abs(x), 1.0e30F) - std::abs(x)) {
      ++off_by_one;
    } else {
      ++worse;
    }
  }
  CHECK(worse == 0);
  std::printf("float round trip cm->m->cm: %llu exact, %llu off by 1 ulp of %llu\n",
              static_cast<unsigned long long>(exact), static_cast<unsigned long long>(off_by_one),
              static_cast<unsigned long long>(n));
}

void stream_a_records() {
  const std::vector<Record> records = convert(build_stream_a());
  using T = RecordType;
  const std::vector<RecordType> expected{
      T::run,           T::blob,      T::blob,      T::blob,          T::event,     T::blob,
      T::blob,          T::telescope, T::photons,   T::telescope_end, T::telescope, T::photons,
      T::telescope_end, T::blob,      T::event_end, T::blob,          T::run_end};
  CHECK(types_of(records) == expected);
  if (types_of(records) != expected) {
    return;
  }

  const Metadata &run = records[0].metadata;
  CHECK(meta<std::int64_t>(run, "run_number") == 12345);
  CHECK(meta<double>(run, "observation_level_m") == 22000.0);
  CHECK(meta<double>(run, "spectral_index") == -2.0);
  CHECK(meta<double>(run, "energy_min_gev") == 10.0);
  CHECK(meta<double>(run, "energy_max_gev") == 100000.0);
  CHECK(meta<std::int64_t>(run, "n_telescopes") == 3);
  CHECK(meta<std::vector<double>>(run, "telescope_x_m") == (std::vector<double>{0.0, 12.0, -5.0}));
  CHECK(meta<std::vector<double>>(run, "telescope_z_m") == (std::vector<double>{18.0, 19.0, 20.0}));
  CHECK(meta<std::int64_t>(run, "atmosphere_id") == 5);
  CHECK(meta<std::int64_t>(run, "refraction_applied") == 1);
  CHECK(meta<std::int64_t>(run, "curvature_applied") == 0);
  CHECK(meta<double>(run, "wavelength_min_nm") == 250.0);
  CHECK(meta<double>(run, "wavelength_max_nm") == 800.0);

  const Metadata &event = records[4].metadata;
  CHECK(meta<std::int64_t>(event, "event_number") == 1);
  CHECK(meta<std::int64_t>(event, "primary_id") == 14);
  CHECK(meta<double>(event, "energy_gev") == 1000.0);
  CHECK(meta<double>(event, "zenith_rad") == static_cast<double>(0.3F));
  CHECK(meta<double>(event, "azimuth_rad") == static_cast<double>(1.2F));
  CHECK(meta<double>(event, "first_interaction_m") == 25000.0);
  CHECK(meta<std::vector<double>>(event, "core_x_m") == (std::vector<double>{0.0, 15.0}));
  CHECK(meta<std::vector<double>>(event, "core_y_m") == (std::vector<double>{-3.0, 0.0}));
  CHECK(meta<double>(event, "time_offset_ns") == 12.5);

  const Metadata &tel = records[7].metadata;
  CHECK(meta<std::int64_t>(tel, "telescope_id") == 1);
  CHECK(meta<std::int64_t>(tel, "array_index") == 0);
  CHECK(meta<std::int64_t>(tel, "frame") == static_cast<std::int64_t>(Frame::ground));
  CHECK(meta<double>(tel, "x_m") == 0.0);
  CHECK(meta<double>(tel, "z_m") == 18.0);

  const PhotonChunk &chunk = records[8].chunk;
  CHECK(chunk.telescope_id == 1);
  CHECK(chunk.frame == Frame::ground);
  CHECK(chunk.first_index == 0);
  CHECK(chunk.rows() == 2);
  CHECK(chunk[Column::x] == (std::vector<float>{1.0F, 12.0F}));
  CHECK(chunk[Column::y] == (std::vector<float>{-2.0F, 3.0F}));
  CHECK(chunk[Column::cx] == (std::vector<float>{0.1F, -0.1F}));
  CHECK(chunk[Column::t] == (std::vector<float>{35.5F, 36.25F}));
  CHECK(chunk[Column::zem] == (std::vector<float>{22000.0F, 22000.0F}));
  CHECK(chunk[Column::w] == (std::vector<float>{1.5F, 1.5F}));
  CHECK(chunk[Column::wavelength] == (std::vector<float>{0.0F, 0.0F}));
  CHECK(!chunk.has(Column::z));

  CHECK(meta<double>(records[9].metadata, "photons") == 3.0);
  CHECK(records[11].chunk[Column::x] == (std::vector<float>{-3.0F}));
  CHECK(records[11].chunk[Column::wavelength] == (std::vector<float>{350.0F}));
  CHECK(records[14].type == RecordType::event_end);
  CHECK(records[16].type == RecordType::run_end);

  const Record &rune = records[15];
  CHECK(rune.blob_eventio_type == iact::rune && rune.blob_eventio_id == 12345);
  const Record &telpos = records[3];
  CHECK(telpos.blob_eventio_type == iact::telpos);
  CHECK(telpos.blob == telpos_block().raw());
}

void curvature_flag_is_read_from_evth() {
  const std::vector<Record> standard = convert(build_stream_a(0.0F));
  CHECK(!standard.empty() && standard.front().type == RecordType::run);
  if (standard.empty() || standard.front().type != RecordType::run) {
    return;
  }
  CHECK(meta<std::int64_t>(standard.front().metadata, "curvature_applied") == 0);

  for (const float flag : {1.0F, 2.0F}) {
    const std::vector<Record> curved = convert(build_stream_a(flag));
    CHECK(!curved.empty() && curved.front().type == RecordType::run);
    if (curved.empty() || curved.front().type != RecordType::run) {
      return;
    }
    CHECK(meta<std::int64_t>(curved.front().metadata, "curvature_applied") == 1);
  }
}

void stream_a_round_trips_byte_identically() {
  const std::string input = build_stream_a();
  std::ostringstream out;
  SimtelWriter writer{out};
  std::istringstream in{input};
  corsika_to_records(in, [&](const Record &r) { writer.write(r); });
  writer.finish();
  CHECK(out.str() == input);
}

void stream_b_records_and_split_mode() {
  const std::vector<Record> records = convert(build_stream_b());
  using T = RecordType;
  // Split-mode markers 1213 and 1214 are blobs in input order (records 21, 25).
  const std::vector<RecordType> expected{
      T::run,           T::blob,      T::blob,    T::blob,          T::event,     T::blob,
      T::blob,          T::telescope, T::photons, T::telescope_end, T::telescope, T::photons,
      T::telescope_end, T::telescope, T::photons, T::telescope_end, T::blob,      T::event_end,
      T::event,         T::blob,      T::blob,    T::blob,          T::telescope, T::photons,
      T::telescope_end, T::blob,      T::blob,    T::event_end,     T::blob,      T::run_end};
  CHECK(types_of(records) == expected);
  if (types_of(records) != expected) {
    return;
  }
  CHECK(records[21].blob_eventio_type == iact::telarray_head);
  CHECK(records[25].blob_eventio_type == iact::telarray_end);
  // Telescope 3 is compact: values are quantised, positions in m.
  const PhotonChunk &compact = records[14].chunk;
  CHECK(compact.rows() == 2);
  if (compact.rows() == 2) {
    CHECK(compact.telescope_id == 3);
    CHECK(std::abs(compact[Column::x][0] - 0.123F) < 1.0e-4F);
    CHECK(std::abs(compact[Column::y][0] + 0.045F) < 1.0e-4F);
    CHECK(std::abs(compact[Column::x][1] + 12.0F) < 1.0e-3F);
    CHECK(std::abs(compact[Column::zem][1] - 1000.0F) < 2.4F);
  }
  // Split-mode telescope: array reuse 1, telescope 1 (ident 1000 * 1 + 0).
  CHECK(meta<std::int64_t>(records[22].metadata, "array_index") == 1);
  CHECK(meta<std::int64_t>(records[22].metadata, "telescope_id") == 1);
  CHECK(records[23].chunk[Column::x] == (std::vector<float>{corsika_cm_to_m(5.0F)}));
}

void split_mode_is_written_as_normal_container() {
  const std::string input = build_stream_b();
  std::ostringstream out;
  SimtelWriter writer{out};
  std::istringstream source{input};
  corsika_to_records(source, [&](const Record &r) { writer.write(r); });
  writer.finish();
  std::istringstream in{out.str()};
  EventioReader reader{in};
  EventioBlock block;
  int containers = 0;
  int photon_blocks = 0;
  int markers = 0;
  while (reader.next(block)) {
    if (block.type == iact::telarray) {
      ++containers;
      for (const EventioBlock &child : subblocks(block)) {
        photon_blocks += child.type == iact::photons ? 1 : 0;
      }
    }
    markers += (block.type == iact::telarray_head || block.type == iact::telarray_end) ? 1 : 0;
  }
  CHECK(containers == 2);
  CHECK(photon_blocks == 4);
  CHECK(markers == 0);
}

void malformed_input_throws() {
  const std::string valid = build_stream_a();
  std::string big_endian = valid;
  big_endian[0] = '\xD4';
  big_endian[1] = '\x1F';
  big_endian[2] = '\x8A';
  big_endian[3] = '\x37';
  CHECK(throws_error([&] { convert(big_endian); }));
  CHECK(throws_error([&] { convert(valid.substr(0, valid.size() / 2)); }));

  // Photons outside an event.
  const std::string orphan = serialize({long_block(0, 1.0F, tel1_bunches)});
  CHECK(throws_error([&] { convert(orphan); }));
}

// Block-aligned truncation (no EVTE, no RUNE) must not count as a complete file.
void truncated_input_throws() {
  CHECK(throws_error([&] { convert(serialize({runh_block(), telpos_block(), evth_block(1)})); }));
  CHECK(throws_error(
      [&] { convert(serialize({runh_block(), telpos_block(), evth_block(1), evte_block(1)})); }));
  CHECK(throws_error([&] { convert(serialize({runh_block()})); }));
  // A run with RUNE but no events is complete.
  CHECK(!throws_error([&] { convert(serialize({runh_block(), rune_block()})); }));
}

// A raw 1205 block (no range checks on write) with the given indices.
EventioBlock raw_photons_block(std::int16_t array, std::int16_t tel) {
  Bytes body;
  ByteWriter out{body};
  out.i16(array);
  out.i16(tel);
  out.f32(0.0F);
  out.i32(0);
  return make_block(iact::photons, 0, 0, std::move(body));
}

// Negative indices are rejected; the particle block 999/999 stays a blob.
void negative_photon_indices_are_rejected() {
  CHECK(throws_error([&] {
    convert(serialize(
        {runh_block(), evth_block(1), raw_photons_block(0, -1), evte_block(1), rune_block()}));
  }));
  CHECK(throws_error([&] {
    convert(serialize(
        {runh_block(), evth_block(1), raw_photons_block(-1, 0), evte_block(1), rune_block()}));
  }));
  TelPhotons particle;
  particle.array = 999;
  particle.tel = 999;
  particle.photons = 1.0F;
  const std::vector<Record> records = convert(serialize(
      {runh_block(), evth_block(1), make_photons_block(particle), evte_block(1), rune_block()}));
  bool particle_blob = false;
  for (const Record &r : records) {
    particle_blob |= r.type == RecordType::blob && r.blob_eventio_type == iact::photons;
  }
  CHECK(particle_blob);
  const std::vector<RecordType> types = types_of(records);
  CHECK(std::find(types.begin(), types.end(), RecordType::telescope) == types.end());
}

PhotonChunk photon_chunk(std::uint32_t telescope_id, std::uint64_t first_index, std::size_t rows) {
  PhotonChunk chunk;
  chunk.telescope_id = telescope_id;
  chunk.frame = Frame::ground;
  chunk.first_index = first_index;
  chunk.resize(rows, column_bit(Column::x) | column_bit(Column::y) | column_bit(Column::cx) |
                         column_bit(Column::cy) | column_bit(Column::t) | column_bit(Column::zem) |
                         column_bit(Column::w));
  return chunk;
}

Record telescope_record(std::int64_t telescope_id) {
  Record r;
  r.type = RecordType::telescope;
  r.metadata["telescope_id"] = telescope_id;
  r.metadata["array_index"] = std::int64_t{0};
  return r;
}

Record photons_record(PhotonChunk chunk) {
  Record r;
  r.type = RecordType::photons;
  r.chunk = std::move(chunk);
  return r;
}

Record telescope_end_record(std::int64_t telescope_id) {
  Record r;
  r.type = RecordType::telescope_end;
  r.metadata["telescope_id"] = telescope_id;
  r.metadata["photons"] = 0.0;
  return r;
}

void feed_writer(const std::vector<Record> &records) {
  std::ostringstream out;
  SimtelWriter writer{out};
  for (const Record &r : records) {
    writer.write(r);
  }
  writer.finish();
}

// An array reuse is another telescope record of the same telescope and event: its photon
// indices continue, so (event, telescope, index) stays unique.
void chunk_indices_continue_across_array_reuses() {
  Record reuse = telescope_record(1);
  reuse.metadata["array_index"] = std::int64_t{1};
  CHECK(!throws_error([&] {
    feed_writer({telescope_record(1), photons_record(photon_chunk(1, 0, 2)),
                 telescope_end_record(1), reuse, photons_record(photon_chunk(1, 2, 3)),
                 telescope_end_record(1)});
  }));
  CHECK(throws_error([&] {
    feed_writer({telescope_record(1), photons_record(photon_chunk(1, 0, 2)),
                 telescope_end_record(1), reuse, photons_record(photon_chunk(1, 0, 3)),
                 telescope_end_record(1)});
  }));
}

// Chunk indices continue per telescope; gaps, overlaps and reordering throw.
void chunk_indices_must_be_contiguous() {
  const Record tel = telescope_record(1);
  CHECK(!throws_error([&] {
    feed_writer({tel, photons_record(photon_chunk(1, 0, 2)), photons_record(photon_chunk(1, 2, 2)),
                 telescope_end_record(1)});
  }));
  CHECK(throws_error([&] {
    feed_writer({tel, photons_record(photon_chunk(1, 0, 2)), photons_record(photon_chunk(1, 3, 2)),
                 telescope_end_record(1)});
  }));
  CHECK(throws_error([&] {
    feed_writer({tel, photons_record(photon_chunk(1, 0, 2)), photons_record(photon_chunk(1, 1, 2)),
                 telescope_end_record(1)});
  }));
  CHECK(throws_error([&] {
    feed_writer({tel, photons_record(photon_chunk(1, 2, 2)), photons_record(photon_chunk(1, 0, 2)),
                 telescope_end_record(1)});
  }));
  CHECK(!throws_error([&] {
    feed_writer({tel, photons_record(photon_chunk(1, 0, 0)), photons_record(photon_chunk(1, 0, 1)),
                 telescope_end_record(1)});
  }));
  // The next telescope starts again at index 0.
  CHECK(!throws_error([&] {
    feed_writer({tel, photons_record(photon_chunk(1, 0, 2)), telescope_end_record(1),
                 telescope_record(2), photons_record(photon_chunk(2, 0, 3)),
                 telescope_end_record(2)});
  }));
}

// telescope_end must name the open telescope.
void telescope_end_must_match_open_telescope() {
  const Record tel = telescope_record(1);
  CHECK(throws_error(
      [&] { feed_writer({tel, photons_record(photon_chunk(1, 0, 2)), telescope_end_record(2)}); }));
  Record no_id;
  no_id.type = RecordType::telescope_end;
  CHECK(throws_error([&] { feed_writer({tel, photons_record(photon_chunk(1, 0, 2)), no_id}); }));
  CHECK(!throws_error(
      [&] { feed_writer({tel, photons_record(photon_chunk(1, 0, 2)), telescope_end_record(1)}); }));
}

// An open telescope must not be dropped by an array flush or finish().
void open_telescope_blocks_array_flush() {
  const Record tel = telescope_record(1);
  Record evte;
  evte.type = RecordType::blob;
  evte.blob_eventio_type = iact::evte;
  CHECK(throws_error([&] { feed_writer({tel, photons_record(photon_chunk(1, 0, 2)), evte}); }));
  CHECK(throws_error([&] { feed_writer({tel, photons_record(photon_chunk(1, 0, 2))}); }));
  // Array indices above 32767 are rejected before anything is written.
  Record big = telescope_record(1);
  big.metadata["array_index"] = std::int64_t{40000};
  CHECK(throws_error([&] { feed_writer({big}); }));
}

// ---------------------------------------------------------------------------
// 3D photon bunches (1217) and sphere runs.

// Bytes of data/photons3d.iact, written by data/make_photons3d_fixture.c with the CORSIKA IACT
// library (provenance in data/photons3d_provenance.txt).
std::string read_fixture(const char *name) {
  const std::string source = __FILE__; // .../cpp/tests/test_corsika.cpp
  const std::string dir = source.substr(0, source.rfind('/'));
  std::ifstream in{dir + "/data/" + name, std::ios::binary};
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

// A 1217 block as CORSIKA writes it (io_simtel.c write_tel_photons3d): short array, short tel,
// real photons, long n, then n bunches of 10 reals (x y z cx cy cz ctime dist photons lambda).
EventioBlock photons3d_block(int tel, const std::vector<Bunch3d> &bunches, float photons = 0.0F,
                             int array = 0) {
  Bytes body;
  ByteWriter out{body};
  out.i16(static_cast<std::int16_t>(array));
  out.i16(static_cast<std::int16_t>(tel));
  out.f32(photons);
  out.i32(static_cast<std::int32_t>(bunches.size()));
  for (const Bunch3d &b : bunches) {
    for (const float v : {b.x, b.y, b.z, b.cx, b.cy, b.cz, b.ctime, b.dist, b.photons, b.lambda}) {
      out.f32(v);
    }
  }
  return make_block(iact::photons3d, 0, static_cast<std::uint32_t>(1000 * array + tel),
                    std::move(body));
}

Bunch3d bunch3(float x, float y, float z, float cx, float cy, float cz, float ctime, float dist,
               float photons, float lambda = 0.0F) {
  Bunch3d b;
  b.x = x;
  b.y = y;
  b.z = z;
  b.cx = cx;
  b.cy = cy;
  b.cz = cz;
  b.ctime = ctime;
  b.dist = dist;
  b.photons = photons;
  b.lambda = lambda;
  return b;
}

// A run of the synthetic file of telpos_block() (telescopes of radius 12 m, 12 m and 12 m) with
// one event holding the given container children.
std::string sphere_stream(const std::vector<EventioBlock> &children) {
  return serialize({runh_block(), telpos_block(), evth_block(1),
                    make_container(iact::telarray, 0, 0, children), evte_block(1), rune_block()});
}

// Expected G sphere row of one bunch of the fixture (CORSIKA bunches, see
// make_photons3d_fixture.c, converted by hand: m, ns, m a.s.l.).
struct ExpectedRow {
  double x, y, z, cx, cy, cz, t_ns, path_m, zem_m, w, wavelength;
};

// Telescope 1: radius 12 m, centre at 18 m above the observation level (2200 m).
const std::vector<ExpectedRow> fixture_tel1_rows{
    // downward, arrival on the sphere top; E = (0, 0, 42) m
    {0.0, 0.0, 12.0, 0.0, 0.0, -1.0, 100.0, 30.0, 2260.0, 2.0, 0.0},
    // downward, emitted inside at (3, 0, 4) m; arrival 8 m below it
    {3.0, 0.0, 4.0, 0.0, 0.0, -1.0, 50.0 - 8.0 / speed_of_light_m_per_ns, 0.0, 2222.0, 5.0, 350.0},
    // upward from 30 m below the centre; enters the sphere from below at s = 18 m
    {0.0, 0.0, -12.0, 0.0, 0.0, 1.0, 200.0 - 12.0 / speed_of_light_m_per_ns, 18.0, 2188.0, 3.0,
     0.0},
    // horizontal from (-20, 0, 0) m; enters at x = -12 m after 8 m
    {-12.0, 0.0, 0.0, 1.0, 0.0, 0.0, 300.0 - 12.0 / speed_of_light_m_per_ns, 8.0, 2218.0, 4.0, 0.0},
    // arrival (9, 0, -12) m outside; E = (-9, 0, 12) m, the line enters at s = 3 m
    {-7.2, 0.0, 9.6, 0.6, 0.0, -0.8, 400.0 - 27.0 / speed_of_light_m_per_ns, 3.0, 2230.0, 6.0, 0.0},
    // downward, arrival on the sphere top; E = (0, 0, 30) m (arrival 18 m after emission)
    {0.0, 0.0, 12.0, 0.0, 0.0, -1.0, 700.0, 18.0, 2248.0, 1.5, 532.0},
};

// Telescope 2: radius 8 m, centre at 20 m above the observation level.
const std::vector<ExpectedRow> fixture_tel2_rows{
    // downward, arrival on the sphere top; E = (0, 0, 28) m
    {0.0, 0.0, 8.0, 0.0, 0.0, -1.0, 120.0, 20.0, 2248.0, 1.0, 0.0},
    // upward (cz = 0.8), E = (-9.4, 0, -8.8) m outside; enters at s = 5 m
    {-6.4, 0.0, -4.8, 0.6, 0.0, 0.8, 150.0 - 5.0 / speed_of_light_m_per_ns, 5.0, 2211.2, 4.0,
     532.0},
    // downward, emitted inside at (0, 0, 3) m
    {0.0, 0.0, 3.0, 0.0, 0.0, -1.0, 60.0 - 5.0 / speed_of_light_m_per_ns, 0.0, 2223.0, 2.0, 0.0},
};

void check_row(const PhotonChunk &chunk, std::size_t i, const ExpectedRow &e) {
  CHECK_NEAR(chunk[Column::x][i], e.x, 1.0e-5);
  CHECK_NEAR(chunk[Column::y][i], e.y, 1.0e-5);
  CHECK_NEAR(chunk[Column::z][i], e.z, 1.0e-5);
  CHECK_NEAR(chunk[Column::cx][i], e.cx, 1.0e-6);
  CHECK_NEAR(chunk[Column::cy][i], e.cy, 1.0e-6);
  CHECK_NEAR(chunk[Column::cz][i], e.cz, 1.0e-6);
  CHECK_NEAR(chunk[Column::t][i], e.t_ns, 1.0e-3);
  CHECK_NEAR(chunk[Column::zem][i], e.zem_m, 1.0e-3);
  CHECK_NEAR(chunk[Column::w][i], e.w, 1.0e-6);
  CHECK_NEAR(chunk[Column::path][i], e.path_m, 1.0e-5);
  CHECK_NEAR(chunk[Column::wavelength][i], e.wavelength, 1.0e-3);
}

// Every crafted bunch of the fixture: row, path, time and zem against hand-derived values; the two
// misses are counted in the run end and not written.
void photons3d_fixture_sphere_rows() {
  const std::vector<Record> records = convert(read_fixture("photons3d.iact"));
  CHECK(!records.empty());
  if (records.empty()) {
    return;
  }
  using T = RecordType;
  // RUNH, TELPOS before the event; EVTH and TELOFF inside it; EVTE and RUNE after the telescopes.
  const std::vector<RecordType> expected{
      T::run,       T::blob,      T::blob,          T::event,     T::blob,    T::blob,
      T::telescope, T::photons,   T::telescope_end, T::telescope, T::photons, T::telescope_end,
      T::blob,      T::event_end, T::blob,          T::run_end};
  CHECK(types_of(records) == expected);
  if (types_of(records) != expected) {
    return;
  }

  const Metadata &run = records[0].metadata;
  CHECK(meta<std::string>(run, "ground_representation") == "sphere");
  CHECK(meta<std::int64_t>(run, "run_number") == 1);
  CHECK(meta<double>(run, "observation_level_m") == 2200.0);
  CHECK(meta<std::int64_t>(run, "curvature_applied") == 0);
  CHECK(meta<std::int64_t>(records[15].metadata, "misses_sphere") == 2);

  CHECK(meta<std::int64_t>(records[6].metadata, "telescope_id") == 1);
  CHECK(meta<double>(records[6].metadata, "r_m") == 12.0);
  CHECK(meta<double>(records[6].metadata, "z_m") == 18.0);
  CHECK(meta<double>(records[8].metadata, "photons") == 37.5);
  CHECK(meta<std::int64_t>(records[9].metadata, "telescope_id") == 2);
  CHECK(meta<double>(records[9].metadata, "r_m") == 8.0);
  CHECK(meta<double>(records[9].metadata, "z_m") == 20.0);
  CHECK(meta<double>(records[11].metadata, "photons") == 7.0);

  const PhotonChunk &tel1 = records[7].chunk;
  const PhotonChunk &tel2 = records[10].chunk;
  CHECK(tel1.telescope_id == 1);
  CHECK(tel2.telescope_id == 2);
  CHECK(tel1.first_index == 0 && tel2.first_index == 0);
  CHECK(tel1.rows() == fixture_tel1_rows.size());
  CHECK(tel2.rows() == fixture_tel2_rows.size());
  if (tel1.rows() != fixture_tel1_rows.size() || tel2.rows() != fixture_tel2_rows.size()) {
    return;
  }
  for (std::size_t i = 0; i < tel1.rows(); ++i) {
    check_row(tel1, i, fixture_tel1_rows[i]);
  }
  for (std::size_t i = 0; i < tel2.rows(); ++i) {
    check_row(tel2, i, fixture_tel2_rows[i]);
  }
  const std::uint32_t sphere_mask =
      column_bit(Column::x) | column_bit(Column::y) | column_bit(Column::z) |
      column_bit(Column::cx) | column_bit(Column::cy) | column_bit(Column::cz) |
      column_bit(Column::t) | column_bit(Column::zem) | column_bit(Column::w) |
      column_bit(Column::path) | column_bit(Column::wavelength);
  CHECK(tel1.mask() == sphere_mask && tel2.mask() == sphere_mask);
  CHECK(!throws_error([&] { validate_sphere_rows(tel1, 12.0); }));
  CHECK(!throws_error([&] { validate_sphere_rows(tel2, 8.0); }));
}

// Two conversions of the same bytes give the same records (no hidden state across runs).
void photons3d_conversion_is_deterministic() {
  const std::string bytes = read_fixture("photons3d.iact");
  std::ostringstream first;
  std::ostringstream second;
  std::istringstream in_first{bytes};
  std::istringstream in_second{bytes};
  corsika_to_stream(in_first, first, "test");
  corsika_to_stream(in_second, second, "test");
  CHECK(!first.str().empty());
  CHECK(first.str() == second.str());
}

// The converted sphere run is a valid format 1.2 stream: it reads back and its rows pass
// validate_sphere_rows against the radius of their telescope.
void photons3d_stream_reads_back_as_sphere_run() {
  std::istringstream in{read_fixture("photons3d.iact")};
  std::ostringstream out;
  corsika_to_stream(in, out, "test");
  std::istringstream back{out.str()};
  Reader reader{back};
  Record record;
  std::vector<Record> records;
  while (reader.next(record)) {
    records.push_back(record);
  }
  CHECK(reader.ended());
  CHECK(!records.empty());
  std::size_t photon_chunks = 0;
  for (const Record &r : records) {
    if (r.type == RecordType::run) {
      CHECK(ground_representation(r.metadata) == GroundRepresentation::sphere);
    }
    if (r.type == RecordType::photons) {
      ++photon_chunks;
      const double radius = r.chunk.telescope_id == 1 ? 12.0 : 8.0;
      CHECK(!throws_error([&] { validate_sphere_rows(r.chunk, radius); }));
    }
  }
  CHECK(photon_chunks == 2);
}

// Pure geometry, analytic cases (radius 12 m, centre at the origin; metres and ns).
void sphere_bunch_geometry() {
  const double radius = 12.0;
  const double altitude = 2200.0;
  // Arrival (9, 0, -12) m, d = (0.6, 0, -0.8), dist 30 m: E = (-9, 0, 12), |E| = 15 > 12, the
  // entry is at s = 3 m (row (-7.2, 0, 9.6)).
  {
    const detail::SphereBunch s = detail::sphere_bunch(
        bunch3(900.0F, 0.0F, -1200.0F, 0.6F, 0.0F, -0.8F, 400.0F, 3000.0F, 1.0F), radius, altitude);
    CHECK(s.hits);
    CHECK_NEAR(s.x, -7.2, 1.0e-5);
    CHECK_NEAR(s.z, 9.6, 1.0e-5);
    CHECK_NEAR(s.path, 3.0, 1.0e-5);
    CHECK_NEAR(s.time_ns, 400.0 - 27.0 / speed_of_light_m_per_ns, 1.0e-5);
    CHECK_NEAR(s.zem_m, altitude + 12.0, 1.0e-5);
  }
  // Tangent: E = (-20, 12, 0) m, d = (1, 0, 0): one root, s = 20 m, row (0, 12, 0).
  {
    const detail::SphereBunch s = detail::sphere_bunch(
        bunch3(2000.0F, 1200.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 4000.0F, 1.0F), radius, altitude);
    CHECK(s.hits);
    CHECK_NEAR(s.x, 0.0, 1.0e-5);
    CHECK_NEAR(s.y, 12.0, 1.0e-5);
    CHECK_NEAR(s.path, 20.0, 1.0e-5);
  }
  // Passes at y = 13 m: no root.
  CHECK(!detail::sphere_bunch(bunch3(2000.0F, 1300.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 4000.0F, 1.0F),
                              radius, altitude)
             .hits);
  // Emitted at (20, 0, 0) m and pointing away (+x): no root with s >= 0.
  CHECK(!detail::sphere_bunch(bunch3(3000.0F, 0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 1000.0F, 1.0F),
                              radius, altitude)
             .hits);
  // Emitted at the sphere surface: counts as inside, path 0.
  {
    const detail::SphereBunch s = detail::sphere_bunch(
        bunch3(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, -1.0F, 10.0F, 1200.0F, 1.0F), radius, altitude);
    CHECK(s.hits);
    CHECK_NEAR(s.z, 12.0, 1.0e-5);
    CHECK_NEAR(s.path, 0.0, 1.0e-5);
  }
  // Invalid bunches are errors, never misses.
  CHECK(throws_error([&] {
    detail::sphere_bunch(bunch3(0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -1.0F, 0.0F, 3000.0F, -1.0F),
                         radius, altitude);
  }));
  CHECK(throws_error([&] {
    detail::sphere_bunch(bunch3(0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -0.5F, 0.0F, 3000.0F, 1.0F),
                         radius, altitude);
  }));
  CHECK(throws_error([&] {
    detail::sphere_bunch(bunch3(0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -1.0F, 0.0F, -3000.0F, 1.0F),
                         radius, altitude);
  }));
  CHECK(throws_error([&] {
    detail::sphere_bunch(bunch3(0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -1.0F, NAN, 3000.0F, 1.0F), radius,
                         altitude);
  }));
}

// Bunches that miss the sphere are counted in misses_sphere and not written; the rows of the other
// bunches keep their order and the photon indices count written rows only.
void photons3d_misses_are_counted() {
  const std::vector<Bunch3d> bunches{
      bunch3(0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -1.0F, 10.0F, 3000.0F, 1.0F),   // hits
      bunch3(-3000.0F, 1500.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 3000.0F, 9.0F), // misses
      bunch3(0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -1.0F, 20.0F, 3000.0F, 2.0F),   // hits
  };
  const std::vector<Record> records = convert(sphere_stream({photons3d_block(0, bunches, 12.0F)}));
  CHECK(!records.empty() && records.front().type == RecordType::run);
  if (records.empty()) {
    return;
  }
  CHECK(meta<std::string>(records.front().metadata, "ground_representation") == "sphere");
  CHECK(meta<std::int64_t>(records.back().metadata, "misses_sphere") == 1);
  const PhotonChunk *chunk = nullptr;
  for (const Record &r : records) {
    if (r.type == RecordType::photons) {
      chunk = &r.chunk;
    }
  }
  CHECK(chunk != nullptr && chunk != nullptr && chunk->rows() == 2);
  if (chunk != nullptr && chunk->rows() == 2) {
    CHECK(chunk->first_index == 0);
    CHECK((*chunk)[Column::w] == (std::vector<float>{1.0F, 2.0F}));
    CHECK_NEAR((*chunk)[Column::t][1], 20.0, 1.0e-3); // path = dist: arrival time kept
  }
}

// One representation per run: 1205 and 1217 blocks may not mix, in an array or across events.
void photon_representations_do_not_mix() {
  const std::vector<Bunch3d> hits{
      bunch3(0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -1.0F, 0.0F, 3000.0F, 1.0F)};
  // 1205 and 1217 in one array container.
  CHECK(throws_error([&] {
    convert(sphere_stream({long_block(0, 1.0F, tel0_bunches), photons3d_block(1, hits)}));
  }));
  CHECK(throws_error([&] {
    convert(sphere_stream({photons3d_block(0, hits), long_block(1, 1.0F, tel1_bunches)}));
  }));
  // 1205 in event 1, 1217 in event 2 (and the other way round).
  const std::string plane_then_sphere = serialize(
      {runh_block(), telpos_block(), evth_block(1),
       make_container(iact::telarray, 0, 0, {long_block(0, 1.0F, tel0_bunches)}), evte_block(1),
       evth_block(2), make_container(iact::telarray, 0, 0, {photons3d_block(0, hits)}),
       evte_block(2), rune_block()});
  CHECK(throws_error([&] { convert(plane_then_sphere); }));
  const std::string sphere_then_plane = serialize(
      {runh_block(), telpos_block(), evth_block(1),
       make_container(iact::telarray, 0, 0, {photons3d_block(0, hits)}), evte_block(1),
       evth_block(2), make_container(iact::telarray, 0, 0, {long_block(0, 1.0F, tel0_bunches)}),
       evte_block(2), rune_block()});
  CHECK(throws_error([&] { convert(sphere_then_plane); }));
  // Two sphere events are fine; the run key is written once, in the run record.
  const std::vector<Record> records = convert(
      serialize({runh_block(), telpos_block(), evth_block(1),
                 make_container(iact::telarray, 0, 0, {photons3d_block(0, hits)}), evte_block(1),
                 evth_block(2), make_container(iact::telarray, 0, 0, {photons3d_block(1, hits)}),
                 evte_block(2), rune_block()}));
  CHECK(!records.empty() && records.front().type == RecordType::run);
  if (!records.empty() && records.front().type == RecordType::run) {
    CHECK(meta<std::string>(records.front().metadata, "ground_representation") == "sphere");
  }
}

// 1205 runs are plane runs: no run key (absent means plane, output as before), no sphere counts.
void plane_runs_are_marked_plane() {
  const std::vector<Record> records = convert(build_stream_a());
  CHECK(!records.empty() && records.front().type == RecordType::run);
  if (records.empty() || records.front().type != RecordType::run) {
    return;
  }
  CHECK(!records.front().metadata.contains("ground_representation"));
  CHECK(ground_representation(records.front().metadata) == GroundRepresentation::plane);
  CHECK(records.back().type == RecordType::run_end && records.back().metadata.empty());
}

// The 1217 representation needs RUNH (observation level) and a TELPOS entry for the telescope.
void photons3d_needs_runh_and_telpos() {
  const std::vector<Bunch3d> hits{
      bunch3(0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -1.0F, 0.0F, 3000.0F, 1.0F)};
  CHECK(throws_error([&] {
    convert(serialize({telpos_block(), evth_block(1),
                       make_container(iact::telarray, 0, 0, {photons3d_block(0, hits)}),
                       evte_block(1), rune_block()}));
  }));
  CHECK(throws_error([&] {
    convert(serialize({runh_block(), evth_block(1),
                       make_container(iact::telarray, 0, 0, {photons3d_block(0, hits)}),
                       evte_block(1), rune_block()}));
  }));
  CHECK(throws_error([&] { convert(sphere_stream({photons3d_block(3, hits)})); }));
}

// A 1217 block with the particle indices (999, 999) is a blob, as for 1205.
void photons3d_particle_block_is_a_blob() {
  const std::vector<Bunch3d> bunches{
      bunch3(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, -1.0F, 0.0F, 0.0F, 1.0F, 9999.0F)};
  const std::vector<Record> records =
      convert(serialize({runh_block(), telpos_block(), evth_block(1),
                         photons3d_block(999, bunches, 1.0F, 999), evte_block(1), rune_block()}));
  bool blob = false;
  for (const Record &r : records) {
    blob |= r.type == RecordType::blob && r.blob_eventio_type == iact::photons3d;
    CHECK(r.type != RecordType::telescope);
  }
  CHECK(blob);
}

// Photons written by the converter validate as sphere rows also when the run has no photons.
void sphere_run_without_photons() {
  const std::vector<Record> records = convert(sphere_stream({photons3d_block(0, {}, 0.0F)}));
  CHECK(!records.empty() && records.front().type == RecordType::run);
  if (!records.empty() && records.front().type == RecordType::run) {
    CHECK(meta<std::string>(records.front().metadata, "ground_representation") == "sphere");
    CHECK(meta<std::int64_t>(records.back().metadata, "misses_sphere") == 0);
  }
}

} // namespace

int main() {
  return run_checks([] {
    units_are_double_then_float();
    float_round_trip_is_within_one_ulp();
    stream_a_records();
    curvature_flag_is_read_from_evth();
    stream_a_round_trips_byte_identically();
    stream_b_records_and_split_mode();
    split_mode_is_written_as_normal_container();
    malformed_input_throws();
    truncated_input_throws();
    negative_photon_indices_are_rejected();
    chunk_indices_must_be_contiguous();
    chunk_indices_continue_across_array_reuses();
    telescope_end_must_match_open_telescope();
    open_telescope_blocks_array_flush();
    photons3d_fixture_sphere_rows();
    photons3d_conversion_is_deterministic();
    photons3d_stream_reads_back_as_sphere_run();
    sphere_bunch_geometry();
    photons3d_misses_are_counted();
    photon_representations_do_not_mix();
    plane_runs_are_marked_plane();
    photons3d_needs_runh_and_telpos();
    photons3d_particle_block_is_a_blob();
    sphere_run_without_photons();
  });
}
