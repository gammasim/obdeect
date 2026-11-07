#include "check.hpp"

#include "gstrickd/corsika.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace gstrickd;

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

EventioBlock evth_block(float event) {
  std::vector<float> w(273, 0.0F);
  w[1] = event;
  w[2] = 14.0F;      // gamma
  w[3] = 1000.0F;    // energy [GeV]
  w[6] = 2500000.0F; // first interaction at 25 km [cm]
  w[10] = 0.3F;      // zenith [rad]
  w[11] = 1.2F;      // azimuth [rad]
  w[76] = 5136.0F;   // atmosphere 5 (bits 11-20), refraction (bit 5)
  w[95] = 250.0F;    // Cherenkov wavelength range [nm]
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

std::string build_stream_a() {
  return serialize(
      {runh_block(), make_block(iact::inputcfg, 0, 0, Bytes{'R', 'U', 'N', '\n'}), telpos_block(),
       evth_block(1), teloff_block(),
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
  const std::vector<RecordType> expected{
      T::run,           T::blob,      T::blob,    T::blob,          T::event,     T::blob,
      T::blob,          T::telescope, T::photons, T::telescope_end, T::telescope, T::photons,
      T::telescope_end, T::telescope, T::photons, T::telescope_end, T::blob,      T::event_end,
      T::event,         T::blob,      T::blob,    T::telescope,     T::photons,   T::telescope_end,
      T::blob,          T::event_end, T::blob,    T::run_end};
  CHECK(types_of(records) == expected);
  if (types_of(records) != expected) {
    return;
  }
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
  CHECK(meta<std::int64_t>(records[21].metadata, "array_index") == 1);
  CHECK(meta<std::int64_t>(records[21].metadata, "telescope_id") == 1);
  CHECK(records[22].chunk[Column::x] == (std::vector<float>{corsika_cm_to_m(5.0F)}));
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

} // namespace

int main() {
  return run_checks([] {
    units_are_double_then_float();
    float_round_trip_is_within_one_ulp();
    stream_a_records();
    stream_a_round_trips_byte_identically();
    stream_b_records_and_split_mode();
    split_mode_is_written_as_normal_container();
    malformed_input_throws();
  });
}
