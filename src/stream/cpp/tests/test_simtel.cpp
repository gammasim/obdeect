#include "check.hpp"

#include "obdeect/stream/simtel.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace obdeect::stream;

namespace {

template <class Action> bool throws_error(Action action) {
  try {
    action();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

// Corsika-style float block with one value per CORSIKA word (word n at index n-1).
std::vector<float> words273() { return std::vector<float>(273, 0.0F); }

EventioBlock runh_block() {
  std::vector<float> w = words273();
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
  std::vector<float> w = words273();
  w[1] = event;
  w[2] = 14.0F;      // gamma
  w[3] = 1000.0F;    // energy [GeV]
  w[6] = 2500000.0F; // first interaction at 25 km [cm]
  w[10] = 0.3F;      // zenith [rad]
  w[11] = 1.2F;      // azimuth [rad]
  w[76] = 5136.0F;   // atmosphere 5, refraction on
  w[95] = 250.0F;    // Cherenkov wavelength range [nm]
  w[96] = 800.0F;
  return make_float_block(iact::evth, static_cast<std::uint32_t>(event), w);
}

EventioBlock evte_block(float event) {
  std::vector<float> w = words273();
  w[1] = event;
  return make_float_block(iact::evte, static_cast<std::uint32_t>(event), w);
}

EventioBlock rune_block() {
  std::vector<float> w = words273();
  w[1] = 12345.0F;
  w[2] = 1.0F;
  return make_float_block(iact::rune, 12345, w);
}

EventioBlock teloff_block() {
  TelOffset off;
  off.xoff = {0.0};
  off.yoff = {0.0};
  return make_teloff_block(off);
}

// Telescope positions in cm: x, y, z, r per telescope.
EventioBlock telpos_block(const std::vector<double> &z_cm) {
  TelPos pos;
  for (std::size_t i = 0; i < z_cm.size(); ++i) {
    pos.x.push_back(100.0 * static_cast<double>(i));
    pos.y.push_back(0.0);
    pos.z.push_back(z_cm[i]);
    pos.r.push_back(1200.0);
  }
  return make_telpos_block(pos);
}

// Long-format 1205 block with integer-cm positions, as CORSIKA writes them.
EventioBlock photons_block(int tel, const std::vector<float> &x_cm,
                           const std::vector<float> &y_cm) {
  TelPhotons p;
  p.array = 0;
  p.tel = tel;
  for (std::size_t i = 0; i < x_cm.size(); ++i) {
    Bunch b;
    b.x = x_cm[i];
    b.y = y_cm[i];
    b.cx = 0.125F * static_cast<float>(i % 8) - 0.25F;
    b.cy = 0.0625F * static_cast<float>(i % 5) - 0.125F;
    b.ctime = 0.5F * static_cast<float>(i);
    b.zem = 2000000.0F; // 20 km
    b.photons = (i % 3 == 0) ? 2.0F : 1.0F;
    b.lambda = 400.0F;
    p.photons += b.photons;
    p.bunches.push_back(b);
  }
  return make_photons_block(p);
}

std::vector<EventioBlock> read_blocks(const std::string &bytes) {
  std::istringstream in(bytes);
  EventioReader reader{in};
  std::vector<EventioBlock> blocks;
  EventioBlock block;
  while (reader.next(block)) {
    blocks.push_back(block);
  }
  return blocks;
}

// True when cm -> m (float) -> cm (float) returns the same float.
bool round_trips(float cm) { return corsika_m_to_cm(corsika_cm_to_m(cm)) == cm; }

// Two telescopes of one event in normal mode (one 1204 container), plus RUNH, TELPOS,
// EVTH, TELOFF and EVTE, as in CORSIKA output. Positions are integer cm that round-trip.
std::string corsika_file(int events) {
  std::vector<float> exact;
  for (int k = -60000; exact.size() < 400 && k <= 60000; k += 7) {
    if (round_trips(static_cast<float>(k))) {
      exact.push_back(static_cast<float>(k));
    }
  }
  std::string bytes;
  const auto add = [&bytes](const EventioBlock &block) {
    const std::vector<std::uint8_t> raw = simtel_detail::block_bytes(block);
    bytes.append(raw.begin(), raw.end());
  };
  add(runh_block());
  add(telpos_block({2000.0, 2100.0}));
  for (int e = 1; e <= events; ++e) {
    add(evth_block(static_cast<float>(e)));
    add(teloff_block());
    std::vector<EventioBlock> children;
    for (int tel = 0; tel < 2; ++tel) {
      std::vector<float> x;
      std::vector<float> y;
      for (std::size_t i = 0; i < 200; ++i) {
        x.push_back(exact[(i * 3 + static_cast<std::size_t>(tel) * 11) % exact.size()]);
        y.push_back(exact[(i * 5 + 1) % exact.size()]);
      }
      children.push_back(photons_block(tel, x, y));
    }
    add(make_container(iact::telarray, 0, 0, children));
    add(evte_block(static_cast<float>(e)));
  }
  add(rune_block());
  return bytes;
}

// Photon-stream records of one telescope in one frame.
PhotonChunk chunk_of(std::uint32_t telescope, Frame frame, const std::vector<double> &x,
                     const std::vector<double> &y, const std::vector<double> &cx,
                     const std::vector<double> &cy) {
  PhotonChunk chunk;
  chunk.telescope_id = telescope;
  chunk.frame = frame;
  const std::size_t n = x.size();
  chunk.resize(n, required_columns | column_bit(Column::zem) | column_bit(Column::wavelength));
  for (std::size_t i = 0; i < n; ++i) {
    chunk[Column::x][i] = static_cast<float>(x[i]);
    chunk[Column::y][i] = static_cast<float>(y[i]);
    chunk[Column::cx][i] = static_cast<float>(cx[i]);
    chunk[Column::cy][i] = static_cast<float>(cy[i]);
    chunk[Column::t][i] = 0.0F;
    chunk[Column::w][i] = 1.0F;
    chunk[Column::zem][i] = 20000.0F; // m
    chunk[Column::wavelength][i] = 400.0F;
  }
  return chunk;
}

Record record_of(RecordType type, Metadata metadata = {}) {
  Record record;
  record.type = type;
  record.metadata = std::move(metadata);
  return record;
}

Record photons_of(const PhotonChunk &chunk) {
  Record record;
  record.type = RecordType::photons;
  record.chunk = chunk;
  return record;
}

Record blob_of(const EventioBlock &block) {
  Record record;
  record.type = RecordType::blob;
  record.blob_eventio_type = block.type;
  record.blob_eventio_id = block.ident;
  record.blob = simtel_detail::block_bytes(block);
  return record;
}

// Run metadata of a stream without CORSIKA blobs: two telescopes at (0, 0) and (100, 0) m.
Metadata light_source_run() {
  Metadata m;
  m["run_number"] = std::int64_t{7};
  m["observation_level_m"] = 2200.0;
  m["telescope_x_m"] = std::vector<double>{0.0, 100.0};
  m["telescope_y_m"] = std::vector<double>{0.0, 0.0};
  m["telescope_z_m"] = std::vector<double>{0.0, 0.0};
  m["telescope_r_m"] = std::vector<double>{12.0, 12.0};
  return m;
}

// Writes records through a PsToSimtel and returns the EventIO bytes.
std::string run_adapter(const std::vector<Record> &records, const CameraMap &camera,
                        const std::set<std::string> &allow_pending = {},
                        const std::set<std::string> &scaled = {}) {
  std::ostringstream out;
  SimtelOptions options;
  options.camera = camera;
  options.allow_pending = allow_pending;
  options.scaled = scaled;
  PsToSimtel converter{out, std::move(options)};
  for (const Record &record : records) {
    converter.write(record);
  }
  converter.finish();
  return out.str();
}

void test_camera_map_parsing() {
  std::istringstream good("# comment\n1 1 0\n\n2 2 30.5  # dual\n");
  const CameraMap map = read_camera_map(good);
  CHECK(map.size() == 2);
  CHECK(map.at(2).mirror_class == 2);
  CHECK_NEAR(map.at(2).rotation_deg, 30.5, 1e-12);
  std::istringstream bad_class("1 3 0\n");
  CHECK(throws_error([&] { read_camera_map(bad_class); }));
  std::istringstream duplicate("1 1 0\n1 2 0\n");
  CHECK(throws_error([&] { read_camera_map(duplicate); }));
  std::istringstream extra("1 1 0 0\n");
  CHECK(throws_error([&] { read_camera_map(extra); }));
  std::istringstream nan("1 1 nan\n");
  CHECK(throws_error([&] { read_camera_map(nan); }));
  CHECK(bypass_optics_level(1) == 1);
  CHECK(bypass_optics_level(2) == 2);
}

void test_invalid_camera_map() {
  const CameraMap zero{{0, CameraTelescope{1, 0.0}}};
  const std::vector<Record> synthetic{record_of(RecordType::run, light_source_run()),
                                      record_of(RecordType::run_end)};
  const std::vector<Record> raw{record_of(RecordType::run, light_source_run()),
                                blob_of(telpos_block({2000.0, 2100.0}))};
  CHECK(throws_error([&] { run_adapter(synthetic, zero); }));
  CHECK(throws_error([&] { run_adapter(raw, zero); }));
  for (double rotation :
       {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
    const CameraMap bad{{1, CameraTelescope{2, rotation}}};
    CHECK(throws_error([&] { run_adapter({}, bad); }));
    PhotonChunk chunk = chunk_of(1, Frame::camera, {0.5}, {0.2}, {0.1}, {-0.2});
    CHECK(throws_error([&] { camera_to_bypass_convention(chunk, CameraTelescope{2, rotation}); }));
  }
}

// Hand-computed, input
// (x, y) = (0.5, 0.2) m, (cx, cy) = (0.1, -0.2), camera rotation 30 deg:
//   single mirror: y -> -y, cy -> -cy                 -> (0.5, -0.2), (0.1, 0.2)
//   dual mirror:   x' = x cos r + y sin r = 0.5 cos30 + 0.2 sin30 = 0.5330127019
//                  y' = y cos r - x sin r = 0.2 cos30 - 0.5 sin30 = -0.0767949192
//                  cx' = cx cos r + cy sin r = -0.0133974596
//                  cy' = cy cos r - cx sin r = -0.2232050808
void test_camera_conversion() {
  PhotonChunk single = chunk_of(1, Frame::camera, {0.5}, {0.2}, {0.1}, {-0.2});
  camera_to_bypass_convention(single, CameraTelescope{1, 30.0});
  CHECK_NEAR(single[Column::x][0], 0.5, 1e-6);
  CHECK_NEAR(single[Column::y][0], -0.2, 1e-6);
  CHECK_NEAR(single[Column::cx][0], 0.1, 1e-6);
  CHECK_NEAR(single[Column::cy][0], 0.2, 1e-6);

  PhotonChunk dual = chunk_of(2, Frame::camera, {0.5}, {0.2}, {0.1}, {-0.2});
  camera_to_bypass_convention(dual, CameraTelescope{2, 30.0});
  CHECK_NEAR(dual[Column::x][0], 0.5330127019, 1e-6);
  CHECK_NEAR(dual[Column::y][0], -0.0767949192, 1e-6);
  CHECK_NEAR(dual[Column::cx][0], -0.0133974596, 1e-6);
  CHECK_NEAR(dual[Column::cy][0], -0.2232050808, 1e-6);

  // Dual mirror with no rotation is the identity.
  PhotonChunk identity = chunk_of(2, Frame::camera, {0.5}, {0.2}, {0.1}, {-0.2});
  camera_to_bypass_convention(identity, CameraTelescope{2, 0.0});
  CHECK_NEAR(identity[Column::y][0], 0.2, 1e-6);
  CHECK_NEAR(identity[Column::cy][0], -0.2, 1e-6);
}

// Adapter output for a camera-frame telescope: cm, rotated, labelled ground for sim_telarray.
void test_camera_frame_through_adapter() {
  const CameraMap camera{{2, CameraTelescope{2, 30.0}}};
  const std::vector<Record> records = {
      record_of(RecordType::run, light_source_run()),
      record_of(RecordType::event, Metadata{{"event_number", std::int64_t{1}}}),
      record_of(RecordType::telescope,
                Metadata{{"telescope_id", std::int64_t{2}}, {"array_index", std::int64_t{0}}}),
      photons_of(chunk_of(2, Frame::camera, {0.5}, {0.2}, {0.1}, {-0.2})),
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{2}}}),
      record_of(RecordType::event_end),
      record_of(RecordType::run_end),
  };
  const std::vector<EventioBlock> blocks = read_blocks(run_adapter(records, camera));
  bool found = false;
  for (const EventioBlock &block : blocks) {
    if (block.type != iact::telarray) {
      continue;
    }
    for (const EventioBlock &child : subblocks(block)) {
      const TelPhotons p = decode_photons(child);
      CHECK(p.tel == 1);
      CHECK(p.bunches.size() == 1);
      if (p.bunches.size() == 1) {
        CHECK_NEAR(p.bunches[0].x, 53.30127, 1e-3);
        CHECK_NEAR(p.bunches[0].y, -7.6794919, 1e-3);
        CHECK_NEAR(p.bunches[0].cx, -0.0133974596, 1e-6);
        CHECK_NEAR(p.bunches[0].cy, -0.2232050808, 1e-6);
        CHECK_NEAR(p.bunches[0].zem, 2000000.0, 1.0);
      }
      found = true;
    }
  }
  CHECK(found);
}

// TELPOS z is zeroed for camera-frame telescopes only (see simtel.hpp).
void test_telpos_zeroed_for_camera_only() {
  const CameraMap camera{{1, CameraTelescope{1, 0.0}}};
  const std::vector<Record> records = {
      record_of(RecordType::run, light_source_run()),
      blob_of(telpos_block({2000.0, 2100.0})),
      record_of(RecordType::event, Metadata{{"event_number", std::int64_t{1}}}),
      record_of(RecordType::telescope,
                Metadata{{"telescope_id", std::int64_t{1}}, {"array_index", std::int64_t{0}}}),
      photons_of(chunk_of(1, Frame::camera, {0.5}, {0.2}, {0.1}, {-0.2})),
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{1}}}),
      record_of(RecordType::telescope,
                Metadata{{"telescope_id", std::int64_t{2}}, {"array_index", std::int64_t{0}}}),
      photons_of(chunk_of(2, Frame::ground, {0.5}, {0.2}, {0.1}, {-0.2})),
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{2}}}),
      record_of(RecordType::event_end),
      record_of(RecordType::run_end),
  };
  const std::vector<EventioBlock> blocks = read_blocks(run_adapter(records, camera));
  bool checked = false;
  for (const EventioBlock &block : blocks) {
    if (block.type == iact::telpos) {
      const TelPos pos = decode_telpos(block);
      CHECK(pos.z.size() == 2);
      if (pos.z.size() == 2) {
        CHECK_NEAR(pos.z[0], 0.0, 1e-9);    // camera-frame telescope 1
        CHECK_NEAR(pos.z[1], 2100.0, 1e-3); // ground-frame telescope 2 unchanged
        CHECK_NEAR(pos.r[0], 1200.0, 1e-3); // other columns kept
      }
      checked = true;
      break;
    }
  }
  CHECK(checked);

  // A telescope in the camera map must not carry ground-frame photons.
  const std::vector<Record> mixed = {
      record_of(RecordType::run, light_source_run()),
      record_of(RecordType::event, Metadata{{"event_number", std::int64_t{1}}}),
      record_of(RecordType::telescope,
                Metadata{{"telescope_id", std::int64_t{1}}, {"array_index", std::int64_t{0}}}),
      photons_of(chunk_of(1, Frame::ground, {0.5}, {0.2}, {0.1}, {-0.2})),
  };
  CHECK(throws_error([&] { run_adapter(mixed, camera); }));
  // ... and a camera-frame telescope needs a camera map entry.
  const std::vector<Record> unmapped = {
      record_of(RecordType::run, light_source_run()),
      record_of(RecordType::event, Metadata{{"event_number", std::int64_t{1}}}),
      record_of(RecordType::telescope,
                Metadata{{"telescope_id", std::int64_t{2}}, {"array_index", std::int64_t{0}}}),
      photons_of(chunk_of(2, Frame::camera, {0.5}, {0.2}, {0.1}, {-0.2})),
  };
  CHECK(throws_error([&] { run_adapter(unmapped, camera); }));
}

// Without CORSIKA blobs, RUNH, TELPOS, EVTH, TELOFF, EVTE and RUNE are synthesised.
// Writes records through a PsToSimtel with options and returns the EventIO bytes.
std::string run_adapter_with(const std::vector<Record> &records, SimtelOptions options) {
  std::ostringstream out;
  PsToSimtel converter{out, std::move(options)};
  for (const Record &record : records) {
    converter.write(record);
  }
  converter.finish();
  return out.str();
}

// Telescope-frame chunk: emission points x y z [m], direction cx cy (cz below, or none),
// times t [ns]; w = 1.
PhotonChunk telescope_chunk(std::uint32_t telescope, const std::vector<double> &x,
                            const std::vector<double> &y, const std::vector<double> &z,
                            const std::vector<double> &cx, const std::vector<double> &cy,
                            const std::vector<double> &cz, const std::vector<double> &t) {
  PhotonChunk chunk;
  chunk.telescope_id = telescope;
  chunk.frame = Frame::telescope;
  const std::size_t n = x.size();
  chunk.resize(n, required_columns | column_bit(Column::z) | column_bit(Column::wavelength) |
                      (cz.empty() ? 0U : column_bit(Column::cz)));
  for (std::size_t i = 0; i < n; ++i) {
    chunk[Column::x][i] = static_cast<float>(x[i]);
    chunk[Column::y][i] = static_cast<float>(y[i]);
    chunk[Column::z][i] = static_cast<float>(z[i]);
    chunk[Column::cx][i] = static_cast<float>(cx[i]);
    chunk[Column::cy][i] = static_cast<float>(cy[i]);
    if (!cz.empty()) {
      chunk[Column::cz][i] = static_cast<float>(cz[i]);
    }
    chunk[Column::t][i] = static_cast<float>(t[i]);
    chunk[Column::w][i] = 1.0F;
    chunk[Column::wavelength][i] = 400.0F;
  }
  return chunk;
}

std::vector<Record> telescope_frame_records(std::uint32_t telescope, const PhotonChunk &chunk) {
  return {
      record_of(RecordType::run, light_source_run()),
      record_of(RecordType::event, Metadata{{"event_number", std::int64_t{1}}}),
      record_of(RecordType::telescope, Metadata{{"telescope_id", std::int64_t{telescope}},
                                                {"array_index", std::int64_t{0}}}),
      photons_of(chunk),
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{telescope}}}),
      record_of(RecordType::event_end),
      record_of(RecordType::run_end),
  };
}

// Bunches of the first telescope photon block of an adapter output.
std::vector<Bunch> first_bunches(const std::string &bytes) {
  for (const EventioBlock &block : read_blocks(bytes)) {
    if (block.type != iact::telarray) {
      continue;
    }
    for (const EventioBlock &child : subblocks(block)) {
      return decode_photons(child).bunches;
    }
  }
  return {};
}

// Telescope frame (dual mirror, Bypass_Optics=1), observation level 2200 m. Emission point
// (0.1, 0.2, 0.5) m, t = 3 ns, direction (0.6, 0, cz):
//   cz = -0.8 (towards the plane z = 0): path 0.625 m, crossing x = 0.1 + 0.625 * 0.6 = 0.475 m
//   cz = +0.8 (away from it):            path -0.625 m, crossing x = 0.1 - 0.375 = -0.275 m
//   both: y = 0.2 m, t = 3 + 0.625 / 0.299792458 ns, zem = 2200.5 m (bunch values in cm)
void test_telescope_frame_through_adapter() {
  const CameraMap camera{{1, CameraTelescope{2, 30.0}}};
  const double t_out = 3.0 + 0.625 / 0.299792458;
  const std::vector<Bunch> down = first_bunches(
      run_adapter(telescope_frame_records(
                      1, telescope_chunk(1, {0.1}, {0.2}, {0.5}, {0.6}, {0.0}, {-0.8}, {3.0})),
                  camera));
  CHECK(down.size() == 1);
  if (down.size() == 1) {
    CHECK_NEAR(down[0].x, 47.5, 1e-3);
    CHECK_NEAR(down[0].y, 20.0, 1e-3);
    CHECK_NEAR(down[0].cx, 0.6, 1e-6);
    CHECK_NEAR(down[0].cy, 0.0, 1e-6);
    CHECK_NEAR(down[0].ctime, t_out, 1e-4);
    CHECK_NEAR(down[0].zem, 220050.0, 0.1);
    CHECK_NEAR(down[0].photons, 1.0, 1e-6);
    CHECK_NEAR(down[0].lambda, 400.0, 1e-6);
  }
  const std::vector<Bunch> up = first_bunches(
      run_adapter(telescope_frame_records(
                      1, telescope_chunk(1, {0.1}, {0.2}, {0.5}, {0.6}, {0.0}, {0.8}, {3.0})),
                  camera));
  CHECK(up.size() == 1);
  if (up.size() == 1) {
    CHECK_NEAR(up[0].x, -27.5, 1e-3);
    CHECK_NEAR(up[0].y, 20.0, 1e-3);
    CHECK_NEAR(up[0].cx, 0.6, 1e-6);
    CHECK_NEAR(up[0].ctime, t_out, 1e-4);
    CHECK_NEAR(up[0].zem, 220050.0, 0.1);
  }
  // Without cz the photon moves towards the plane: same bunch as cz = -0.8.
  const std::vector<Bunch> implied = first_bunches(run_adapter(
      telescope_frame_records(1, telescope_chunk(1, {0.1}, {0.2}, {0.5}, {0.6}, {0.0}, {}, {3.0})),
      camera));
  CHECK(implied.size() == 1);
  if (implied.size() == 1 && down.size() == 1) {
    CHECK_NEAR(implied[0].x, down[0].x, 1e-4);
    CHECK_NEAR(implied[0].ctime, down[0].ctime, 1e-4);
  }
}

// Bypass_Optics=1 is reported for a telescope-frame telescope; the camera rotation is unused.
void test_telescope_frame_report_and_telpos() {
  const CameraMap camera{{1, CameraTelescope{2, 30.0}}};
  const std::vector<Record> records = telescope_frame_records(
      1, telescope_chunk(1, {0.1}, {0.2}, {0.5}, {0.6}, {0.0}, {-0.8}, {3.0}));
  std::ostringstream warnings;
  SimtelOptions options;
  options.camera = camera;
  options.warnings = &warnings;
  const std::string bytes = run_adapter_with(records, options);
  CHECK(warnings.str().find("obdeect stream to-simtel: telescope 1: telescope frame, mirror class 2, "
                            "Bypass_Optics=1\n") != std::string::npos);
  for (const EventioBlock &block : read_blocks(bytes)) {
    if (block.type == iact::telpos) {
      const TelPos pos = decode_telpos(block);
      CHECK(pos.z.size() == 2);
      if (pos.z.size() == 2) {
        CHECK_NEAR(pos.z[0], 0.0, 1e-9);
      }
    }
  }
  // Rotation 0 and 30 deg give the same bunches.
  const std::vector<Bunch> other =
      first_bunches(run_adapter(records, {{1, CameraTelescope{2, 0.0}}}));
  const std::vector<Bunch> rotated = first_bunches(bytes);
  CHECK(other.size() == 1 && rotated.size() == 1);
  if (other.size() == 1 && rotated.size() == 1) {
    CHECK_NEAR(other[0].x, rotated[0].x, 1e-6);
    CHECK_NEAR(other[0].cx, rotated[0].cx, 1e-6);
  }
}

// Telescope-frame photons need a dual-mirror camera map entry and a single frame per telescope.
void test_telescope_frame_rejected() {
  const PhotonChunk tframe = telescope_chunk(1, {0.1}, {0.2}, {0.5}, {0.6}, {0.0}, {-0.8}, {3.0});
  const std::vector<Record> records = telescope_frame_records(1, tframe);
  CHECK(throws_error([&] { run_adapter(records, {}); }));                             // no entry
  CHECK(throws_error([&] { run_adapter(records, {{1, CameraTelescope{1, 0.0}}}); })); // single
  CHECK(throws_error([&] { run_adapter(records, {{2, CameraTelescope{2, 0.0}}}); })); // other id
  CHECK(!throws_error([&] { run_adapter(records, {{1, CameraTelescope{2, 0.0}}}); }));

  // A photon parallel to the plane z = 0, and a chunk without z, are errors.
  const CameraMap dual{{1, CameraTelescope{2, 0.0}}};
  CHECK(throws_error([&] {
    run_adapter(telescope_frame_records(
                    1, telescope_chunk(1, {0.1}, {0.2}, {0.5}, {1.0}, {0.0}, {0.0}, {3.0})),
                dual);
  }));
  PhotonChunk no_z = tframe;
  no_z[Column::z].clear();
  CHECK(throws_error([&] { telescope_to_bypass_convention(no_z, 2200.0); }));
  for (const float bad : {std::nanf(""), INFINITY}) {
    PhotonChunk bad_direction = tframe; // with an explicit cz column
    bad_direction[Column::cx][0] = bad;
    CHECK(throws_error([&] { telescope_to_bypass_convention(bad_direction, 2200.0); }));
    PhotonChunk bad_cy = tframe;
    bad_cy[Column::cy][0] = bad;
    CHECK(throws_error([&] { telescope_to_bypass_convention(bad_cy, 2200.0); }));
  }
  // A direction that is not a unit vector, and results out of float range, are errors; the
  // chunk is left unchanged.
  PhotonChunk not_unit = tframe;
  not_unit[Column::cx][0] = 0.0F;
  not_unit[Column::cz][0] = -2.0F;
  CHECK(throws_error([&] { telescope_to_bypass_convention(not_unit, 2200.0); }));
  CHECK_NEAR(not_unit[Column::x][0], 0.1, 1e-6);
  PhotonChunk implied_not_unit = telescope_chunk(1, {0.1}, {0.2}, {0.5}, {0.9}, {0.9}, {}, {3.0});
  try {
    telescope_to_bypass_convention(implied_not_unit, 2200.0);
    CHECK(false);
  } catch (const std::invalid_argument &error) {
    CHECK(std::string(error.what()).find("unit vector") != std::string::npos);
  }
  PhotonChunk tiny_cz = telescope_chunk(1, {0.1}, {0.2}, {5.0}, {1.0}, {0.0}, {-1e-38}, {3.0});
  CHECK(throws_error([&] { telescope_to_bypass_convention(tiny_cz, 2200.0); }));
  CHECK_NEAR(tiny_cz[Column::t][0], 3.0, 1e-6);
  PhotonChunk bad_level = tframe;
  CHECK(throws_error([&] { telescope_to_bypass_convention(bad_level, std::nan("")); }));

  // Camera and telescope frame for one telescope in one run.
  for (const Frame first : {Frame::camera, Frame::telescope}) {
    std::vector<Record> mixed = telescope_frame_records(1, tframe);
    const PhotonChunk camera_chunk = chunk_of(1, Frame::camera, {0.5}, {0.2}, {0.1}, {-0.2});
    const PhotonChunk a = first == Frame::camera ? camera_chunk : tframe;
    const PhotonChunk b = first == Frame::camera ? tframe : camera_chunk;
    mixed[3] = photons_of(a);
    mixed.insert(mixed.begin() + 4, photons_of(b));
    CHECK(throws_error([&] { run_adapter(mixed, dual); }));
  }
  // Ground frame for a camera-map telescope is still an error.
  CHECK(throws_error([&] {
    run_adapter(telescope_frame_records(1, chunk_of(1, Frame::ground, {0.5}, {0.2}, {0.1}, {-0.2})),
                dual);
  }));
}

// Atmospheric table of three levels (a stand-in for a CORSIKA atmprof file).
AtmosphereTable three_level_table() {
  AtmosphereTable table;
  table.profile_id = custom_profile_id;
  table.file_name = "test_atmprof.dat";
  table.alt_km = {0.0, 1.0, 5.0};
  table.rho = {1.2, 1.0, 0.5};
  table.thick = {1000.0, 900.0, 500.0};
  table.refidx_m1 = {2.9e-4, 2.7e-4, 2.2e-4};
  return table;
}

// Light-source records with one telescope: run, one event, end.
std::vector<Record> light_source_records() {
  return {
      record_of(RecordType::run, light_source_run()),
      record_of(RecordType::event, Metadata{{"event_number", std::int64_t{3}}}),
      record_of(RecordType::telescope,
                Metadata{{"telescope_id", std::int64_t{1}}, {"array_index", std::int64_t{0}}}),
      photons_of(chunk_of(1, Frame::ground, {0.5, 1.5}, {0.2, 0.0}, {0.1, 0.0}, {-0.2, 0.0})),
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{1}}}),
      record_of(RecordType::event_end),
      record_of(RecordType::run_end),
  };
}

// Table text, parsing and the rejected variants.
void test_atmprof_text() {
  std::istringstream text("# altitude_km density thickness n-1\n"
                          "0.0  1.2  1000.0 2.9e-4 extra ignored\n"
                          "\n"
                          "  1.0\t1.0\t900.0\t2.7e-4  # trailing comment\n"
                          "5.0 0.5 500.0 2.2e-4\n");
  const AtmosphereTable table = read_atmprof_text(text, "t.dat");
  CHECK(table.profile_id == custom_profile_id);
  CHECK(table.file_name == "t.dat");
  CHECK(table.alt_km == std::vector<double>({0.0, 1.0, 5.0}));
  CHECK(table.rho == std::vector<double>({1.2, 1.0, 0.5}));
  CHECK(table.thick == std::vector<double>({1000.0, 900.0, 500.0}));
  CHECK(table.refidx_m1 == std::vector<double>({2.9e-4, 2.7e-4, 2.2e-4}));

  const auto parse = [](const char *content) {
    std::istringstream in(content);
    return read_atmprof_text(in, "bad.dat");
  };
  CHECK(throws_error([&] { parse("0 1.2 1000\n1 1.0 900 2.7e-4\n"); }));  // 3 columns
  CHECK(throws_error([&] { parse("0 1.2 x 2.9e-4\n1 1 900 2e-4\n"); }));  // not a number
  CHECK(throws_error([&] { parse("0 1.2 1000 nan\n1 1 900 2e-4\n"); }));  // not finite
  CHECK(throws_error([&] { parse("1 1.2 1000 2e-4\n0 1 900 2e-4\n"); })); // descending
  CHECK(throws_error([&] { parse("0 1.2 1000 2e-4\n"); }));               // one level
  CHECK(throws_error([&] { parse("0 -1 1000 2e-4\n1 1 900 2e-4\n"); }));  // negative density
  std::string many;
  for (int i = 0; i < 121; ++i) {
    many += std::to_string(i) + " 1 1 1e-4\n";
  }
  CHECK(throws_error([&] { parse(many.c_str()); })); // 121 levels
}

// Without CORSIKA headers the table is written as ATMPROF 99 before the first EVTH, with the
// run's observation level; EVTH word 77 names the same profile.
void test_atmprof_synthesised() {
  SimtelOptions options;
  options.atmprof = three_level_table();
  std::ostringstream warnings;
  options.warnings = &warnings;
  const std::vector<EventioBlock> blocks =
      read_blocks(run_adapter_with(light_source_records(), options));
  std::vector<std::uint32_t> types;
  for (const EventioBlock &block : blocks) {
    types.push_back(block.type);
  }
  const std::vector<std::uint32_t> expected = {iact::runh, iact::telpos,   iact::atmprof,
                                               iact::evth, iact::teloff,   iact::telarray,
                                               iact::evte, iact::inputcfg, iact::rune};
  CHECK(types == expected);
  if (types != expected) {
    return;
  }
  const EventioBlock &atmprof = blocks[2];
  CHECK(atmprof.ident == custom_profile_id);
  const AtmosphereTable table = decode_atmprof(atmprof);
  CHECK(table.observation_level_cm == 220000.0); // run metadata, 2200 m
  CHECK(table.alt_km == std::vector<double>({0.0, 1.0, 5.0}));
  CHECK(table.refidx_m1 == std::vector<double>({2.9e-4, 2.7e-4, 2.2e-4}));
  const std::vector<float> evth = decode_float_block(blocks[3]);
  const auto flags = static_cast<std::int64_t>(std::llround(evth[76]));
  CHECK(((flags >> 10) & 0x3FF) == custom_profile_id); // profile of EVTH word 77
  CHECK(warnings.str().find("EVTH") != std::string::npos);
  CHECK(warnings.str().find("ATMPROF") == std::string::npos); // table given: no note
}

// Without a table the output is unchanged, and each synthesised header is reported.
void test_synthesised_warnings() {
  std::ostringstream warnings;
  SimtelOptions options;
  options.warnings = &warnings;
  const std::vector<EventioBlock> blocks =
      read_blocks(run_adapter_with(light_source_records(), options));
  for (const EventioBlock &block : blocks) {
    CHECK(block.type != iact::atmprof);
  }
  const std::string text = warnings.str();
  for (const char *name : {"RUNH", "TELPOS", "EVTH", "TELOFF", "EVTE", "INPUTCFG", "RUNE"}) {
    CHECK(text.find(std::string("no CORSIKA ") + name + " block") != std::string::npos);
  }
  CHECK(text.find("--atmprof FILE") != std::string::npos);
  // The same stream without a warning stream gives the same bytes.
  CHECK(run_adapter_with(light_source_records(), SimtelOptions{}) ==
        run_adapter(light_source_records(), {}));
}

// A CORSIKA EVTH names its own profile: the ATMPROF block takes that number and comes before it.
void test_atmprof_takes_evth_profile() {
  std::vector<Record> records = light_source_records();
  records.insert(records.begin() + 2, blob_of(evth_block(3.0F)));
  records.insert(records.begin() + 2, blob_of(teloff_block()));
  SimtelOptions options;
  options.atmprof = three_level_table();
  const std::vector<EventioBlock> blocks = read_blocks(run_adapter_with(records, options));
  std::size_t atmprof_at = blocks.size();
  std::size_t evth_at = blocks.size();
  for (std::size_t i = 0; i < blocks.size(); ++i) {
    if (blocks[i].type == iact::atmprof) {
      atmprof_at = i;
    }
    if (blocks[i].type == iact::evth) {
      evth_at = i;
    }
  }
  CHECK(atmprof_at < evth_at);
  if (atmprof_at < blocks.size()) {
    CHECK(blocks[atmprof_at].ident == 5); // evth_block() names atmosphere 5 (word 77)
  }
}

// A stream with its own ATMPROF block cannot take --atmprof as well.
void test_atmprof_option_conflicts_with_stream() {
  std::vector<Record> records = light_source_records();
  records.insert(records.begin() + 1, blob_of(make_atmprof_block(three_level_table())));
  SimtelOptions options;
  options.atmprof = three_level_table();
  CHECK(throws_error([&] { run_adapter_with(records, options); }));
  // Without the option the block passes through unchanged.
  CHECK(!throws_error([&] { run_adapter_with(records, SimtelOptions{}); }));
}

// CORSIKA output carries all headers: no warnings, and the bytes are the input's.
void test_no_warnings_with_corsika_headers() {
  const std::string input = corsika_file(1);
  std::ostringstream stream_bytes;
  {
    Writer writer{stream_bytes, "test_simtel"};
    std::istringstream in(input);
    corsika_to_records(in, [&writer](const Record &record) { write_record(writer, record); });
    writer.finish();
  }
  std::istringstream stream_in(stream_bytes.str());
  std::ostringstream output;
  std::ostringstream warnings;
  stream_to_simtel(stream_in, output, SimtelOptions{CameraMap{}, std::nullopt, &warnings, {}, {}});
  CHECK(output.str() == input);
  CHECK(warnings.str().empty());
}

void test_synthesised_headers() {
  const std::vector<Record> records = {
      record_of(RecordType::run, light_source_run()),
      record_of(RecordType::event, Metadata{{"event_number", std::int64_t{3}}}),
      record_of(RecordType::telescope,
                Metadata{{"telescope_id", std::int64_t{1}}, {"array_index", std::int64_t{0}}}),
      photons_of(chunk_of(1, Frame::ground, {0.5, 1.5}, {0.2, 0.0}, {0.1, 0.0}, {-0.2, 0.0})),
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{1}}}),
      record_of(RecordType::event_end),
      record_of(RecordType::run_end),
  };
  const std::vector<EventioBlock> blocks = read_blocks(run_adapter(records, {}));
  std::vector<std::uint32_t> types;
  for (const EventioBlock &block : blocks) {
    types.push_back(block.type);
  }
  const std::vector<std::uint32_t> expected = {iact::runh,     iact::telpos,   iact::evth,
                                               iact::teloff,   iact::telarray, iact::evte,
                                               iact::inputcfg, iact::rune};
  CHECK(types == expected);
  if (types == expected) {
    const std::vector<float> runh = decode_float_block(blocks[0]);
    CHECK_NEAR(runh[1], 7.0, 0.0);       // run number, word 2
    CHECK_NEAR(runh[5], 220000.0, 1e-3); // observation level, word 6 [cm]
    const TelPos pos = decode_telpos(blocks[1]);
    CHECK(pos.x.size() == 2);
    CHECK_NEAR(pos.x[1], 10000.0, 1e-3); // second telescope at 100 m
    CHECK_NEAR(pos.r[0], 1200.0, 1e-3);
    const std::vector<float> evth = decode_float_block(blocks[2]);
    CHECK_NEAR(evth[1], 3.0, 0.0);       // event number, word 2
    CHECK_NEAR(evth[43], 7.0, 0.0);      // run number, word 44
    CHECK_NEAR(evth[4], 220000.0, 1e-3); // word 5 [cm]
    const std::vector<float> evte = decode_float_block(blocks[5]);
    CHECK_NEAR(evte[1], 3.0, 0.0);
    // RUNE body: count 3, "RUNE", run number, number of showers.
    const std::vector<std::uint8_t> &body = blocks[7].body;
    CHECK(body.size() == 16);
    if (body.size() == 16) {
      CHECK(std::string(body.begin() + 4, body.begin() + 8) == "RUNE");
      float run_value = 0.0F;
      float showers = 0.0F;
      std::memcpy(&run_value, body.data() + 8, 4);
      std::memcpy(&showers, body.data() + 12, 4);
      CHECK_NEAR(run_value, 7.0, 0.0);
      CHECK_NEAR(showers, 1.0, 0.0);
    }
  }
}

// Records of a light-source run with `events` events of one telescope each.
std::vector<Record> light_source_events(int events) {
  std::vector<Record> records = {record_of(RecordType::run, light_source_run())};
  for (int i = 1; i <= events; ++i) {
    records.push_back(record_of(RecordType::event, Metadata{{"event_number", std::int64_t{i}}}));
    records.push_back(record_of(RecordType::telescope, Metadata{{"telescope_id", std::int64_t{1}},
                                                                {"array_index", std::int64_t{0}}}));
    records.push_back(
        photons_of(chunk_of(1, Frame::ground, {0.5, 1.5}, {0.2, 0.0}, {0.1, 0.0}, {-0.2, 0.0})));
    records.push_back(
        record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{1}}}));
    records.push_back(record_of(RecordType::event_end));
  }
  records.push_back(record_of(RecordType::run_end));
  return records;
}

// The synthesised run end counts the events written, and the input line NSHOW
// that precedes it gives the same number; a CORSIKA input block is not duplicated.
void test_synthesised_run_end_counts_events() {
  for (const int events : {1, 4}) {
    const std::vector<EventioBlock> blocks =
        read_blocks(run_adapter(light_source_events(events), {}));
    const EventioBlock &rune = blocks.back();
    const EventioBlock &inputcfg = blocks[blocks.size() - 2];
    CHECK(rune.type == iact::rune);
    CHECK(inputcfg.type == iact::inputcfg);
    if (rune.type != iact::rune || inputcfg.type != iact::inputcfg) {
      return;
    }
    // RUNE: count 3, "RUNE", run number, number of showers.
    CHECK(rune.body.size() == 16);
    float showers = 0.0F;
    std::memcpy(&showers, rune.body.data() + 12, 4);
    CHECK_NEAR(showers, static_cast<double>(events), 0.0);
    // INPUTCFG: int32 line count, then int16 length and text per line.
    const std::vector<std::uint8_t> &body = inputcfg.body;
    CHECK(body.size() >= 6);
    std::int32_t lines = 0;
    std::int16_t length = 0;
    std::memcpy(&lines, body.data(), 4);
    std::memcpy(&length, body.data() + 4, 2);
    CHECK(lines == 1);
    const std::string text(body.begin() + 6, body.begin() + 6 + length);
    CHECK(text == "NSHOW " + std::to_string(events));
  }
  // A stream with its own INPUTCFG keeps it and gets no second one.
  std::vector<Record> records = light_source_events(2);
  std::vector<std::uint8_t> body = {1, 0, 0, 0, 7, 0};
  const std::string line = "NSHOW 2";
  body.insert(body.end(), line.begin(), line.end());
  records.insert(records.begin() + 1, blob_of(make_block(iact::inputcfg, 0, 0, body)));
  int count = 0;
  for (const EventioBlock &block : read_blocks(run_adapter(records, {}))) {
    count += block.type == iact::inputcfg ? 1 : 0;
  }
  CHECK(count == 1);
}

// A stream without the metadata sim_telarray needs is rejected, not padded with defaults.
void test_missing_metadata_rejected() {
  Metadata run;
  run["run_number"] = std::int64_t{1};
  const std::vector<Record> records = {
      record_of(RecordType::run, run),
      record_of(RecordType::event, Metadata{{"event_number", std::int64_t{1}}}),
  };
  CHECK(throws_error([&] { run_adapter(records, {}); }));
}

// Sky and telescope frames must not reach sim_telarray.
void test_frames_rejected() {
  const std::vector<Record> sky = {
      record_of(RecordType::run, light_source_run()),
      record_of(RecordType::event, Metadata{{"event_number", std::int64_t{1}}}),
      record_of(RecordType::telescope,
                Metadata{{"telescope_id", std::int64_t{1}}, {"array_index", std::int64_t{0}}}),
      photons_of(chunk_of(1, Frame::sky, {0.5}, {0.2}, {0.1}, {-0.2})),
  };
  CHECK(throws_error([&] { run_adapter(sky, {}); }));
  std::vector<Record> telescope_frame = sky;
  telescope_frame.back() = photons_of(chunk_of(1, Frame::telescope, {0.5}, {0.2}, {0.1}, {-0.2}));
  CHECK(throws_error([&] { run_adapter(telescope_frame, {}); }));
}

Metadata efficiency_status_of(const std::string &factor, const std::string &state) {
  return Metadata{{"factors", factor}, {factor + ".state", state}};
}

// One event with telescope 2 in `frame`; the efficiency status records come before the event (run)
// or after the telescope record (telescope).
std::vector<Record>
efficiency_status_run(Frame frame, const std::vector<Metadata> &run_efficiency_statuses,
                      const std::vector<Metadata> &telescope_efficiency_statuses) {
  std::vector<Record> records = {record_of(RecordType::run, light_source_run())};
  for (const Metadata &efficiency_status : run_efficiency_statuses) {
    records.push_back(record_of(RecordType::efficiency_status, efficiency_status));
  }
  records.push_back(record_of(RecordType::event, Metadata{{"event_number", std::int64_t{1}}}));
  records.push_back(record_of(RecordType::telescope, Metadata{{"telescope_id", std::int64_t{2}},
                                                              {"array_index", std::int64_t{0}}}));
  for (const Metadata &efficiency_status : telescope_efficiency_statuses) {
    records.push_back(record_of(RecordType::efficiency_status, efficiency_status));
  }
  records.push_back(photons_of(chunk_of(2, frame, {0.5}, {0.2}, {0.1}, {-0.2})));
  records.push_back(
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{2}}}));
  records.push_back(record_of(RecordType::event_end));
  records.push_back(record_of(RecordType::run_end));
  return records;
}

// Efficiency status rules (docs/usage.md): which factors sim_telarray may leave to
// itself in each frame, and which ones must be applied upstream.
void test_efficiency_status_rules() {
  const CameraMap camera{{2, CameraTelescope{1, 0.0}}};
  const std::set<std::string> constant_allowed{"telescope_transmission_constant"};
  const std::string constant = "telescope_transmission_constant";

  // Pending constant part of the telescope transmission (angular mode) in camera frame.
  const auto constant_run =
      efficiency_status_run(Frame::camera, {efficiency_status_of(constant, "pending")}, {});
  CHECK(throws_error([&] { run_adapter(constant_run, camera); }));
  CHECK(!throws_error([&] { run_adapter(constant_run, camera, constant_allowed); }));
  // Ground frame: the same factor is not in the sim_telarray list either.
  CHECK(throws_error([&] {
    run_adapter(
        efficiency_status_run(Frame::ground, {efficiency_status_of(constant, "pending")}, {}), {});
  }));

  // sim_telarray applies the camera efficiencies in camera frame: pending is accepted.
  CHECK(!throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {efficiency_status_of("qe", "pending")}, {}),
                camera);
  }));
  // An envelope on a factor sim_telarray applies in full would make the result too dim.
  CHECK(throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {efficiency_status_of("qe", "envelope")}, {}),
                camera);
  }));
  // Ground frame: atmosphere is accepted in any state (sim_telarray's table must then be zero).
  CHECK(!throws_error([&] {
    run_adapter(
        efficiency_status_run(Frame::ground, {efficiency_status_of("atmosphere", "applied")}, {}),
        {});
  }));
  CHECK(!throws_error([&] {
    run_adapter(
        efficiency_status_run(Frame::ground, {efficiency_status_of("atmosphere", "pending")}, {}),
        {});
  }));
  // Mirror reflectivity is applied by sim_telarray in ground frame only.
  CHECK(!throws_error([&] {
    run_adapter(efficiency_status_run(Frame::ground,
                                      {efficiency_status_of("mirror_reflectivity", "pending")}, {}),
                {});
  }));
  CHECK(throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera,
                                      {efficiency_status_of("mirror_reflectivity", "pending")}, {}),
                camera);
  }));

  // Applied upstream and not applied by sim_telarray in the frame: accepted.
  CHECK(!throws_error([&] {
    run_adapter(
        efficiency_status_run(Frame::camera, {efficiency_status_of("atmosphere", "applied")}, {}),
        camera);
  }));
  // dlufd assigns wavelengths and may add scattered photons; both are upstream-only factors.
  for (const char *factor : {"wavelength", "scattering"}) {
    CHECK(!throws_error([&] {
      run_adapter(
          efficiency_status_run(Frame::ground, {efficiency_status_of(factor, "applied")}, {}), {});
    }));
    CHECK(throws_error([&] {
      run_adapter(
          efficiency_status_run(Frame::ground, {efficiency_status_of(factor, "pending")}, {}), {});
    }));
  }
  // Status records describe ordered stage updates: pending -> envelope -> applied is valid;
  // moving backwards is not.
  CHECK(!throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera,
                                      {efficiency_status_of("wavelength", "pending"),
                                       efficiency_status_of("wavelength", "applied")},
                                      {}),
                camera);
  }));
  CHECK(throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera,
                                      {efficiency_status_of("wavelength", "applied"),
                                       efficiency_status_of("wavelength", "pending")},
                                      {}),
                camera);
  }));
  CHECK(!throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {},
                                      {efficiency_status_of("qe", "pending"),
                                       efficiency_status_of("qe", "envelope")}),
                camera, {}, {"qe"});
  }));
  CHECK(throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {},
                                      {efficiency_status_of("qe", "envelope"),
                                       efficiency_status_of("qe", "pending")}),
                camera, {}, {"qe"});
  }));
  // Applied upstream and applied again by sim_telarray: rejected.
  CHECK(throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {efficiency_status_of("qe", "applied")}, {}),
                camera);
  }));
  // Atmosphere pending is not applied by sim_telarray: rejected.
  CHECK(throws_error([&] {
    run_adapter(
        efficiency_status_run(Frame::camera, {efficiency_status_of("atmosphere", "pending")}, {}),
        camera);
  }));

  // Unknown factor names are rejected, also when applied; --allow-pending accepts them.
  CHECK(throws_error([&] {
    run_adapter(
        efficiency_status_run(Frame::camera, {efficiency_status_of("made_up", "applied")}, {}),
        camera);
  }));
  CHECK(!throws_error([&] {
    run_adapter(
        efficiency_status_run(Frame::camera, {efficiency_status_of("made_up", "applied")}, {}),
        camera, {"made_up"});
  }));
  CHECK(throws_error(
      [&] { run_adapter(efficiency_status_run(Frame::camera, {}, {}), camera, {}, {"qe"}); }));

  // Scaled factors (sim_telarray table holds f/f_env): only the envelope state of a factor
  // that sim_telarray applies in the frame is accepted.
  CHECK(throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {efficiency_status_of("qe", "envelope")}, {}),
                camera);
  }));
  CHECK(!throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {efficiency_status_of("qe", "envelope")}, {}),
                camera, {}, {"qe"});
  }));
  CHECK(throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {efficiency_status_of("qe", "pending")}, {}),
                camera, {}, {"qe"});
  }));
  CHECK(throws_error([&] {
    run_adapter(efficiency_status_run(
                    Frame::camera, {efficiency_status_of("mirror_reflectivity", "envelope")}, {}),
                camera, {}, {"mirror_reflectivity"});
  }));
  // Ground frame: the atmosphere cannot be a scaled factor, in any state, although the
  // unscaled atmosphere is accepted in any state.
  for (const char *state : {"applied", "pending", "envelope"}) {
    CHECK(throws_error([&] {
      run_adapter(
          efficiency_status_run(Frame::ground, {efficiency_status_of("atmosphere", state)}, {}), {},
          {}, {"atmosphere"});
    }));
  }

  // Malformed efficiency status records: no factors list, no state, unknown state.
  CHECK(throws_error([&] {
    run_adapter(
        efficiency_status_run(Frame::camera, {Metadata{{"qe.state", std::string("pending")}}}, {}),
        camera);
  }));
  CHECK(throws_error([&] {
    run_adapter(
        efficiency_status_run(Frame::camera, {Metadata{{"factors", std::string("qe")}}}, {}),
        camera);
  }));
  CHECK(throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {efficiency_status_of("qe", "maybe")}, {}),
                camera);
  }));

  // A telescope efficiency status overrides the run efficiency status for its factor.
  CHECK(!throws_error([&] {
    run_adapter(efficiency_status_run(Frame::camera, {efficiency_status_of(constant, "pending")},
                                      {efficiency_status_of(constant, "applied")}),
                camera);
  }));

  // A telescope efficiency status after the photons of its telescope is rejected.
  std::vector<Record> late = efficiency_status_run(Frame::camera, {}, {});
  late.insert(late.end() - 3,
              record_of(RecordType::efficiency_status, efficiency_status_of("qe", "pending")));
  CHECK(throws_error([&] { run_adapter(late, camera); }));
}

// CORSIKA-style file -> corsika_to_records -> Writer -> Reader -> SimtelWriter gives the
// input back byte for byte, for positions that round-trip cm -> m -> cm.
void test_round_trip_byte_identical() {
  const std::string input = corsika_file(2);
  std::istringstream in(input);
  std::ostringstream stream_bytes;
  {
    Writer writer{stream_bytes, "test_simtel"};
    corsika_to_records(in, [&writer](const Record &record) { write_record(writer, record); });
    writer.finish();
  }
  std::istringstream stream_in(stream_bytes.str());
  Reader reader{stream_in};
  std::ostringstream output;
  SimtelWriter sim{output};
  Record record;
  while (reader.next(record)) {
    sim.write(record);
  }
  sim.finish();
  CHECK(reader.ended());
  CHECK(output.str() == input);
  CHECK(!input.empty());

  // Fraction of random float positions that survive cm -> m -> cm (not a pass criterion).
  std::mt19937 rng{20261107};
  std::uniform_real_distribution<float> uniform{-1.0e5F, 1.0e5F};
  int survived = 0;
  const int trials = 100000;
  for (int i = 0; i < trials; ++i) {
    survived += round_trips(uniform(rng)) ? 1 : 0;
  }
  std::printf("round trip cm->m->cm: %.2f%% of %d random float positions in +-1e5 cm\n",
              100.0 * survived / trials, trials);
}

} // namespace

int main() {
  return run_checks([] {
    test_camera_map_parsing();
    test_invalid_camera_map();
    test_camera_conversion();
    test_camera_frame_through_adapter();
    test_telpos_zeroed_for_camera_only();
    test_telescope_frame_through_adapter();
    test_telescope_frame_report_and_telpos();
    test_telescope_frame_rejected();
    test_synthesised_headers();
    test_synthesised_run_end_counts_events();
    test_missing_metadata_rejected();
    test_frames_rejected();
    test_efficiency_status_rules();
    test_round_trip_byte_identical();
    test_atmprof_text();
    test_atmprof_synthesised();
    test_synthesised_warnings();
    test_atmprof_takes_evth_profile();
    test_atmprof_option_conflicts_with_stream();
    test_no_warnings_with_corsika_headers();
  });
}
