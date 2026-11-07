#pragma once

// Adapters between the photon stream and the CORSIKA IACT EventIO world:
// CORSIKA output to a photon stream (corsika2ps) and a photon stream to
// sim_telarray input (ps2simtel). Header-only, standard library only.
//
// Camera-frame (C) photons are converted per telescope to the sim_telarray
// bypass convention exactly as simtools' camera_bunches does
// (simtools/src/simtools/simtel/photon_bunch_file.py:105-156):
//   - single mirror (Bypass_Optics=1): (x, y) -> (x, -y), (cx, cy) -> (cx, -cy);
//   - dual mirror (Bypass_Optics=2): the camera rotation by camera_rotation_deg,
//     x' = x cos r + y sin r, y' = y cos r - x sin r (same for cx, cy).
// The bypass value is an option of the sim_telarray run and is printed by
// ps2simtel for each converted telescope.
//
// Telescopes with C photons get z = 0 in TELPOS (1201). sim_telarray moves each
// bunch from the observation level to the telescope by ztel * airmass along its
// direction (sim_telarray.c:2740-2790), which applies to ground-frame bunches
// only. Camera-frame bunches are already at the camera, so the shift must be
// zero (simtools shower_route._zero_traced_telescope_heights, :594).
//
// A telescope is either camera-frame (listed in the camera map) or ground-frame
// (not listed). Mixing the two for one telescope is rejected, because the TELPOS
// height is per telescope for the whole run.
//
// Streams without CORSIKA blobs (e.g. from light sources) get the sim_telarray
// headers RUNH, TELPOS, EVTH, TELOFF, EVTE and RUNE synthesised from the run and
// event metadata with the words simtools photon_bunch_file.py fills
// (write_photon_bunches, :379-466; _EVENT_WORDS :49-98). ATMPROF is not written:
// the stream does not carry an atmospheric profile.

#include "gstrickd/corsika.hpp"
#include "gstrickd/eventio.hpp"
#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/stream.hpp"

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <limits>
#include <map>
#include <numbers>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace gstrickd {

struct CameraTelescope {
  int mirror_class{1};      // 1: single mirror (Bypass_Optics=1), 2: dual mirror (Bypass_Optics=2)
  double rotation_deg{0.0}; // camera rotation, used for dual mirrors only
};

// Camera-frame telescopes by telescope id (1-based, sim_telarray numbering).
using CameraMap = std::map<std::uint32_t, CameraTelescope>;

// sim_telarray Bypass_Optics value for a mirror class (simtools bypass_optics_level).
[[nodiscard]] inline int bypass_optics_level(int mirror_class) {
  if (mirror_class != 1 && mirror_class != 2) {
    throw std::invalid_argument("mirror class must be 1 (single) or 2 (dual)");
  }
  return mirror_class == 2 ? 2 : 1;
}

// Parses the camera map: one line per telescope, "telescope_id mirror_class
// camera_rotation_deg"; '#' starts a comment. Throws std::invalid_argument with
// the line number on malformed input or a repeated telescope.
namespace simtel_detail {

// Whole-string number parsers; false instead of an exception.
inline bool parse_integer(const std::string &text, long long &value) {
  try {
    std::size_t used = 0;
    value = std::stoll(text, &used);
    return used == text.size();
  } catch (const std::exception &) {
    return false;
  }
}

inline bool parse_real(const std::string &text, double &value) {
  try {
    std::size_t used = 0;
    value = std::stod(text, &used);
    return used == text.size() && std::isfinite(value);
  } catch (const std::exception &) {
    return false;
  }
}

} // namespace simtel_detail

inline CameraMap read_camera_map(std::istream &in) {
  CameraMap map;
  std::string line;
  std::size_t number = 0;
  const auto fail = [&number](const std::string &what) {
    return std::invalid_argument("camera map line " + std::to_string(number) + ": " + what);
  };
  while (std::getline(in, line)) {
    ++number;
    const auto hash = line.find('#');
    if (hash != std::string::npos) {
      line.erase(hash);
    }
    std::istringstream fields(line);
    std::string id_text, class_text, rotation_text, extra;
    if (!(fields >> id_text)) {
      continue; // blank or comment line
    }
    if (!(fields >> class_text >> rotation_text) || (fields >> extra)) {
      throw fail("expected 'telescope_id mirror_class camera_rotation_deg'");
    }
    long long id_value = 0;
    long long class_value = 0;
    double rotation = 0.0;
    if (!simtel_detail::parse_integer(id_text, id_value) || id_value < 1 ||
        id_value > std::numeric_limits<std::uint32_t>::max()) {
      throw fail("telescope_id must be an integer of 1 or larger");
    }
    if (!simtel_detail::parse_integer(class_text, class_value) ||
        (class_value != 1 && class_value != 2)) {
      throw fail("mirror_class must be 1 (single) or 2 (dual)");
    }
    if (!simtel_detail::parse_real(rotation_text, rotation)) {
      throw fail("camera_rotation_deg must be a finite number");
    }
    const auto id = static_cast<std::uint32_t>(id_value);
    if (!map.emplace(id, CameraTelescope{static_cast<int>(class_value), rotation}).second) {
      throw fail("telescope " + std::to_string(id) + " listed twice");
    }
  }
  return map;
}

// Converts camera-frame positions and directions in place to the sim_telarray
// bypass convention (metres, as in the stream). The frame label is not changed.
// Arithmetic in double, rounded to float once, as simtools does in float64.
inline void camera_to_bypass_convention(PhotonChunk &chunk, const CameraTelescope &camera) {
  if (!(chunk.has(Column::x) && chunk.has(Column::y) && chunk.has(Column::cx) &&
        chunk.has(Column::cy))) {
    throw std::invalid_argument("camera-frame chunk lacks x, y, cx or cy");
  }
  if (!std::isfinite(camera.rotation_deg)) {
    throw std::invalid_argument("camera rotation must be finite");
  }
  const int level = bypass_optics_level(camera.mirror_class);
  std::vector<float> &x = chunk[Column::x];
  std::vector<float> &y = chunk[Column::y];
  std::vector<float> &cx = chunk[Column::cx];
  std::vector<float> &cy = chunk[Column::cy];
  if (level == 1) {
    for (std::size_t i = 0; i < x.size(); ++i) {
      y[i] = -y[i];
      cy[i] = -cy[i];
    }
    return;
  }
  const double angle = camera.rotation_deg * std::numbers::pi / 180.0;
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double xi = x[i];
    const double yi = y[i];
    const double cxi = cx[i];
    const double cyi = cy[i];
    x[i] = static_cast<float>(xi * c + yi * s);
    y[i] = static_cast<float>(yi * c - xi * s);
    cx[i] = static_cast<float>(cxi * c + cyi * s);
    cy[i] = static_cast<float>(cyi * c - cxi * s);
  }
}

namespace simtel_detail {

// Words of the simtools header blocks (1-based CORSIKA word numbers).
using Words = std::vector<std::pair<int, double>>;

inline constexpr std::size_t block_words = 273;
inline constexpr double energy_gev = 1e-9; // simtools _ENERGY_GEV
inline constexpr double date = 261022.0;   // simtools _DATE
inline constexpr double corsika_version = 6.999;
inline constexpr std::int64_t profile_id = 99;   // simtools _PROFILE_ID
inline constexpr std::int64_t event_flags = 187; // simtools _EVENT_FLAGS

// simtools _EVENT_WORDS: constant words of the event header.
inline const Words event_constants = {{3, 1.0},
                                      {4, energy_gev},
                                      {10, -1.0},
                                      {12, std::numbers::pi},
                                      {46, corsika_version},
                                      {45, date},
                                      {47, 1.0},
                                      {59, energy_gev},
                                      {60, energy_gev},
                                      {61, 1.0},
                                      {62, 0.106},
                                      {63, 0.000512},
                                      {64, energy_gev},
                                      {71, 20.0},
                                      {75, 2.0},
                                      {76, 3.0},
                                      {80, 3.0},
                                      {83, 180.0},
                                      {84, 180.0},
                                      {85, 1.0},
                                      {86, 1.0},
                                      {87, 1.0},
                                      {95, 1.0},
                                      {96, 200.0},
                                      {97, 1000.0},
                                      {98, 1.0}};

// Parses a raw top-level EventIO block (as stored in a blob record).
inline EventioBlock parse_block(const std::vector<std::uint8_t> &raw) {
  std::istringstream in(std::string(raw.begin(), raw.end()));
  EventioReader reader{in};
  EventioBlock block;
  if (!reader.next(block)) {
    throw std::runtime_error("gstrickd: blob holds no EventIO block");
  }
  return block;
}

// Bytes of a top-level block as written to a file (sync tag and header included).
inline std::vector<std::uint8_t> block_bytes(const EventioBlock &block) {
  std::ostringstream out;
  write_block(out, block);
  const std::string text = out.str();
  return {text.begin(), text.end()};
}

// A 273-word float block with the block name in word 0 (write_tel_block layout).
inline EventioBlock named_block(std::uint32_t type, std::uint32_t ident, const char *name,
                                const Words &words) {
  std::vector<float> values(block_words, 0.0F);
  std::uint32_t name_bits = 0;
  for (std::size_t i = 0; i < 4; ++i) {
    name_bits |= std::uint32_t{static_cast<std::uint8_t>(name[i])} << (8 * i);
  }
  values[0] = std::bit_cast<float>(name_bits);
  for (const auto &[word, value] : words) {
    values[static_cast<std::size_t>(word - 1)] = static_cast<float>(value);
  }
  return make_float_block(type, ident, values);
}

inline double number_field(const Metadata &metadata, const char *key) {
  const auto it = metadata.find(key);
  if (it == metadata.end()) {
    throw std::runtime_error(std::string("stream lacks '") + key +
                             "', needed for sim_telarray headers of a stream without CORSIKA "
                             "blobs");
  }
  if (const auto *integer = std::get_if<std::int64_t>(&it->second)) {
    return static_cast<double>(*integer);
  }
  if (const auto *real = std::get_if<double>(&it->second)) {
    return *real;
  }
  throw std::runtime_error(std::string("stream key '") + key + "' is not a number");
}

inline std::vector<double> array_field(const Metadata &metadata, const char *key) {
  const auto it = metadata.find(key);
  if (it == metadata.end()) {
    throw std::runtime_error(std::string("stream lacks '") + key +
                             "', needed for sim_telarray TELPOS");
  }
  const auto *values = std::get_if<std::vector<double>>(&it->second);
  if (values == nullptr) {
    throw std::runtime_error(std::string("stream key '") + key + "' is not a float array");
  }
  return *values;
}

// Telescope positions and radii from metres (run metadata) to the cm of TELPOS.
inline TelPos telpos_from_metres(const std::vector<double> &x, const std::vector<double> &y,
                                 const std::vector<double> &z, const std::vector<double> &r) {
  const std::size_t n = x.size();
  if (y.size() != n || z.size() != n || r.size() != n) {
    throw std::runtime_error("stream telescope position arrays differ in length");
  }
  TelPos pos;
  for (std::size_t i = 0; i < n; ++i) {
    pos.x.push_back(100.0 * x[i]);
    pos.y.push_back(100.0 * y[i]);
    pos.z.push_back(100.0 * z[i]);
    pos.r.push_back(100.0 * r[i]);
  }
  return pos;
}

} // namespace simtel_detail

// Writes a photon stream as sim_telarray input (EventIO, normal mode: one 1204
// container per array). Records are read in stream order; see the file comment
// for the frame and header rules.
class PsToSimtel {
public:
  PsToSimtel(std::ostream &out, CameraMap camera) : writer_{out}, camera_{std::move(camera)} {
    for (const auto &[id, telescope] : camera_) {
      if (id == 0) {
        throw std::invalid_argument("camera map telescope id must be at least 1");
      }
      static_cast<void>(bypass_optics_level(telescope.mirror_class));
      if (!std::isfinite(telescope.rotation_deg)) {
        throw std::invalid_argument("camera rotation must be finite");
      }
    }
  }

  void write(const Record &record) {
    switch (record.type) {
    case RecordType::stream:
    case RecordType::stream_end:
    case RecordType::ledger:
      return; // not part of sim_telarray input
    case RecordType::run:
      end_run();
      begin_run(record.metadata);
      return;
    case RecordType::event:
      ensure_run_headers();
      begin_event(record.metadata);
      return;
    case RecordType::telescope:
      require_event("telescope record");
      ensure_event_headers();
      writer_.write(record);
      return;
    case RecordType::photons:
      photons(record.chunk);
      return;
    case RecordType::telescope_end:
      writer_.write(record);
      return;
    case RecordType::event_end:
      end_event();
      return;
    case RecordType::run_end:
      end_run();
      return;
    case RecordType::blob:
      blob(record);
      return;
    }
  }

  // Closes the open run and array container, and writes the stream's end.
  void finish() {
    end_run();
    writer_.finish();
  }

private:
  void photons(const PhotonChunk &chunk) {
    switch (chunk.frame) {
    case Frame::sky:
    case Frame::telescope:
      throw std::runtime_error(
          "photons in the sky or telescope frame cannot go to sim_telarray: convert them to "
          "the ground or camera frame first");
    case Frame::ground:
      if (camera_.contains(chunk.telescope_id)) {
        throw std::runtime_error("telescope " + std::to_string(chunk.telescope_id) +
                                 " is in the camera map but has ground-frame photons");
      }
      writer_.write_chunk(chunk);
      return;
    case Frame::camera: {
      const auto it = camera_.find(chunk.telescope_id);
      if (it == camera_.end()) {
        throw std::runtime_error("telescope " + std::to_string(chunk.telescope_id) +
                                 " has camera-frame photons but no camera map entry");
      }
      PhotonChunk converted = chunk;
      camera_to_bypass_convention(converted, it->second);
      converted.frame = Frame::ground; // sim_telarray input is written as ground-frame bunches
      writer_.write_chunk(converted);
      return;
    }
    }
  }

  void blob(const Record &record) {
    switch (record.blob_eventio_type) {
    case iact::runh:
      have_runh_ = true;
      break;
    case iact::telpos:
      have_telpos_ = true;
      if (!camera_.empty()) {
        const EventioBlock telpos = simtel_detail::parse_block(record.blob);
        TelPos pos = decode_telpos(telpos);
        for (const auto &entry : camera_) {
          if (entry.first == 0 || entry.first > pos.z.size()) {
            throw std::runtime_error("camera map telescope " + std::to_string(entry.first) +
                                     " is not in the TELPOS block");
          }
          pos.z[entry.first - 1] = 0.0;
        }
        Record zeroed = record;
        zeroed.blob = simtel_detail::block_bytes(make_telpos_block(pos));
        writer_.write(zeroed);
        return;
      }
      break;
    case iact::evth:
      have_evth_ = true;
      break;
    case iact::teloff:
      have_teloff_ = true;
      break;
    case iact::evte:
      have_evte_ = true;
      break;
    case iact::rune:
      have_rune_ = true;
      break;
    default:
      break;
    }
    writer_.write(record);
  }

  void require_event(const char *what) const {
    if (!event_open_) {
      throw std::runtime_error(std::string(what) + " outside an event");
    }
  }

  void begin_run(const Metadata &metadata) {
    run_open_ = true;
    run_metadata_ = metadata;
    run_headers_done_ = false;
    have_runh_ = have_telpos_ = have_rune_ = false;
    events_ = 0;
  }

  void begin_event(const Metadata &metadata) {
    event_open_ = true;
    event_headers_done_ = false;
    have_evth_ = have_teloff_ = have_evte_ = false;
    ++events_;
    const auto it = metadata.find("event_number");
    event_number_ = it != metadata.end() && std::holds_alternative<std::int64_t>(it->second)
                        ? std::get<std::int64_t>(it->second)
                        : static_cast<std::int64_t>(events_);
  }

  // Run headers come before the first event; CORSIKA blobs have arrived by then.
  void ensure_run_headers() {
    if (!run_open_ || run_headers_done_) {
      return;
    }
    run_headers_done_ = true;
    if (!have_runh_) {
      emit(run_header());
    }
    if (!have_telpos_) {
      emit(telescope_positions());
    }
  }

  // Event headers come before the first telescope; the blobs of CORSIKA output
  // have arrived by then, and an event without them gets synthesised ones.
  void ensure_event_headers() {
    if (!event_open_ || event_headers_done_) {
      return;
    }
    event_headers_done_ = true;
    if (!have_evth_) {
      emit(event_header());
    }
    if (!have_teloff_) {
      emit(array_offsets());
    }
  }

  void end_event() {
    if (!event_open_) {
      return;
    }
    ensure_event_headers();
    if (!have_evte_) {
      emit(event_end_block());
    }
    event_open_ = false;
  }

  void end_run() {
    if (!run_open_) {
      return;
    }
    end_event();
    ensure_run_headers();
    if (!have_rune_) {
      emit(run_end_block());
    }
    run_open_ = false;
  }

  void emit(const EventioBlock &block) {
    Record record;
    record.type = RecordType::blob;
    record.blob_eventio_type = block.type;
    record.blob_eventio_id = block.ident;
    record.blob = simtel_detail::block_bytes(block);
    writer_.write(record);
  }

  // RUNH: simtools photon_bunch_file.py write_photon_bunches, run = {2, 3, 4, 5, 6, 17, 18}.
  EventioBlock run_header() const {
    using namespace simtel_detail;
    const auto run = static_cast<std::int64_t>(number_field(run_metadata_, "run_number"));
    const double level_cm = 100.0 * number_field(run_metadata_, "observation_level_m");
    const Words words = {{2, static_cast<double>(run)},
                         {3, date},
                         {4, corsika_version},
                         {5, 1.0},
                         {6, level_cm},
                         {17, energy_gev},
                         {18, energy_gev}};
    return named_block(iact::runh, static_cast<std::uint32_t>(run), "RUNH", words);
  }

  // TELPOS: one entry per telescope id; the position in metres from run metadata.
  EventioBlock telescope_positions() const {
    using namespace simtel_detail;
    TelPos pos = telpos_from_metres(
        array_field(run_metadata_, "telescope_x_m"), array_field(run_metadata_, "telescope_y_m"),
        array_field(run_metadata_, "telescope_z_m"), array_field(run_metadata_, "telescope_r_m"));
    for (const auto &entry : camera_) {
      if (entry.first == 0 || entry.first > pos.z.size()) {
        throw std::runtime_error("camera map telescope " + std::to_string(entry.first) +
                                 " is not in the run metadata");
      }
      pos.z[entry.first - 1] = 0.0;
    }
    return make_telpos_block(pos);
  }

  // EVTH: simtools _EVENT_WORDS plus the event number, observation level, run number,
  // and the profile word (77) that refers to the ATMPROF block, which is not written.
  EventioBlock event_header() const {
    using namespace simtel_detail;
    const double level_cm = 100.0 * number_field(run_metadata_, "observation_level_m");
    const auto run = static_cast<std::int64_t>(number_field(run_metadata_, "run_number"));
    Words words = event_constants;
    words.push_back({2, static_cast<double>(event_number_)});
    words.push_back({5, level_cm});
    words.push_back({7, -level_cm});
    words.push_back({44, static_cast<double>(run)});
    words.push_back({48, level_cm});
    words.push_back({77, static_cast<double>(1024 * profile_id + event_flags)});
    return named_block(iact::evth, static_cast<std::uint32_t>(event_number_), "EVTH", words);
  }

  // TELOFF: one array at zero offset (simtools _telescope_offsets_payload with one reuse).
  EventioBlock array_offsets() const {
    TelOffset off;
    off.toff = 0.0;
    off.xoff = {0.0};
    off.yoff = {0.0};
    return make_teloff_block(off);
  }

  EventioBlock event_end_block() const {
    using namespace simtel_detail;
    return named_block(iact::evte, static_cast<std::uint32_t>(event_number_), "EVTE",
                       {{2, static_cast<double>(event_number_)}});
  }

  // RUNE: "RUNE", run number and number of showers (simtools write_photon_bunches).
  EventioBlock run_end_block() const {
    using namespace simtel_detail;
    const auto run = static_cast<std::int64_t>(number_field(run_metadata_, "run_number"));
    std::vector<std::uint8_t> body;
    ByteWriter out{body};
    out.i32(3);
    for (const char c : std::string("RUNE")) {
      body.push_back(static_cast<std::uint8_t>(c));
    }
    out.f32(static_cast<float>(run));
    out.f32(static_cast<float>(events_));
    return make_block(iact::rune, 0, static_cast<std::uint32_t>(run), std::move(body));
  }

  SimtelWriter writer_;
  CameraMap camera_;
  Metadata run_metadata_;
  std::int64_t event_number_{};
  std::uint64_t events_{};
  bool run_open_{};
  bool run_headers_done_{};
  bool have_runh_{};
  bool have_telpos_{};
  bool have_rune_{};
  bool event_open_{};
  bool event_headers_done_{};
  bool have_evth_{};
  bool have_teloff_{};
  bool have_evte_{};
};

// CORSIKA IACT EventIO (long, compact or split photon blocks) to a photon stream
// with producer `producer`, read from `in` and written to `out`.
inline void corsika_to_stream(std::istream &in, std::ostream &out, const std::string &producer) {
  Writer writer{out, producer};
  corsika_to_records(in, [&writer](const Record &record) { write_record(writer, record); });
  writer.finish();
}

// Photon stream from `in` to sim_telarray EventIO on `out`. Throws on a truncated
// stream or on a record the adapter cannot convert.
inline void stream_to_simtel(std::istream &in, std::ostream &out, CameraMap camera) {
  PsToSimtel converter{out, std::move(camera)};
  Reader reader{in};
  Record record;
  while (reader.next(record)) {
    converter.write(record);
  }
  if (!reader.ended()) {
    throw std::runtime_error("gstrickd: no stream_end (truncated stream)");
  }
  if (in.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("gstrickd: data after stream_end");
  }
  converter.finish();
}

} // namespace gstrickd
