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
// Telescope-frame (T) photons are for dual-mirror telescopes only (camera map
// mirror_class 2), e.g. flashers at the camera firing onto the secondary mirror. Each photon
// is written as a bunch at the plane z = 0 of the T frame (through its origin): the
// crossing point of the photon's line with that plane (x, y), the direction (cx, cy), the
// arrival time t + |s| / c with s the path from the emission point to the plane, and the
// emission altitude zem = observation_level_m + z. The sign of cz is not written. The
// sim_telarray run needs Bypass_Optics=1. The camera rotation is not used.
//
// Telescopes with C photons get z = 0 in TELPOS (1201). sim_telarray moves each
// bunch from the observation level to the telescope by ztel * airmass along its
// direction (sim_telarray.c:2740-2790), which applies to ground-frame bunches
// only. Camera-frame bunches are already at the camera, so the shift must be
// zero (simtools shower_route._zero_traced_telescope_heights, :594).
//
// A telescope is either listed in the camera map (camera-frame or telescope-frame
// photons) or ground-frame (not listed). Mixing ground frame with the others, or camera
// with telescope frame, for one telescope is rejected, because the TELPOS height and the
// Bypass_Optics value are per telescope for the whole run.
//
// Streams without CORSIKA blobs (e.g. from light sources) get the sim_telarray
// headers RUNH, TELPOS, EVTH, TELOFF, EVTE, INPUTCFG (one line, NSHOW <events>, before
// RUNE) and RUNE synthesised from the run and
// event metadata with the words simtools photon_bunch_file.py fills
// (write_photon_bunches, :379-466; _EVENT_WORDS :49-98). Each synthesised header
// is reported as a warning. ATMPROF (1216) is written only when a table is given
// (ps2simtel --atmprof): before the first EVTH, with the profile number that EVTH
// names (99 for synthesised headers). sim_telarray takes the table over the
// profile number only when the block precedes the first event.
//
// Efficiency ledger (docs/format.md): run ledgers apply to every telescope, a
// telescope ledger to its telescope and overrides the run ledgers per factor.
// A factor is accepted when its state is `applied` and sim_telarray does not
// apply it in the telescope's frame, or when its state is `pending` and
// sim_telarray applies it in that frame. `envelope` is an error for factors
// sim_telarray applies, because sim_telarray applies the full factor, not
// f/f_env. With ground-frame photons `atmosphere` is accepted in any state:
// sim_telarray applies its transmission table, which must then be all-zero when
// the atmosphere was applied upstream. A factor listed in `scaled` must be in
// state `envelope` and applied by sim_telarray: its sim_telarray table then
// holds f/f_env (e.g. a QE table divided by the QE envelope). Any other factor is
// an error, unless it is listed in `allow_pending`. Streams without a ledger are
// converted as before.

#include "gstrickd/corsika.hpp"
#include "gstrickd/eventio.hpp"
#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/stream.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <ostream>
#include <set>
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

// Profile number of the custom table (CORSIKA ATMOFILE, simtools _PROFILE_ID).
inline constexpr std::uint32_t custom_profile_id = 99;

// Parses an atmospheric table as CORSIKA and sim_telarray read it (dlufd profile.hpp):
// '#' starts a comment, blank lines are skipped, and the first four columns are altitude
// [km a.s.l.], density [g/cm^3], thickness [g/cm^2] and n-1. Further columns are ignored.
// Altitudes must ascend; 2 to 120 levels (atmo.c MAX_PROFILE). The profile number is
// custom_profile_id; the observation level is set by the caller. Throws
// std::invalid_argument with the line number.
inline AtmosphereTable read_atmprof_text(std::istream &in, const std::string &file_name) {
  AtmosphereTable table;
  table.profile_id = custom_profile_id;
  table.file_name = file_name;
  std::string line;
  std::size_t number = 0;
  const auto fail = [&number](const std::string &what) {
    return std::invalid_argument("atmospheric table line " + std::to_string(number) + ": " + what);
  };
  while (std::getline(in, line)) {
    ++number;
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') {
      continue;
    }
    std::istringstream fields(line);
    std::array<std::string, 4> tokens;
    for (std::string &token : tokens) {
      if (!(fields >> token)) {
        throw fail("expected 4 columns: altitude_km density thickness n_minus_1");
      }
    }
    std::array<double, 4> v{};
    for (std::size_t k = 0; k < v.size(); ++k) {
      if (!simtel_detail::parse_real(tokens[k], v[k])) {
        throw fail("not a finite number: " + tokens[k]);
      }
    }
    if (v[1] < 0.0 || v[2] < 0.0) {
      throw fail("density and thickness must not be negative");
    }
    if (!table.alt_km.empty() && !(v[0] > table.alt_km.back())) {
      throw fail("altitudes must ascend");
    }
    table.alt_km.push_back(v[0]);
    table.rho.push_back(v[1]);
    table.thick.push_back(v[2]);
    table.refidx_m1.push_back(v[3]);
  }
  if (table.alt_km.size() < 2) {
    throw std::invalid_argument("atmospheric table needs at least 2 levels");
  }
  if (table.alt_km.size() > 120) {
    throw std::invalid_argument("atmospheric table has more than 120 levels");
  }
  return table;
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

// Converts telescope-frame (T) emission points and directions in place to bunches at the
// plane z = 0 (see the file comment): position (x, y) of the crossing, direction (cx, cy),
// time t + |s| / c and emission altitude zem = observation_level_m + z, all in metres and
// ns as in the stream. Without a cz column cz = -sqrt(1 - cx^2 - cy^2). Arithmetic in
// double, rounded to float once. Throws std::invalid_argument, leaving the chunk unchanged,
// for a photon with a non-finite value or result, with a direction that is not a unit vector,
// or parallel to the plane.
inline void telescope_to_bypass_convention(PhotonChunk &chunk, double observation_level_m) {
  if (!(chunk.has(Column::x) && chunk.has(Column::y) && chunk.has(Column::z) &&
        chunk.has(Column::cx) && chunk.has(Column::cy) && chunk.has(Column::t))) {
    throw std::invalid_argument("telescope-frame chunk lacks x, y, z, cx, cy or t");
  }
  if (!std::isfinite(observation_level_m)) {
    throw std::invalid_argument("observation level must be finite");
  }
  constexpr double speed_of_light_m_per_ns = 0.299792458;
  const std::size_t n = chunk.rows();
  const std::vector<float> &x = chunk[Column::x];
  const std::vector<float> &y = chunk[Column::y];
  const std::vector<float> &t = chunk[Column::t];
  const std::vector<float> &z = chunk[Column::z];
  const std::vector<float> &cx = chunk[Column::cx];
  const std::vector<float> &cy = chunk[Column::cy];
  const bool has_cz = chunk.has(Column::cz);
  // Direction cosines are unit vectors to within the float32 rounding of the stream.
  constexpr double unit_tolerance = 1e-4;
  // The results are collected first, so a rejected chunk is left unchanged.
  std::vector<float> new_x(n);
  std::vector<float> new_y(n);
  std::vector<float> new_t(n);
  std::vector<float> zem(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double cxi = cx[i];
    const double cyi = cy[i];
    const double norm2 = cxi * cxi + cyi * cyi;
    const double czi =
        has_cz ? static_cast<double>(chunk[Column::cz][i]) : -std::sqrt(std::max(0.0, 1.0 - norm2));
    if (!std::isfinite(cxi) || !std::isfinite(cyi) || !std::isfinite(czi) || !std::isfinite(x[i]) ||
        !std::isfinite(y[i]) || !std::isfinite(z[i]) || !std::isfinite(t[i])) {
      throw std::invalid_argument("telescope-frame photon with a non-finite value");
    }
    if (std::abs(norm2 + czi * czi - 1.0) > unit_tolerance) {
      throw std::invalid_argument("telescope-frame photon direction (cx, cy, cz) is not a unit "
                                  "vector");
    }
    if (czi == 0.0) {
      throw std::invalid_argument("telescope-frame photon parallel to the plane z = 0");
    }
    const double path = -static_cast<double>(z[i]) / czi;
    new_x[i] = static_cast<float>(x[i] + path * cxi);
    new_y[i] = static_cast<float>(y[i] + path * cyi);
    new_t[i] = static_cast<float>(t[i] + std::abs(path) / speed_of_light_m_per_ns);
    zem[i] = static_cast<float>(observation_level_m + z[i]);
    if (!std::isfinite(new_x[i]) || !std::isfinite(new_y[i]) || !std::isfinite(new_t[i]) ||
        !std::isfinite(zem[i])) {
      throw std::invalid_argument("telescope-frame photon converts to a value out of range");
    }
  }
  chunk[Column::x] = std::move(new_x);
  chunk[Column::y] = std::move(new_y);
  chunk[Column::t] = std::move(new_t);
  chunk[Column::zem] = std::move(zem);
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

// Efficiency ledger rules for ps2simtel (docs/usage.md, "Efficiency ledger").
// sim_telarray applies these factors itself, by the frame of the telescope.
// Camera frame: the camera efficiencies. Ground frame: also mirror reflectivity
// and telescope transmission. Atmosphere is not in either list.
inline bool sim_telarray_applies(const std::string &factor, bool camera_frame) {
  static const std::set<std::string> camera_factors{"qe", "camera_transmission", "camera_filter"};
  static const std::set<std::string> ground_factors{"qe", "camera_transmission", "camera_filter",
                                                    "mirror_reflectivity",
                                                    "telescope_transmission"};
  return (camera_frame ? camera_factors : ground_factors).contains(factor);
}

// Factor names of the ledger vocabulary (docs/format.md). Any other name is unknown.
inline bool known_factor(const std::string &factor) {
  return sim_telarray_applies(factor, false) || factor == "atmosphere" ||
         factor == "telescope_transmission_constant";
}

// Factor -> state of one ledger record.
using LedgerStates = std::map<std::string, std::string>;

inline LedgerStates ledger_states(const Metadata &ledger) {
  const auto factors = ledger.find("factors");
  if (factors == ledger.end()) {
    throw std::runtime_error("efficiency ledger lacks 'factors'");
  }
  const auto *list = std::get_if<std::string>(&factors->second);
  if (list == nullptr) {
    throw std::runtime_error("efficiency ledger 'factors' is not a string");
  }
  LedgerStates states;
  std::istringstream in(*list);
  std::string name;
  while (std::getline(in, name, ',')) {
    const auto first = name.find_first_not_of(" \t");
    if (first == std::string::npos) {
      continue;
    }
    name = name.substr(first, name.find_last_not_of(" \t") - first + 1);
    const auto state_it = ledger.find(name + ".state");
    if (state_it == ledger.end()) {
      throw std::runtime_error("efficiency ledger lacks '" + name + ".state'");
    }
    const auto *state = std::get_if<std::string>(&state_it->second);
    if (state == nullptr || (*state != "applied" && *state != "envelope" && *state != "pending")) {
      throw std::runtime_error("efficiency ledger factor '" + name +
                               "' has an unknown state (applied, envelope or pending)");
    }
    if (!states.emplace(name, *state).second) {
      throw std::runtime_error("efficiency ledger lists factor '" + name + "' twice");
    }
  }
  return states;
}

// Checks the ledger of one telescope: run ledgers, then its telescope ledgers,
// which override the run ledgers per factor. Throws std::runtime_error naming the
// first factor that sim_telarray would apply twice, or would not apply at all.
inline void check_ledger(const std::vector<Metadata> &run_ledgers,
                         const std::vector<Metadata> &telescope_ledgers, std::uint32_t telescope_id,
                         bool camera_frame, const std::set<std::string> &allow_pending,
                         const std::set<std::string> &scaled = {}) {
  LedgerStates states;
  for (const Metadata &ledger : run_ledgers) {
    for (const auto &[name, state] : ledger_states(ledger)) {
      if (!states.emplace(name, state).second) {
        throw std::runtime_error("efficiency ledger lists factor '" + name +
                                 "' in more than one run ledger");
      }
    }
  }
  LedgerStates own;
  for (const Metadata &ledger : telescope_ledgers) {
    for (const auto &[name, state] : ledger_states(ledger)) {
      if (!own.emplace(name, state).second) {
        throw std::runtime_error("efficiency ledger lists factor '" + name +
                                 "' twice for telescope " + std::to_string(telescope_id));
      }
      states[name] = state;
    }
  }
  const std::string where = "telescope " + std::to_string(telescope_id) +
                            (camera_frame ? " (camera frame)" : " (ground frame)");
  for (const auto &[name, state] : states) {
    // Scaled factors are checked first, so no shortcut below accepts them in another state.
    if (scaled.contains(name)) {
      if (!sim_telarray_applies(name, camera_frame) || state != "envelope") {
        throw std::runtime_error("efficiency ledger of " + where + ": factor '" + name +
                                 "' is given as scaled (sim_telarray applies f/f_env), so it "
                                 "must be in state envelope and applied by sim_telarray in "
                                 "this frame; it is " +
                                 state);
      }
      continue;
    }
    if (allow_pending.contains(name) || (name == "atmosphere" && !camera_frame)) {
      continue;
    }
    const std::string option = "; pass --allow-pending " + name + " to accept it as it is";
    if (!known_factor(name)) {
      throw std::runtime_error("efficiency ledger of " + where + ": unknown factor '" + name + "'" +
                               option);
    }
    const bool sim_telarray = sim_telarray_applies(name, camera_frame);
    if (state == "applied") {
      if (sim_telarray) {
        throw std::runtime_error("efficiency ledger of " + where + ": factor '" + name +
                                 "' is already applied upstream, and sim_telarray applies it "
                                 "again in this frame" +
                                 option);
      }
      continue;
    }
    if (sim_telarray && state == "envelope") {
      throw std::runtime_error("efficiency ledger of " + where + ": factor '" + name +
                               "' has only its envelope applied, but sim_telarray applies the "
                               "full factor, so the result would be too dim; pass "
                               "--scaled-factor " +
                               name + " if its sim_telarray table holds f/f_env" + option);
    }
    if (!sim_telarray) {
      throw std::runtime_error("efficiency ledger of " + where + ": factor '" + name + "' is " +
                               state +
                               ", but sim_telarray does not apply it in this frame; "
                               "apply it upstream (state applied)" +
                               option);
    }
  }
}

} // namespace simtel_detail

// Options of the photon stream to sim_telarray conversion.
struct SimtelOptions {
  CameraMap camera;                       // camera-frame telescopes
  std::optional<AtmosphereTable> atmprof; // atmospheric table (ATMPROF block), see the file comment
  std::ostream *warnings{nullptr};        // reports synthesised headers; null: silent
  std::set<std::string> allow_pending;    // ledger factors accepted as they are
  std::set<std::string> scaled;           // ledger factors whose camera table holds f/f_env
};

// Writes a photon stream as sim_telarray input (EventIO, normal mode: one 1204
// container per array). Records are read in stream order; see the file comment
// for the frame and header rules.
class PsToSimtel {
public:
  PsToSimtel(std::ostream &out, CameraMap camera)
      : PsToSimtel(out, SimtelOptions{std::move(camera), std::nullopt, nullptr, {}, {}}) {}

  PsToSimtel(std::ostream &out, SimtelOptions options)
      : writer_{out}, camera_{std::move(options.camera)}, atmprof_{std::move(options.atmprof)},
        warnings_{options.warnings}, allow_pending_{std::move(options.allow_pending)},
        scaled_{std::move(options.scaled)} {
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
      return; // not part of sim_telarray input
    case RecordType::ledger:
      ledger(record.metadata);
      return;
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
      begin_telescope(record.metadata);
      writer_.write(record);
      return;
    case RecordType::photons:
      if (telescope_open_) {
        telescope_photons_ = true;
      }
      photons(record.chunk);
      return;
    case RecordType::telescope_end:
      end_telescope();
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
      throw std::runtime_error(
          "photons in the sky frame cannot go to sim_telarray: convert them to the ground, "
          "telescope or camera frame first");
    case Frame::ground:
      if (camera_.contains(chunk.telescope_id)) {
        throw std::runtime_error("telescope " + std::to_string(chunk.telescope_id) +
                                 " is in the camera map but has ground-frame photons");
      }
      writer_.write_chunk(chunk);
      return;
    case Frame::telescope: {
      const auto it = camera_.find(chunk.telescope_id);
      if (it == camera_.end()) {
        throw std::runtime_error("telescope " + std::to_string(chunk.telescope_id) +
                                 " has telescope-frame photons but no camera map entry");
      }
      if (it->second.mirror_class != 2) {
        throw std::runtime_error("telescope " + std::to_string(chunk.telescope_id) +
                                 " has telescope-frame photons but is not a dual-mirror "
                                 "telescope (mirror_class 2 in the camera map)");
      }
      note_frame(chunk, it->second);
      PhotonChunk converted = chunk;
      telescope_to_bypass_convention(
          converted, simtel_detail::number_field(run_metadata_, "observation_level_m"));
      converted.frame = Frame::ground;
      writer_.write_chunk(converted);
      return;
    }
    case Frame::camera: {
      const auto it = camera_.find(chunk.telescope_id);
      if (it == camera_.end()) {
        throw std::runtime_error("telescope " + std::to_string(chunk.telescope_id) +
                                 " has camera-frame photons but no camera map entry");
      }
      note_frame(chunk, it->second);
      PhotonChunk converted = chunk;
      camera_to_bypass_convention(converted, it->second);
      converted.frame = Frame::ground; // sim_telarray input is written as ground-frame bunches
      writer_.write_chunk(converted);
      return;
    }
    }
  }

  // Remembers the frame of a camera-map telescope, rejects a second frame, and reports the
  // Bypass_Optics value the first time.
  void note_frame(const PhotonChunk &chunk, const CameraTelescope &camera) {
    const auto [it, inserted] = frames_.emplace(chunk.telescope_id, chunk.frame);
    if (!inserted) {
      if (it->second != chunk.frame) {
        throw std::runtime_error("telescope " + std::to_string(chunk.telescope_id) +
                                 " has photons in the camera and the telescope frame");
      }
      return;
    }
    if (warnings_ == nullptr) {
      return;
    }
    if (chunk.frame == Frame::telescope) {
      *warnings_ << "ps2simtel: telescope " << chunk.telescope_id
                 << ": telescope frame, mirror class 2, Bypass_Optics=1\n";
    } else {
      *warnings_ << "ps2simtel: telescope " << chunk.telescope_id << ": camera frame, mirror class "
                 << camera.mirror_class
                 << ", Bypass_Optics=" << bypass_optics_level(camera.mirror_class)
                 << ", camera rotation " << camera.rotation_deg << " deg\n";
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
    case iact::atmprof:
      if (atmprof_) {
        throw std::runtime_error("the stream has its own ATMPROF block; drop --atmprof");
      }
      break;
    case iact::evth:
      have_evth_ = true;
      if (atmprof_pending_) {
        emit_atmprof(evth_profile(record.blob));
      }
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
    case iact::inputcfg:
      have_inputcfg_ = true;
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
    atmprof_pending_ = atmprof_.has_value();
    have_runh_ = have_telpos_ = have_rune_ = have_inputcfg_ = false;
    events_ = 0;
    run_ledgers_.clear();
    telescope_ledgers_.clear();
    telescope_open_ = false;
    frames_.clear();
  }

  // Ledger of the run (before the first event) or of the open telescope (before its photons).
  void ledger(const Metadata &metadata) {
    if (!run_open_) {
      throw std::runtime_error("ledger record outside a run");
    }
    if (telescope_open_) {
      if (telescope_photons_) {
        throw std::runtime_error("efficiency ledger for telescope " +
                                 std::to_string(telescope_id_) + " comes after its photons");
      }
      telescope_ledgers_.push_back(metadata);
      return;
    }
    if (events_ > 0) {
      throw std::runtime_error("run efficiency ledger comes after the first event");
    }
    run_ledgers_.push_back(metadata);
  }

  void begin_telescope(const Metadata &metadata) {
    telescope_open_ = true;
    telescope_photons_ = false;
    telescope_ledgers_.clear();
    const auto it = metadata.find("telescope_id");
    telescope_id_known_ = it != metadata.end() && std::holds_alternative<std::int64_t>(it->second);
    telescope_id_ =
        telescope_id_known_ ? static_cast<std::uint32_t>(std::get<std::int64_t>(it->second)) : 0;
  }

  // Checks the efficiency ledger of the telescope that ends here (no-op without ledgers).
  void end_telescope() {
    if (!telescope_open_) {
      return;
    }
    telescope_open_ = false;
    if (run_ledgers_.empty() && telescope_ledgers_.empty()) {
      return;
    }
    if (!telescope_id_known_ || telescope_id_ == 0) {
      throw std::runtime_error("telescope record lacks telescope_id, needed for the efficiency "
                               "ledger");
    }
    simtel_detail::check_ledger(run_ledgers_, telescope_ledgers_, telescope_id_,
                                camera_.contains(telescope_id_), allow_pending_, scaled_);
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
      note_synthesised("RUNH");
      emit(run_header());
    }
    if (!have_telpos_) {
      note_synthesised("TELPOS");
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
      // Without a table the EVTH profile is not in the stream: say so.
      note_synthesised("EVTH");
      if (!atmprof_) {
        note_missing_atmprof();
      }
      emit_atmprof(custom_profile_id);
      emit(event_header());
    }
    if (!have_teloff_) {
      note_synthesised("TELOFF");
      emit(array_offsets());
    }
  }

  void end_event() {
    if (!event_open_) {
      return;
    }
    end_telescope();
    ensure_event_headers();
    if (!have_evte_) {
      note_synthesised("EVTE");
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
      if (!have_inputcfg_) {
        // The receiving run checks the number of showers it reads against the NSHOW
        // line of the CORSIKA input, which a stream without CORSIKA blobs lacks.
        note_synthesised("INPUTCFG");
        emit(input_lines_block());
      }
      note_synthesised("RUNE");
      emit(run_end_block());
    }
    run_open_ = false;
  }

  void note_synthesised(const char *name) {
    if (warnings_ == nullptr || !warned_.insert(name).second) {
      return;
    }
    *warnings_ << "ps2simtel: warning: no CORSIKA " << name
               << " block in the stream; synthesised it from the run metadata\n";
  }

  void note_missing_atmprof() {
    if (warnings_ == nullptr || !warned_.insert("ATMPROF").second) {
      return;
    }
    *warnings_ << "ps2simtel: warning: no atmospheric table; sim_telarray cannot find profile "
               << custom_profile_id << " (use --atmprof FILE)\n";
  }

  // Profile number of an EVTH block (CORSIKA word 77, bits 10-19).
  static std::uint32_t evth_profile(const std::vector<std::uint8_t> &raw) {
    const std::vector<float> words = decode_float_block(simtel_detail::parse_block(raw));
    if (words.size() < corsika_word::evth_flags + 1) {
      throw std::runtime_error("EVTH block too short for its profile word");
    }
    const auto flags = static_cast<std::int64_t>(std::llround(words[corsika_word::evth_flags]));
    return static_cast<std::uint32_t>((flags >> corsika_word::evth_flag_atmosphere_shift) &
                                      corsika_word::evth_flag_atmosphere_mask);
  }

  // Writes the ATMPROF block once per run, when a table was given, with `profile`
  // as its number and the run's observation level.
  void emit_atmprof(std::uint32_t profile) {
    if (!atmprof_pending_) {
      return;
    }
    atmprof_pending_ = false;
    AtmosphereTable table = *atmprof_;
    table.profile_id = profile;
    table.observation_level_cm =
        100.0 * simtel_detail::number_field(run_metadata_, "observation_level_m");
    emit(make_atmprof_block(table));
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

  // INPUTCFG: the single input line `NSHOW n` with the number of events written.
  EventioBlock input_lines_block() const {
    std::vector<std::uint8_t> body;
    ByteWriter out{body};
    const std::string line = "NSHOW " + std::to_string(events_);
    out.i32(1);
    out.i16(static_cast<std::int16_t>(line.size()));
    body.insert(body.end(), line.begin(), line.end());
    return make_block(iact::inputcfg, 0, 0, std::move(body));
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
  std::map<std::uint32_t, Frame> frames_; // frame of each camera-map telescope, per run
  std::optional<AtmosphereTable> atmprof_;
  std::ostream *warnings_{};
  std::set<std::string> warned_;
  bool atmprof_pending_{};
  Metadata run_metadata_;
  std::int64_t event_number_{};
  std::uint64_t events_{};
  bool run_open_{};
  bool run_headers_done_{};
  bool have_runh_{};
  bool have_telpos_{};
  bool have_rune_{};
  bool have_inputcfg_{};
  bool event_open_{};
  bool event_headers_done_{};
  bool have_evth_{};
  bool have_teloff_{};
  bool have_evte_{};
  std::set<std::string> allow_pending_;
  std::set<std::string> scaled_;
  std::vector<Metadata> run_ledgers_;
  std::vector<Metadata> telescope_ledgers_;
  bool telescope_open_{};
  bool telescope_photons_{};
  bool telescope_id_known_{};
  std::uint32_t telescope_id_{};
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
inline void stream_to_simtel(std::istream &in, std::ostream &out, SimtelOptions options) {
  PsToSimtel converter{out, std::move(options)};
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

inline void stream_to_simtel(std::istream &in, std::ostream &out, CameraMap camera) {
  stream_to_simtel(in, out, SimtelOptions{std::move(camera), std::nullopt, nullptr, {}, {}});
}

} // namespace gstrickd
