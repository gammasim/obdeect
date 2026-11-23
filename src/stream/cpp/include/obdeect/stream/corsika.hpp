#pragma once

// CORSIKA IACT output (EventIO 1200-1217) to photon-stream records, and back to
// sim_telarray input. Uses only eventio.hpp, telpos.hpp, constants.hpp and the stream types; the
// stream Writer is not needed.

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/eventio.hpp"
#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/telpos.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <istream>
#include <map>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace obdeect::stream {

// Word indices into the CORSIKA RUNH (1200) and EVTH (1202) float arrays: C
// index = Fortran word - 1 (e.g. EVTH word 3 is the primary particle, word 7 the
// first interaction height).
namespace corsika_word {
// RUNH
inline constexpr std::size_t runh_run_number = 1;  // word 2
inline constexpr std::size_t runh_n_levels = 4;    // word 5
inline constexpr std::size_t runh_level_1 = 5;     // word 6: observation level 1 [cm]
inline constexpr std::size_t runh_slope = 15;      // word 16
inline constexpr std::size_t runh_energy_min = 16; // word 17 [GeV]
inline constexpr std::size_t runh_energy_max = 17; // word 18 [GeV]
// EVTH
inline constexpr std::size_t evth_event_number = 1; // word 2
inline constexpr std::size_t evth_particle = 2;     // word 3: particle code
inline constexpr std::size_t evth_energy = 3;       // word 4 [GeV]
inline constexpr std::size_t evth_first_height = 6; // word 7 [cm]
inline constexpr std::size_t evth_zenith = 10;      // word 11 [rad]
inline constexpr std::size_t evth_azimuth = 11;     // word 12 [rad, CORSIKA convention]
inline constexpr std::size_t evth_flags = 76;       // word 77: integer flags (see below)
// Word 79: CURVED flag (0 = standard, 1 = curved, 2 = curved with upgoing primaries).
inline constexpr std::size_t evth_curved_flag = 78;
inline constexpr std::size_t evth_wavelength_min = 95; // word 96 [nm]
inline constexpr std::size_t evth_wavelength_max = 96; // word 97 [nm]
// Word 77 as an integer (corsika.tex:8624-8637, bit n = 1 << (n - 1)): bit 5 =
// atmospheric refraction enabled; bits 11-20 = external atmosphere table.
inline constexpr int evth_flag_refraction_bit = 4;
inline constexpr int evth_flag_atmosphere_shift = 10;
inline constexpr std::int64_t evth_flag_atmosphere_mask = 1023;
} // namespace corsika_word

// Lengths: CORSIKA cm -> stream m. The division is done in double and rounded
// to float once. A round trip back to cm is not always exact (see test_corsika).
[[nodiscard]] inline float corsika_cm_to_m(float value_cm) {
  return static_cast<float>(static_cast<double>(value_cm) / 100.0);
}
[[nodiscard]] inline float corsika_m_to_cm(float value_m) {
  return static_cast<float>(static_cast<double>(value_m) * 100.0);
}

namespace detail {

inline float header_word(const std::vector<float> &words, std::size_t index, const char *block) {
  if (index >= words.size()) {
    throw std::runtime_error(std::string("CORSIKA ") + block + " block too short");
  }
  return words[index];
}

inline std::int64_t round_to_int(double value) {
  return static_cast<std::int64_t>(std::llround(value));
}

inline std::int64_t metadata_int(const Metadata &metadata, const std::string &key) {
  const auto it = metadata.find(key);
  if (it == metadata.end()) {
    throw std::runtime_error("record lacks metadata key " + key);
  }
  const auto *value = std::get_if<std::int64_t>(&it->second);
  if (value == nullptr) {
    throw std::runtime_error("metadata key " + key + " is not an int");
  }
  return *value;
}

// Geometry of one 3D photon bunch (1217) of a sphere run, in m and ns. The bunch gives the arrival
// point A (x, y, z, relative to the telescope sphere centre) and the direction d = (cx, cy, cz);
// the emission point is E = A - dist d. The G sphere row (docs/format.md) is E itself if |E| <= R
// (path 0), else the entry point E + s d with the smallest s >= 0 where |E + s d| = R (path s). A
// bunch without such s misses the sphere. The time at the row is ctime + (s - dist) / c, with the
// vacuum speed of light: ctime refers to the arrival point, at s = dist.
struct SphereBunch {
  bool hits{};
  double x{}, y{}, z{}; // row [m], relative to the sphere centre
  double path{};        // distance from the emission point to the row [m]
  double time_ns{};     // arrival time at the row [ns]
  double zem_m{};       // emission altitude a.s.l. [m]
};

// `radius_m` is the sphere radius of the telescope, `altitude_of_centre_m` the altitude a.s.l. of
// its centre (observation level + TELPOS z), so that zem = altitude + E_z (flat geometry). Throws
// std::runtime_error for non-finite values, a negative distance or photons, and a direction that is
// not a unit vector (cz_is_consistent).
inline SphereBunch sphere_bunch(const Bunch3d &b, double radius_m, double altitude_of_centre_m) {
  const double finite_values[] = {b.x,  b.y,     b.z,    b.cx,      b.cy,
                                  b.cz, b.ctime, b.dist, b.photons, b.lambda};
  for (const double value : finite_values) {
    if (!std::isfinite(value)) {
      throw std::runtime_error("3D photon bunch with a non-finite value");
    }
  }
  if (!cz_is_consistent(b.cz, b.cx, b.cy)) {
    throw std::runtime_error("3D photon bunch direction (cx, cy, cz) is not a unit vector");
  }
  if (b.dist < 0.0F) {
    throw std::runtime_error("3D photon bunch with a negative distance");
  }
  if (b.photons < 0.0F) {
    throw std::runtime_error("3D photon bunch with negative photons (not supported)");
  }

  const double d[3] = {static_cast<double>(b.cx), static_cast<double>(b.cy),
                       static_cast<double>(b.cz)};
  const double dist_m = cm_to_m(static_cast<double>(b.dist));
  const double a[3] = {cm_to_m(static_cast<double>(b.x)), cm_to_m(static_cast<double>(b.y)),
                       cm_to_m(static_cast<double>(b.z))};
  const double e[3] = {a[0] - dist_m * d[0], a[1] - dist_m * d[1], a[2] - dist_m * d[2]};
  const double e2 = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
  const double r2 = radius_m * radius_m;

  SphereBunch out;
  double s = 0.0;
  out.hits = true;
  if (e2 > r2) {
    // |E + s d|^2 = R^2 with |d| = 1 (within the direction tolerance): s^2 + 2 (E.d) s + (|E|^2 -
    // R^2) = 0. Both roots have the same sign (product > 0), positive iff E.d < 0.
    const double edotd = e[0] * d[0] + e[1] * d[1] + e[2] * d[2];
    const double disc = edotd * edotd - (e2 - r2);
    if (edotd >= 0.0 || disc < 0.0) {
      out.hits = false;
      return out;
    }
    s = -edotd - std::sqrt(disc);
  }
  out.x = e[0] + s * d[0];
  out.y = e[1] + s * d[1];
  out.z = e[2] + s * d[2];
  out.path = s;
  out.time_ns = static_cast<double>(b.ctime) + (s - dist_m) / speed_of_light_m_per_ns;
  out.zem_m = altitude_of_centre_m + e[2];
  return out;
}

// One G sphere row of a 1217 bunch as stored in a photon chunk (m, ns, float32).
struct SphereRowValues {
  float x, y, z, cx, cy, cz, t, zem, w, path, wavelength;
};

// Converts one CORSIKA IACT run (blocks in file order) into stream records.
// Order of records follows docs/format.md: run, then per event: event, blobs
// and telescopes in file order, event_end; then run_end. The run record is
// written at the first event header, so the blocks before it (RUNH, INPUTCFG,
// TELPOS, ...) are held back and written after it.
class CorsikaToRecords {
public:
  CorsikaToRecords(std::istream &in, const std::function<void(const Record &)> &sink)
      : reader_{in}, sink_{sink} {}

  void run() {
    while (std::optional<EventioBlock> block = pop()) {
      switch (block->type) {
      case iact::runh:
        close_run();
        rune_seen_ = false;
        runh_ = decode_float_block(*block);
        blob(std::move(*block));
        break;
      case iact::telpos:
        telpos_ = decode_telpos(*block);
        blob(std::move(*block));
        break;
      case iact::evth:
        begin_event(std::move(*block));
        break;
      case iact::teloff:
        blob(std::move(*block));
        break;
      case iact::photons: {
        TelPhotons photons = decode_photons(*block);
        if (photons.array == 999 && photons.tel == 999) {
          blob(std::move(*block)); // particle block: not photons of a telescope
        } else {
          require_event("photon block");
          emit_photons(photons);
        }
        break;
      }
      case iact::telarray:
        require_event("array container");
        for (const EventioBlock &child : subblocks(*block)) {
          if (child.type == iact::photons) {
            emit_photons(decode_photons(child));
          } else if (child.type == iact::photons3d) {
            const TelPhotons3d photons = decode_photons3d(child);
            if (photons.array == 999 && photons.tel == 999) {
              throw std::runtime_error("particle block (999, 999) inside an array container");
            }
            emit_photons3d(photons);
          } else {
            throw std::runtime_error("unexpected block inside an array container");
          }
        }
        break;
      case iact::photons3d: {
        TelPhotons3d photons = decode_photons3d(*block);
        if (photons.array == 999 && photons.tel == 999) {
          blob(std::move(*block)); // particle block: not photons of a telescope
        } else {
          require_event("3D photon block");
          emit_photons3d(photons);
        }
        break;
      }
      case iact::telarray_head: // split-mode markers: kept as blobs, in input order
      case iact::telarray_end:
        blob(std::move(*block));
        break;
      case iact::evte:
        blob(std::move(*block));
        close_event();
        break;
      case iact::rune:
        blob(std::move(*block));
        rune_seen_ = true;
        close_run();
        break;
      default:
        blob(std::move(*block));
        break;
      }
    }
    check_complete();
  }

private:
  std::optional<EventioBlock> pop() {
    if (!queue_.empty()) {
      EventioBlock block = std::move(queue_.front());
      queue_.pop_front();
      return block;
    }
    return read_one();
  }

  std::optional<EventioBlock> read_one() {
    EventioBlock block;
    if (!reader_.next(block)) {
      return std::nullopt;
    }
    return block;
  }

  // Records of a run are held back until its ground representation is known (the first photon
  // block, or the end of the run). Only the header blocks of the first event are held.
  void out(Record &&record) {
    if (representation_) {
      sink_(record);
    } else {
      held_.push_back(std::move(record));
    }
  }

  void emit(RecordType type, Metadata metadata = {}) {
    Record record;
    record.type = type;
    record.metadata = std::move(metadata);
    out(std::move(record));
  }

  void emit_blob(const EventioBlock &block) {
    Record record;
    record.type = RecordType::blob;
    record.blob_eventio_type = block.type;
    record.blob_eventio_id = block.ident;
    record.blob = block.raw();
    out(std::move(record));
  }

  // The first photon block of the run fixes its representation; a block of the other kind
  // afterwards is an error (one representation per run).
  void claim(GroundRepresentation representation) {
    if (!representation_) {
      decide(representation);
    } else if (*representation_ != representation) {
      throw std::runtime_error("CORSIKA run mixes 1205 and 1217 photon blocks; a run has one "
                               "ground representation");
    }
  }

  // Sets the run key and writes the held records of the run.
  void decide(GroundRepresentation representation) {
    representation_ = representation;
    // Plane runs keep their records as before (an absent key means plane).
    if (representation == GroundRepresentation::sphere && !held_.empty() &&
        held_.front().type == RecordType::run) {
      held_.front().metadata["ground_representation"] = std::string("sphere");
    }
    for (const Record &record : held_) {
      sink_(record);
    }
    held_.clear();
  }

  // Blocks outside an event (run header and friends) wait for the run record.
  void blob(EventioBlock &&block) {
    if (run_open_) {
      emit_blob(block);
    } else {
      pending_.push_back(std::move(block));
    }
  }

  void require_event(const char *what) const {
    if (!event_open_) {
      throw std::runtime_error(std::string(what) + " outside an event");
    }
  }

  // Look ahead for TELOFF (written after EVTH). Blocks read here go back into
  // the queue, so the main loop still sees them in file order.
  std::optional<TelOffset> find_teloff() {
    for (std::size_t i = 0;; ++i) {
      if (i == queue_.size()) {
        std::optional<EventioBlock> block = read_one();
        if (!block) {
          return std::nullopt;
        }
        queue_.push_back(std::move(*block));
      }
      const EventioBlock &next = queue_[i];
      if (next.type == iact::teloff) {
        return decode_teloff(next);
      }
      switch (next.type) {
      case iact::evth:
      case iact::runh:
      case iact::telarray:
      case iact::photons:
      case iact::photons3d:
      case iact::evte:
      case iact::rune:
      case iact::telarray_head:
        return std::nullopt;
      default:
        break;
      }
    }
  }

  void begin_event(EventioBlock &&block) {
    const std::vector<float> evth = decode_float_block(block);
    close_event();
    if (!run_open_) {
      open_run(&evth);
    }
    const std::optional<TelOffset> teloff = find_teloff();

    Metadata metadata;
    metadata["event_number"] =
        round_to_int(header_word(evth, corsika_word::evth_event_number, "EVTH"));
    metadata["primary_id"] = round_to_int(header_word(evth, corsika_word::evth_particle, "EVTH"));
    metadata["energy_gev"] =
        static_cast<double>(header_word(evth, corsika_word::evth_energy, "EVTH"));
    metadata["zenith_rad"] =
        static_cast<double>(header_word(evth, corsika_word::evth_zenith, "EVTH"));
    metadata["azimuth_rad"] =
        static_cast<double>(header_word(evth, corsika_word::evth_azimuth, "EVTH"));
    metadata["first_interaction_m"] =
        static_cast<double>(header_word(evth, corsika_word::evth_first_height, "EVTH")) / 100.0;
    if (teloff) {
      std::vector<double> core_x;
      std::vector<double> core_y;
      for (std::size_t i = 0; i < teloff->xoff.size(); ++i) {
        core_x.push_back(teloff->xoff[i] / 100.0);
        core_y.push_back(teloff->yoff[i] / 100.0);
      }
      metadata["core_x_m"] = std::move(core_x);
      metadata["core_y_m"] = std::move(core_y);
      metadata["time_offset_ns"] = teloff->toff;
    }
    emit(RecordType::event, std::move(metadata));
    event_open_ = true;
    next_photon_index_.clear();
    blob(std::move(block));
  }

  void close_event() {
    if (event_open_) {
      emit(RecordType::event_end);
      event_open_ = false;
    }
  }

  // EOF is valid only after the explicit end markers. An open event (no EVTE),
  // an open run (no RUNE) or run blocks without a RUNE mean a truncated file.
  void check_complete() const {
    if (event_open_) {
      throw std::runtime_error("CORSIKA input ends inside an event (no EVTE block)");
    }
    if (run_open_) {
      throw std::runtime_error("CORSIKA input ends inside a run (no RUNE block)");
    }
    if (!pending_.empty() && !rune_seen_) {
      throw std::runtime_error("CORSIKA input ends without a RUNE block");
    }
  }

  Metadata run_metadata(const std::vector<float> *evth) const {
    Metadata metadata;
    if (runh_) {
      const std::vector<float> &runh = *runh_;
      metadata["run_number"] =
          round_to_int(header_word(runh, corsika_word::runh_run_number, "RUNH"));
      metadata["observation_level_m"] =
          static_cast<double>(header_word(runh, corsika_word::runh_level_1, "RUNH")) / 100.0;
      metadata["spectral_index"] =
          static_cast<double>(header_word(runh, corsika_word::runh_slope, "RUNH"));
      metadata["energy_min_gev"] =
          static_cast<double>(header_word(runh, corsika_word::runh_energy_min, "RUNH"));
      metadata["energy_max_gev"] =
          static_cast<double>(header_word(runh, corsika_word::runh_energy_max, "RUNH"));
    }
    if (telpos_) {
      const TelPos &pos = *telpos_;
      std::vector<double> x, y, z, r;
      for (std::size_t i = 0; i < pos.x.size(); ++i) {
        x.push_back(cm_to_m(pos.x[i]));
        y.push_back(cm_to_m(pos.y[i]));
        z.push_back(cm_to_m(pos.z[i]));
        r.push_back(cm_to_m(pos.r[i]));
      }
      metadata["n_telescopes"] = static_cast<std::int64_t>(pos.x.size());
      metadata["telescope_x_m"] = std::move(x);
      metadata["telescope_y_m"] = std::move(y);
      metadata["telescope_z_m"] = std::move(z);
      metadata["telescope_r_m"] = std::move(r);
    }
    if (evth != nullptr) {
      const auto flags = round_to_int(header_word(*evth, corsika_word::evth_flags, "EVTH"));
      metadata["atmosphere_id"] = (flags >> corsika_word::evth_flag_atmosphere_shift) &
                                  corsika_word::evth_flag_atmosphere_mask;
      metadata["refraction_applied"] = (flags >> corsika_word::evth_flag_refraction_bit) & 1;
      // CURVED: 1 = downgoing primaries, 2 = upgoing primaries also allowed.
      const int curved = round_to_int(header_word(*evth, corsika_word::evth_curved_flag, "EVTH"));
      metadata["curvature_applied"] =
          (curved == 1 || curved == 2) ? std::int64_t{1} : std::int64_t{0};
      metadata["wavelength_min_nm"] =
          static_cast<double>(header_word(*evth, corsika_word::evth_wavelength_min, "EVTH"));
      metadata["wavelength_max_nm"] =
          static_cast<double>(header_word(*evth, corsika_word::evth_wavelength_max, "EVTH"));
    }
    return metadata;
  }

  // Writes the run record and the blocks held back for it.
  void open_run(const std::vector<float> *evth) {
    representation_.reset();
    misses_sphere_ = 0;
    emit(RecordType::run, run_metadata(evth));
    run_open_ = true;
    for (const EventioBlock &block : pending_) {
      emit_blob(block);
    }
    pending_.clear();
  }

  // A run without photon blocks is a plane run. Sphere runs count the bunches that miss the
  // telescope sphere (misses_sphere in the run end).
  void close_run() {
    close_event();
    if (!run_open_ && pending_.empty()) {
      runh_.reset();
      telpos_.reset();
      return;
    }
    if (!run_open_) {
      open_run(nullptr);
    }
    if (!representation_) {
      decide(GroundRepresentation::plane);
    }
    Metadata end_metadata;
    if (*representation_ == GroundRepresentation::sphere) {
      end_metadata["misses_sphere"] = static_cast<std::int64_t>(misses_sphere_);
    }
    emit(RecordType::run_end, std::move(end_metadata));
    run_open_ = false;
    runh_.reset();
    telpos_.reset();
  }

  void check_telescope_indices(int array, int tel) const {
    if (array < 0 || tel < 0) {
      throw std::runtime_error("photon block with a negative array or telescope index");
    }
  }

  // Telescope record of the current event; returns the telescope id (1-based).
  std::int64_t begin_telescope(int array, int tel) {
    const auto telescope_id = static_cast<std::int64_t>(tel) + 1;
    Metadata tel_meta;
    tel_meta["telescope_id"] = telescope_id;
    tel_meta["frame"] = static_cast<std::int64_t>(Frame::ground);
    tel_meta["array_index"] = static_cast<std::int64_t>(array);
    if (telpos_ && tel >= 0 && static_cast<std::size_t>(tel) < telpos_->x.size()) {
      const auto i = static_cast<std::size_t>(tel);
      tel_meta["x_m"] = cm_to_m(telpos_->x[i]);
      tel_meta["y_m"] = cm_to_m(telpos_->y[i]);
      tel_meta["z_m"] = cm_to_m(telpos_->z[i]);
      tel_meta["r_m"] = cm_to_m(telpos_->r[i]);
    }
    emit(RecordType::telescope, std::move(tel_meta));
    return telescope_id;
  }

  void end_telescope(std::int64_t telescope_id, float header_photons) {
    Metadata end_meta;
    end_meta["telescope_id"] = telescope_id;
    // Photon sum from the block header (float, as CORSIKA wrote it).
    end_meta["photons"] = static_cast<double>(header_photons);
    emit(RecordType::telescope_end, std::move(end_meta));
  }

  // Photon chunks of one telescope: `total` rows in index order, row src of the input goes to
  // row dst of the chunk through fill(chunk, dst, src). Indices continue across array reuses.
  template <class Fill>
  void emit_rows(std::int64_t telescope_id, int tel, std::size_t total, std::uint32_t mask,
                 Fill fill) {
    std::uint64_t &first_index = next_photon_index_[tel];
    for (std::size_t start = 0; start < total; start += max_chunk_rows) {
      const std::size_t rows = std::min<std::size_t>(max_chunk_rows, total - start);
      Record record;
      record.type = RecordType::photons;
      record.chunk.telescope_id = static_cast<std::uint32_t>(telescope_id);
      record.chunk.frame = Frame::ground;
      record.chunk.first_index = first_index;
      record.chunk.resize(rows, mask);
      for (std::size_t i = 0; i < rows; ++i) {
        fill(record.chunk, i, start + i);
      }
      out(std::move(record));
      first_index += rows;
    }
  }

  // One telescope of the current event (1205): telescope record, photon chunks, end.
  void emit_photons(const TelPhotons &p) {
    check_telescope_indices(p.array, p.tel);
    claim(GroundRepresentation::plane);
    const std::int64_t telescope_id = begin_telescope(p.array, p.tel);
    constexpr std::uint32_t mask = column_bit(Column::x) | column_bit(Column::y) |
                                   column_bit(Column::cx) | column_bit(Column::cy) |
                                   column_bit(Column::t) | column_bit(Column::zem) |
                                   column_bit(Column::w) | column_bit(Column::wavelength);
    emit_rows(telescope_id, p.tel, p.bunches.size(), mask,
              [&p](PhotonChunk &chunk, std::size_t dst, std::size_t src) {
                const Bunch &b = p.bunches[src];
                chunk[Column::x][dst] = corsika_cm_to_m(b.x);
                chunk[Column::y][dst] = corsika_cm_to_m(b.y);
                chunk[Column::cx][dst] = b.cx;
                chunk[Column::cy][dst] = b.cy;
                chunk[Column::t][dst] = b.ctime;
                chunk[Column::zem][dst] = corsika_cm_to_m(b.zem);
                chunk[Column::w][dst] = b.photons;
                chunk[Column::wavelength][dst] = b.lambda;
              });
    end_telescope(telescope_id, p.photons);
  }

  // One telescope of the current event (1217) in a sphere run: each bunch becomes its G sphere row
  // (sphere_bunch) or counts as a miss. The telescope's radius and altitude come from TELPOS and
  // RUNH; the zem column is the flat-geometry emission altitude.
  void emit_photons3d(const TelPhotons3d &p) {
    check_telescope_indices(p.array, p.tel);
    const auto i = static_cast<std::size_t>(p.tel);
    if (!telpos_ || i >= telpos_->r.size()) {
      throw std::runtime_error("3D photon bunches of telescope " + std::to_string(p.tel + 1) +
                               " have no TELPOS entry (sphere radius)");
    }
    if (!runh_) {
      throw std::runtime_error("3D photon bunches need the RUNH block (observation level)");
    }
    claim(GroundRepresentation::sphere);
    const double radius_m = cm_to_m(telpos_->r[i]);
    const double altitude_m =
        cm_to_m(static_cast<double>(header_word(*runh_, corsika_word::runh_level_1, "RUNH"))) +
        cm_to_m(telpos_->z[i]);

    std::vector<detail::SphereRowValues> rows;
    rows.reserve(p.bunches.size());
    for (const Bunch3d &b : p.bunches) {
      const SphereBunch s = sphere_bunch(b, radius_m, altitude_m);
      if (!s.hits) {
        ++misses_sphere_;
        continue;
      }
      rows.push_back({static_cast<float>(s.x), static_cast<float>(s.y), static_cast<float>(s.z),
                      b.cx, b.cy, b.cz, static_cast<float>(s.time_ns), static_cast<float>(s.zem_m),
                      b.photons, static_cast<float>(s.path), b.lambda});
    }

    constexpr std::uint32_t mask =
        column_bit(Column::x) | column_bit(Column::y) | column_bit(Column::z) |
        column_bit(Column::cx) | column_bit(Column::cy) | column_bit(Column::cz) |
        column_bit(Column::t) | column_bit(Column::zem) | column_bit(Column::w) |
        column_bit(Column::path) | column_bit(Column::wavelength);
    const std::int64_t telescope_id = begin_telescope(p.array, p.tel);
    emit_rows(telescope_id, p.tel, rows.size(), mask,
              [&rows](PhotonChunk &chunk, std::size_t dst, std::size_t src) {
                const detail::SphereRowValues &r = rows[src];
                chunk[Column::x][dst] = r.x;
                chunk[Column::y][dst] = r.y;
                chunk[Column::z][dst] = r.z;
                chunk[Column::cx][dst] = r.cx;
                chunk[Column::cy][dst] = r.cy;
                chunk[Column::cz][dst] = r.cz;
                chunk[Column::t][dst] = r.t;
                chunk[Column::zem][dst] = r.zem;
                chunk[Column::w][dst] = r.w;
                chunk[Column::path][dst] = r.path;
                chunk[Column::wavelength][dst] = r.wavelength;
              });
    end_telescope(telescope_id, p.photons);
  }

  EventioReader reader_;
  const std::function<void(const Record &)> &sink_;
  std::deque<EventioBlock> queue_;
  std::vector<EventioBlock> pending_;
  std::optional<std::vector<float>> runh_;
  std::optional<TelPos> telpos_;
  bool run_open_{};
  bool event_open_{};
  std::map<int, std::uint64_t> next_photon_index_;     // per telescope, reset at each event
  bool rune_seen_{};                                   // a RUNE block was read since the last RUNH
  std::optional<GroundRepresentation> representation_; // of the open run, once known
  std::vector<Record> held_;                           // run records waiting for representation_
  std::uint64_t misses_sphere_{}; // 1217 bunches of the open run that miss their sphere
};

} // namespace detail

// Convert a CORSIKA IACT EventIO stream (long, compact or split photon blocks, or 3D photon
// bunches 1217) to photon-stream records. The stream and end records are not emitted; the
// caller writes them around this. Every non-photon top-level block becomes a blob record with its
// raw bytes. Throws std::runtime_error on malformed input.
//
// Representation: a run with 1205 blocks is a plane run (no ground_representation key, as before);
// a run with 1217 blocks is a sphere run ("sphere", G rows on or inside the telescope sphere, with
// columns x y z cx cy cz t zem w path wavelength; see SphereBunch). A run may not mix the two.
// Sphere runs count the bunches missing their sphere in the run end (misses_sphere). Heights: zem
// is observation level + TELPOS z + E_z, flat geometry, as for all CORSIKA heights here (no
// curvature conversion; the run key curvature_applied tells whether CORSIKA applied it).
inline void corsika_to_records(std::istream &in, const std::function<void(const Record &)> &sink) {
  detail::CorsikaToRecords converter{in, sink};
  converter.run();
}

// Writes records as an IACT EventIO stream for sim_telarray: blobs verbatim,
// photon chunks of each telescope as one long-format 1205 block, one 1204
// container per array reuse (normal mode). Photon chunks must be in the ground
// frame and in index order (first_index continues the previous chunk of the
// same telescope); the stream units (m) are converted back to cm in double.
//
// Split-mode markers (1213, 1214) are dropped on purpose: the photons are
// written in normal mode, which has no markers. They are blobs in the stream.
class SimtelWriter {
public:
  explicit SimtelWriter(std::ostream &out) : out_{out} {}

  void write(const Record &record) {
    switch (record.type) {
    case RecordType::telescope: {
      if (tel_open_) {
        throw std::runtime_error("telescope record before telescope_end");
      }
      const auto array_key = record.metadata.find("array_index");
      const std::int64_t array64 = array_key == record.metadata.end()
                                       ? 0
                                       : detail::metadata_int(record.metadata, "array_index");
      if (array64 < 0 || array64 > 32767) {
        throw std::out_of_range("array_index must be 0..32767");
      }
      const auto array = static_cast<int>(array64);
      if (array != array_) {
        flush_array();
      }
      array_ = array;
      const std::int64_t id64 = detail::metadata_int(record.metadata, "telescope_id");
      if (id64 < 1 || id64 > 32768) {
        throw std::out_of_range("telescope_id must be 1..32768 (sim_telarray numbering)");
      }
      telescope_id_ = static_cast<int>(id64);
      tel_open_ = true;
      // The 1205 header (12 bytes) is filled in when the telescope ends.
      bunch_bytes_.assign(photons_header_size, 0);
      photon_sum_ = 0.0;
      return;
    }
    case RecordType::event:
      begin_event();
      flush_array();
      break;
    case RecordType::photons:
      append(record.chunk);
      return;
    case RecordType::telescope_end:
      finish_telescope(record.metadata);
      return;
    case RecordType::efficiency_status:
      return; // efficiency status: not written to sim_telarray input
    case RecordType::blob:
      if (record.blob_eventio_type == iact::telarray_head ||
          record.blob_eventio_type == iact::telarray_end) {
        return; // split-mode marker, omitted in normal-mode output; keeps the container open
      }
      flush_array();
      break;
    default:
      flush_array();
      break;
    }
    if (record.type == RecordType::blob) {
      out_.write(reinterpret_cast<const char *>(record.blob.data()),
                 static_cast<std::streamsize>(record.blob.size()));
      if (!out_) {
        throw std::runtime_error("EventIO write failed");
      }
    }
  }

  // Close the open array container. Any non-photon record closes it too;
  // call this after the last record if the stream does not end with one.
  void finish() { flush_array(); }

  // A new event: photon indices of each telescope start again at 0.
  void begin_event() { next_index_.clear(); }

  // Photons of the open telescope record (ground frame), without a Record copy.
  void write_chunk(const PhotonChunk &chunk) { append(chunk); }

private:
  void append(const PhotonChunk &chunk) {
    if (!tel_open_ || chunk.telescope_id != static_cast<std::uint32_t>(telescope_id_)) {
      throw std::runtime_error("photon chunk outside its telescope record");
    }
    if (chunk.frame != Frame::ground) {
      throw std::runtime_error("sim_telarray input needs ground-frame photons");
    }
    chunk.validate();
    std::uint64_t &next_index = next_index_[telescope_id_];
    if (chunk.first_index != next_index) {
      throw std::runtime_error("photon chunk first_index " + std::to_string(chunk.first_index) +
                               " where " + std::to_string(next_index) +
                               " is expected (gap, overlap or reordering)");
    }
    const std::size_t n = chunk.rows();
    const auto &x = chunk[Column::x];
    const auto &y = chunk[Column::y];
    const auto &cx = chunk[Column::cx];
    const auto &cy = chunk[Column::cy];
    const auto &t = chunk[Column::t];
    const auto &zem = chunk[Column::zem];
    const auto &w = chunk[Column::w];
    const auto &lambda = chunk[Column::wavelength];
    if (n > 0 && zem.empty()) {
      throw std::runtime_error("sim_telarray input needs the zem column");
    }
    // Each bunch is encoded in place as its 1205 record (32 bytes):
    // x y cx cy t zem photons lambda, float32 little-endian, cm and ns as CORSIKA.
    const std::size_t at = bunch_bytes_.size();
    bunch_bytes_.resize(at + n * photon_bunch_size);
    std::uint8_t *dst = bunch_bytes_.data() + at;
    for (std::size_t i = 0; i < n; ++i, dst += photon_bunch_size) {
      detail::store_f32_le(dst, corsika_m_to_cm(x[i]));
      detail::store_f32_le(dst + 4, corsika_m_to_cm(y[i]));
      detail::store_f32_le(dst + 8, cx[i]);
      detail::store_f32_le(dst + 12, cy[i]);
      detail::store_f32_le(dst + 16, t[i]);
      detail::store_f32_le(dst + 20, corsika_m_to_cm(zem[i]));
      detail::store_f32_le(dst + 24, w[i]);
      detail::store_f32_le(dst + 28, lambda.empty() ? 0.0F : lambda[i]);
      photon_sum_ += w[i];
    }
    next_index += n;
  }

  void finish_telescope(const Metadata &metadata) {
    if (!tel_open_) {
      throw std::runtime_error("telescope_end without telescope record");
    }
    const std::int64_t id = detail::metadata_int(metadata, "telescope_id");
    if (id != telescope_id_) {
      throw std::runtime_error("telescope_end for telescope " + std::to_string(id) +
                               " while telescope " + std::to_string(telescope_id_) + " is open");
    }
    float photons = 0.0F;
    const auto it = metadata.find("photons");
    if (it != metadata.end()) {
      photons = static_cast<float>(std::get<double>(it->second));
    } else {
      photons = static_cast<float>(photon_sum_);
    }
    // Header of the 1205 block (make_photons_block layout): array and telescope as int16,
    // photon sum as float32, bunch count as int32; the bunches follow.
    const auto count =
        static_cast<std::int32_t>((bunch_bytes_.size() - photons_header_size) / photon_bunch_size);
    std::uint8_t *head = bunch_bytes_.data();
    const auto array16 = static_cast<std::uint16_t>(array_);
    const auto tel16 = static_cast<std::uint16_t>(telescope_id_ - 1);
    const auto count32 = static_cast<std::uint32_t>(count);
    for (std::size_t i = 0; i < 2; ++i) {
      head[i] = static_cast<std::uint8_t>(array16 >> (8 * i));
      head[2 + i] = static_cast<std::uint8_t>(tel16 >> (8 * i));
    }
    detail::store_f32_le(head + 4, photons);
    for (std::size_t i = 0; i < 4; ++i) {
      head[8 + i] = static_cast<std::uint8_t>(count32 >> (8 * i));
    }
    children_.push_back(make_block(iact::photons, 0,
                                   static_cast<std::uint32_t>(1000 * array_ + telescope_id_ - 1),
                                   std::move(bunch_bytes_)));
    bunch_bytes_.clear();
    tel_open_ = false;
  }

  void flush_array() {
    if (tel_open_) {
      throw std::runtime_error("array closed while a telescope is open (no telescope_end)");
    }
    if (children_.empty()) {
      return;
    }
    write_container(out_, iact::telarray, 0, static_cast<std::uint32_t>(array_), children_);
    children_.clear();
  }

  static constexpr std::size_t photons_header_size = 12; // int16, int16, float32, int32
  static constexpr std::size_t photon_bunch_size = 32;   // 8 float32 (long format)

  std::ostream &out_;
  int array_{};
  int telescope_id_{};
  bool tel_open_{};
  // first_index the next chunk of each telescope must have in the current event; it
  // continues across array reuses.
  std::map<int, std::uint64_t> next_index_;
  std::vector<std::uint8_t> bunch_bytes_; // 1205 block of the open telescope (header first)
  double photon_sum_{};                   // sum of the weights of the open telescope
  std::vector<EventioBlock> children_;
};

} // namespace obdeect::stream
