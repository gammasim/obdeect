#pragma once

// Atmosphere stage of a photon stream (docs/propagation/propagation.md). Records are read with
// obdeect::stream::Reader, passed to Propagator::process, and written with obdeect::stream::Writer.
//
// Per photon chunk, in this order (Propagator::photons):
//   1. sky-frame photons only: propagation to the telescope plane, bending, air delay, frame G
//   (sky.hpp)
//   2. wavelength assignment for rows without a wavelength (wavelength.hpp)
//   3. atmospheric extinction p_atm(lambda, zem, airmass) (extinction.hpp)
//   4. envelope E(lambda), the product of the envelope factors, then thinning or weighting
//      (thinning.hpp)
// Optional single scattering (scattering.hpp) runs on the primaries of an event before step 3.
// Its scattered photons go through step 1 (sky to G), then step 4 with the envelope and the
// in-telescope factor exp(-od_focus) only (Propagator::emit_scattered). An event is held until its
// end, see Propagator::flush_event. Each effect of the efficiency status is listed in
// make_efficiency_status(); a new effect adds a step in Propagator::photons and an entry in
// make_efficiency_status().

#include "obdeect/propagation/curved_refraction.hpp"
#include "obdeect/propagation/dispersion.hpp"
#include "obdeect/propagation/extinction.hpp"
#include "obdeect/propagation/pipeline.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/scattering.hpp"
#include "obdeect/propagation/sky.hpp"
#include "obdeect/propagation/sphere.hpp"
#include "obdeect/propagation/streams.hpp"
#include "obdeect/propagation/thinning.hpp"
#include "obdeect/propagation/timing.hpp"
#include "obdeect/propagation/transmission.hpp"
#include "obdeect/propagation/wavelength.hpp"

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <istream>
#include <iterator>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace obdeect::propagation {

enum class OutputMode { thinned, weighted };

// Earth curvature of the atmosphere (extinction only; see docs/propagation/extinction.md). The
// modes are explicit: the input run key curvature_applied says which geometry the photons have.
//   none:      plane-parallel atmosphere. Input with curvature_applied = 1 is rejected.
//   corsika:   photons from CORSIKA with the CURVED option (curvature_applied = 1, frame G, zem
//              above sea level). Only the airmass is replaced by the curved slant depth.
//   spherical: photons not from CORSIKA CURVED (curvature_applied = 0). The stage applies the
//              curvature to the sky-to-ground geometry (emission altitudes, local vertical).
enum class Curvature { none, corsika, spherical };

// "none", "corsika" or "spherical"; throws std::invalid_argument for any other text.
inline Curvature parse_curvature(std::string_view text) {
  if (text == "none") {
    return Curvature::none;
  }
  if (text == "corsika") {
    return Curvature::corsika;
  }
  if (text == "spherical") {
    return Curvature::spherical;
  }
  throw std::invalid_argument("curvature must be none, corsika or spherical: '" +
                              std::string(text) + "'");
}

// Ground representation of the output (docs/propagation/propagation.md, "Sphere representation"):
// plane (default, rows at z = 0 relative to the telescope) or sphere (rows on the telescope sphere,
// columns path and cz). "plane" or "sphere"; throws std::invalid_argument for other text.
inline obdeect::stream::GroundRepresentation parse_ground_representation(std::string_view text) {
  if (text == "plane") {
    return obdeect::stream::GroundRepresentation::plane;
  }
  if (text == "sphere") {
    return obdeect::stream::GroundRepresentation::sphere;
  }
  throw std::invalid_argument("ground representation must be plane or sphere: '" +
                              std::string(text) + "'");
}

// An envelope factor f^env(lambda) of a later efficiency (efficiency_status state `envelope`).
// The stage applies the product of all envelope factors. name: [A-Za-z0-9_-]+.
struct EnvelopeFactor {
  std::string name;
  std::vector<double> wavelength_nm;
  std::vector<double> value;
  std::string source; // provenance, e.g. the file name
};

struct PropagatorConfig {
  // Required unless geometry_only. Its H2 must match site_altitude_m (within 25 m).
  std::optional<TransmissionTable> transmission;
  // Geometry only: sky to G (refraction, air delay) and wavelength assignment, but no
  // extinction, envelope or thinning. All rows keep their weight; the efficiency status records
  // `atmosphere` as pending. Needs weighted mode, no transmission table, no envelope, no
  // scattering, no store_p and no in-telescope path.
  bool geometry_only{false};
  std::string transmission_source; // efficiency_status provenance; default "transmission table"
  double site_altitude_m{};        // observation level [m a.s.l.]
  // Optical path inside a telescope [m] (the od_focus part of the telescope depths).
  // default_pathlen_m applies to every telescope; telescope_pathlen_m (key: telescope_id,
  // 1-based) overrides it for one telescope. 0: no path inside the telescope.
  double default_pathlen_m{0.0};
  std::map<std::uint32_t, double> telescope_pathlen_m;
  // Optional: air delay for sky-frame photons, and refraction (bending) unless the run says that
  // refraction is already applied.
  std::optional<AtmosphereProfile> profile;
  std::string profile_source;
  // Dispersion (docs/propagation/refraction.md): reference wavelength [nm] of the profile's n-1.
  // Needs the profile; the bent sky-frame photons of a run with refraction_applied = 0 are scaled
  // to their wavelength (200 to 2000 nm). Rejected for runs that are already refracted.
  std::optional<double> dispersion_reference_nm;
  std::vector<EnvelopeFactor> envelope;
  // Wavelengths drawn from 1/lambda^2 for rows without one. Default: the run keys
  // wavelength_min_nm and wavelength_max_nm.
  std::optional<WavelengthRange> wavelength_range;
  std::uint64_t seed{};
  OutputMode mode{OutputMode::thinned};
  bool store_p{false}; // keep the acceptance column p in the output
  // Optional single scattering of the sky-frame photons (docs/propagation/scattering.md), off by
  // default. The table switches it on: scattering is on if scattering_table is set
  // (scattering.enabled is set from it; a true flag without a table is rejected). The run needs
  // telescope_r_m (sphere radii).
  ScatteringConfig scattering;
  std::optional<ScatteringTable> scattering_table;
  std::string scattering_source; // efficiency_status provenance, e.g. the table file name
  // Earth curvature (see Curvature). earth_radius_m: Earth radius R [m], used by corsika (the
  // efficiency status) and spherical. Must be finite and above 1e6 m.
  Curvature curvature{Curvature::none};
  double earth_radius_m{obdeect::stream::earth_radius_m};
  // Ground representation of the output (see parse_ground_representation). Sphere needs the run key
  // telescope_r_m (one radius per telescope). Input sphere runs are kept as sphere rows and need
  // this setting too. Sphere rows take extinction and scattering (docs/propagation/extinction.md,
  // "Sphere rows"; docs/propagation/scattering.md, "Sphere runs").
  obdeect::stream::GroundRepresentation ground_representation{
      obdeect::stream::GroundRepresentation::plane};
  // Threads for the extinction of each chunk and for the curved tables (1: in place). The output
  // does not depend on it.
  unsigned threads{1};
  // Telescopes to keep (1-based ids). Blocks of other telescopes (telescope, photons,
  // telescope_end and their efficiency status) are dropped; their photons are not decoded when the
  // stream is read with the same selection (propagate_stream). With scattering, the kept
  // telescopes are the scattering targets; the blocks of other telescopes are still read for
  // primaries (scattering_primaries) but not written. Default: all telescopes.
  std::optional<std::vector<std::uint32_t>> telescopes;
  // Telescopes whose photons are scattered (1-based ids; needs scattering). Default: the blocks
  // of all telescopes.
  std::optional<std::vector<std::uint32_t>> scattering_primaries;
};

// Counts of a run. Every photon of the input has one fate (docs/propagation/propagation.md,
// "Accounting"): the primaries are written, dropped by reason, or removed by the acceptance
// (thinned_out); the scattered photons are written, dropped by reason or removed by the acceptance
// (scattered_thinned_out).
struct PropagationSummary {
  std::uint64_t photons_in{};  // rows of the input photon records (selected telescopes only)
  std::uint64_t photons_out{}; // rows written (thinned: surviving photons, weighted: w > 0)
  std::uint64_t dropped{};     // primaries that do not reach the telescope plane (sphere)
  // The dropped primaries by reason (DropCounts): dropped == dropped_by_reason.total().
  DropCounts dropped_by_reason{};
  double weight_out{}; // sum of w of the output
  // Weight removed by the acceptance draw of the primaries: thinned, the photons that do not
  // survive (sum of w - survivors over the reaching rows; equal to a count for rows of w = 1, where
  // photons_in = photons_out + dropped + thinned_out); weighted, the weight w (1 - p) removed.
  double thinned_out{};
  // Scattered photons (single scattering) that do not reach the telescope, by reason, and the
  // weight removed by their acceptance (as thinned_out).
  DropCounts scattered_dropped_by_reason{};
  std::uint64_t scattered_dropped{}; // == scattered_dropped_by_reason.total()
  double scattered_thinned_out{};
};

// Reads a two-column text table (wavelength nm, value; '#' starts a comment).
inline EnvelopeFactor parse_envelope_text(std::string_view text, std::string name,
                                          std::string source) {
  EnvelopeFactor factor{std::move(name), {}, {}, std::move(source)};
  std::size_t pos = 0;
  while (pos < text.size()) {
    std::size_t end = text.find('\n', pos);
    if (end == std::string_view::npos) {
      end = text.size();
    }
    const std::string_view line = text.substr(pos, end - pos);
    pos = std::min(end + 1, text.size());
    const auto tokens = detail::split_whitespace(detail::trim_comment(line));
    if (tokens.empty()) {
      continue;
    }
    if (tokens.size() != 2) {
      throw std::invalid_argument("efficiency bound table: each line needs wavelength and value");
    }
    factor.wavelength_nm.push_back(detail::parse_number(tokens[0], "wavelength"));
    factor.value.push_back(detail::parse_number(tokens[1], "value"));
  }
  // Envelope's constructor validates the table (order, range of values).
  (void)Envelope(factor.wavelength_nm, factor.value);
  return factor;
}

// Reads an envelope file (see parse_envelope_text).
inline EnvelopeFactor read_envelope_file(const std::string &name, const std::string &path) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    throw std::runtime_error("cannot read envelope file: " + path);
  }
  std::ostringstream buffer;
  buffer << file.rdbuf();
  return parse_envelope_text(buffer.str(), name, path);
}

// True if both names are the same file (same spelling, symlink, hard link, or the same
// canonical path when the output does not exist yet). Used to refuse an output that would
// truncate the input.
inline bool same_file(const std::string &a, const std::string &b) {
  namespace fs = std::filesystem;
  if (a == b) {
    return true;
  }
  std::error_code ec_a;
  std::error_code ec_b;
  const fs::path canonical_a = fs::weakly_canonical(a, ec_a);
  const fs::path canonical_b = fs::weakly_canonical(b, ec_b);
  if (!ec_a && !ec_b && canonical_a == canonical_b) {
    return true;
  }
  std::error_code ec;
  return fs::equivalent(a, b, ec); // true only if both exist and are the same file
}

// Telescope ids of a selection (1-based, at most max_telescope); throws std::invalid_argument.
// what names the selection in the message.
inline void check_telescope_ids(const std::vector<std::uint32_t> &ids, const char *what) {
  for (const auto id : ids) {
    if (id < 1 || id > obdeect::stream::max_telescope) {
      throw std::invalid_argument(std::string(what) + ": telescope ids are 1-based");
    }
  }
}

// Sorts a list of telescope ids and removes the repeats.
inline void sort_unique_telescope_ids(std::vector<std::uint32_t> &ids) {
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
}

class Propagator {
public:
  explicit Propagator(PropagatorConfig config)
      : config_(std::move(config)), pool_(config_.threads) {
    if (config_.geometry_only) {
      check_geometry_only();
    } else if (!config_.transmission) {
      throw std::invalid_argument("propagator needs a transmission table");
    }
    if (!std::isfinite(config_.site_altitude_m)) {
      throw std::invalid_argument("site altitude must be finite");
    }
    if (config_.transmission && std::abs(config_.transmission->observation_level_m() -
                                         config_.site_altitude_m) > altitude_tolerance_m) {
      throw std::invalid_argument("transmission table H2 does not match the site altitude");
    }
    if (config_.dispersion_reference_nm) {
      if (!config_.profile) {
        throw std::invalid_argument("dispersion needs an atmospheric profile");
      }
      (void)refractivity_e8(*config_.dispersion_reference_nm); // validates 200 to 2000 nm
    }
    if (config_.wavelength_range) {
      config_.wavelength_range->validate();
    }
    check_pathlen(config_.default_pathlen_m);
    for (const auto &[id, pathlen] : config_.telescope_pathlen_m) {
      if (id < 1 || id > obdeect::stream::max_telescope) {
        throw std::invalid_argument("telescope_pathlen_m: telescope ids are 1-based");
      }
      check_pathlen(pathlen);
    }
    if (config_.scattering.enabled && !config_.scattering_table) {
      throw std::invalid_argument("scattering needs a scattering table");
    }
    config_.scattering.enabled = scattering_on(); // the scatter kernel reads the flag
    if (scattering_on() && (config_.scattering.points_per_photon < 1 ||
                            config_.scattering.points_per_photon > max_points_per_photon)) {
      throw std::invalid_argument("scattering needs 1 to 1024 points per photon");
    }
    check_curvature();
    if (config_.telescopes) {
      check_telescopes();
    }
    if (config_.scattering_primaries) {
      check_scattering_primaries();
    }
    // Air delay of the scattered photons: the run's profile (points at config_.profile, which
    // lives as long as this Propagator).
    if (config_.profile) {
      config_.scattering.air_profile = &*config_.profile;
    }
    for (std::size_t i = 0; i < config_.envelope.size(); ++i) {
      const auto &name = config_.envelope[i].name;
      if (!valid_factor_name(name)) {
        throw std::invalid_argument("invalid envelope factor name: '" + name + "'");
      }
      for (std::size_t j = 0; j < i; ++j) {
        if (config_.envelope[j].name == name) {
          throw std::invalid_argument("duplicate envelope factor: " + name);
        }
      }
      (void)Envelope(config_.envelope[i].wavelength_nm, config_.envelope[i].value);
    }
  }

  // The propagator is not copied or moved: config_.scattering.air_profile points into config_, and
  // pool_ owns its threads.
  Propagator(const Propagator &) = delete;
  Propagator &operator=(const Propagator &) = delete;
  Propagator(Propagator &&) = delete;
  Propagator &operator=(Propagator &&) = delete;
  ~Propagator() = default;

  // Processes one input record. Records of the output go to sink(const obdeect::stream::Record &).
  // The efficiency status of a run is written before the first record that is not an efficiency
  // status after it. Throws std::invalid_argument if an incoming efficiency_status shows a factor
  // applied already.
  // Processes a record that the caller keeps (its photons are copied).
  template <class Sink> void process(const obdeect::stream::Record &record, Sink &&sink) {
    obdeect::stream::Record copy = record;
    process(copy, sink);
  }

  // Processes a record that may be changed: its photons are moved into the propagator.
  template <class Sink> void process(obdeect::stream::Record &record, Sink &&sink) {
    using obdeect::stream::RecordType;
    if (!scattering_on()) {
      process_record(record, sink);
      return;
    }
    // With scattering, an event is held until its end: the scattered photons of a telescope come
    // from the primaries of all telescopes of the event.
    if (in_event_) {
      if (record.type == RecordType::event_end) {
        event_records_.push_back(std::move(record));
        flush_event(sink);
        return;
      }
      const bool ends_event = record.type == RecordType::event || record.type == RecordType::run ||
                              record.type == RecordType::run_end ||
                              record.type == RecordType::stream ||
                              record.type == RecordType::stream_end;
      if (!ends_event) {
        event_records_.push_back(std::move(record));
        return;
      }
      flush_event(sink);
    }
    if (record.type == RecordType::event) {
      in_event_ = true;
      event_records_.assign(1, record);
      return;
    }
    process_record(record, sink);
  }

  // Writes the held event (if any) and the efficiency status still pending.
  template <class Sink> void finish(Sink &&sink) {
    if (in_event_) {
      flush_event(sink);
    }
    flush_efficiency_status(sink);
  }

  [[nodiscard]] const PropagationSummary &summary() const noexcept { return summary_; }

private:
  // Sum of the bunch weights of a chunk (0 without a w column).
  static double weight_sum(const obdeect::stream::PhotonChunk &chunk) {
    double sum = 0.0;
    for (const float weight : chunk[obdeect::stream::Column::w]) {
      sum += static_cast<double>(weight);
    }
    return sum;
  }

  static double weight_sum(const std::vector<obdeect::stream::PhotonChunk> &chunks) {
    double sum = 0.0;
    for (const auto &chunk : chunks) {
      sum += weight_sum(chunk);
    }
    return sum;
  }

  // With scattering: one telescope block per telescope and event. The scattered photons of a
  // telescope are written in its block; a second block would need its own scattered photons.
  static void check_one_block_per_telescope(const std::vector<obdeect::stream::Record> &records) {
    std::vector<std::int64_t> ids;
    for (const auto &record : records) {
      if (record.type == obdeect::stream::RecordType::telescope) {
        ids.push_back(
            checked_integer(record.metadata, "telescope_id", obdeect::stream::max_telescope));
      }
    }
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) {
      throw std::invalid_argument("single scattering needs one telescope block per telescope and "
                                  "event (CORSIKA array reuse is not supported with scattering)");
    }
  }

  // The photon blocks of a telescope in an event continue the index: first_index counts the
  // photons of the telescope in the event, and the random keys of the stages use it. A block must
  // start where the previous one of the telescope ended (in_index_ is reset at each event).
  void check_first_index(const obdeect::stream::PhotonChunk &input) {
    std::uint64_t &next = in_index_[telescope_.id];
    if (input.first_index != next) {
      throw std::invalid_argument("photon block of telescope " + std::to_string(telescope_.id) +
                                  " in event " + std::to_string(event_number_) +
                                  " starts at first_index " + std::to_string(input.first_index) +
                                  ", expected " + std::to_string(next) +
                                  " (the blocks of a telescope must continue the index)");
    }
    next += input.rows();
  }

  static constexpr std::int64_t max_u32 = 0xFFFFFFFFLL;
  static constexpr double min_earth_radius_m = 1e6;

  // One record of the stream (process without the event hold).
  template <class Sink> void process_record(obdeect::stream::Record &record, Sink &sink) {
    using obdeect::stream::RecordType;
    if (record.type == RecordType::efficiency_status) {
      check_incoming_efficiency_status(record.metadata);
      if (!telescope_.open && efficiency_status_pending_) {
        // Run status of the input: merged into the one run status of this stage (flush).
        incoming_run_status_ = merge_status(incoming_run_status_, record.metadata);
        return;
      }
      // A telescope status (after the telescope record) belongs to the dropped block if its
      // telescope is not selected.
      if (!telescope_.open || selected(telescope_.id)) {
        sink(record);
      }
      return;
    }
    flush_efficiency_status(sink);
    switch (record.type) {
    case RecordType::stream:
    case RecordType::event_end:
    case RecordType::run_end:
    case RecordType::blob:
      sink(record);
      break;
    case RecordType::run: {
      obdeect::stream::Record out = record;
      begin_run(record.metadata, out.metadata);
      sink(out);
      efficiency_status_pending_ = true;
      break;
    }
    case RecordType::event:
      event_number_ =
          static_cast<std::uint32_t>(checked_integer(record.metadata, "event_number", max_u32));
      have_event_ = true;
      out_index_.clear();
      in_index_.clear();
      sink(record);
      break;
    case RecordType::telescope: {
      obdeect::stream::Record out = record;
      open_telescope(record.metadata, out.metadata);
      if (selected(telescope_.id)) {
        sink(out);
      }
      break;
    }
    case RecordType::photons:
      photons(record.chunk, sink); // moves the chunk
      break;
    case RecordType::telescope_end: {
      const auto id =
          checked_integer(record.metadata, "telescope_id", obdeect::stream::max_telescope);
      if (!telescope_.open || static_cast<std::int64_t>(telescope_.id) != id) {
        throw std::invalid_argument("telescope_end of a telescope that is not open");
      }
      obdeect::stream::Record out = record;
      out.metadata["photons"] = telescope_.weight_out;
      out_index_[telescope_.id] = telescope_.out_index;
      telescope_.open = false;
      if (selected(telescope_.id)) {
        sink(out);
      }
      break;
    }
    case RecordType::efficiency_status:
    case RecordType::stream_end:
      break;
    }
  }

  struct Run {
    double observation_m{};
    bool refraction_applied{};
    bool sphere{};                // output in the sphere representation (config)
    bool sphere_input{};          // the input run is in the sphere representation
    std::vector<double> radius_m; // telescope sphere radii [m] (sphere mode)
    std::vector<double> x_m, y_m, z_m;
    std::vector<double> pathlen_m; // path inside each telescope [m], index = telescope_id - 1
    WavelengthRange wavelengths{};
    Envelope envelope;
    obdeect::stream::Metadata efficiency_status;
  };

  struct Telescope {
    bool open{};
    std::uint32_t id{};
    double altitude_m{}; // a.s.l. of the telescope plane
    double pathlen_m{};  // path inside the telescope [m]
    double x_m{}, y_m{}, z_m{};
    double radius_m{};                // sphere radius [m] (sphere mode)
    double nx{0.0}, ny{0.0}, nz{1.0}; // local vertical of the telescope (frame G)
    const TelescopeDepths *depths{};  // cached per (altitude, pathlen) (depths_cache_)
    std::uint64_t out_index{};        // output rows so far (first_index of the next chunk)
    double weight_out{};
  };

  // Geometry-only runs apply no extinction: reject every setting that belongs to it.
  void check_geometry_only() const {
    if (config_.transmission) {
      throw std::invalid_argument("geometry-only run takes no transmission table");
    }
    if (!config_.envelope.empty()) {
      throw std::invalid_argument("geometry-only run takes no envelope (no thinning)");
    }
    if (config_.mode != OutputMode::weighted) {
      throw std::invalid_argument("geometry-only run keeps all photons: mode must be weighted");
    }
    if (scattering_on()) {
      throw std::invalid_argument("geometry-only run cannot be combined with scattering");
    }
    if (config_.store_p) {
      throw std::invalid_argument("geometry-only run has no acceptance column p to store");
    }
    if (config_.default_pathlen_m != 0.0 || !config_.telescope_pathlen_m.empty()) {
      throw std::invalid_argument("geometry-only run takes no in-telescope path");
    }
  }

  static void check_pathlen(double pathlen_m) {
    if (!std::isfinite(pathlen_m) || pathlen_m < 0.0) {
      throw std::invalid_argument("pathlen must be finite and not negative");
    }
  }

  // Telescope selection: at least one id, 1-based; sorted and without repeats.
  void check_telescopes() {
    auto &ids = *config_.telescopes;
    if (ids.empty()) {
      throw std::invalid_argument("telescope selection is empty");
    }
    check_telescope_ids(ids, "telescope selection");
    sort_unique_telescope_ids(ids);
  }

  // Scattering primaries: needs scattering; at least one id, 1-based; sorted and without repeats.
  // Ids of the run are checked in begin_run.
  void check_scattering_primaries() {
    auto &ids = *config_.scattering_primaries;
    if (!scattering_on()) {
      throw std::invalid_argument("scattering_primaries needs scattering");
    }
    if (ids.empty()) {
      throw std::invalid_argument("scattering_primaries is empty");
    }
    check_telescope_ids(ids, "scattering_primaries");
    sort_unique_telescope_ids(ids);
  }

  // Scattering is on if the config has a scattering table (the table is the switch).
  [[nodiscard]] bool scattering_on() const noexcept { return config_.scattering_table.has_value(); }

  // Primary photons of telescope id are scattered (all blocks by default).
  [[nodiscard]] bool is_scattering_primary(std::uint32_t telescope_id) const {
    return !config_.scattering_primaries ||
           std::binary_search(config_.scattering_primaries->begin(),
                              config_.scattering_primaries->end(), telescope_id);
  }

  // Array reuse: the scattering needs one array_index per event (the primaries of an event are
  // scattered on the same core positions).
  static void check_single_array_index(const std::vector<obdeect::stream::Record> &records) {
    std::optional<std::int64_t> first;
    for (const auto &record : records) {
      if (record.type != obdeect::stream::RecordType::telescope) {
        continue;
      }
      const auto it = record.metadata.find("array_index");
      if (it == record.metadata.end()) {
        continue;
      }
      const auto *value = std::get_if<std::int64_t>(&it->second);
      if (value == nullptr) {
        throw std::invalid_argument("telescope key array_index is not an integer");
      }
      if (!first) {
        first = *value;
      } else if (*value != *first) {
        throw std::invalid_argument("single scattering needs one array_index per event; run "
                                    "CORSIKA without array reuse (CSCAT 1 ...)");
      }
    }
  }

  [[nodiscard]] bool selected(std::uint32_t telescope_id) const {
    return !config_.telescopes ||
           std::binary_search(config_.telescopes->begin(), config_.telescopes->end(), telescope_id);
  }

  // Config-only checks of the curvature modes (the input is checked in begin_run).
  void check_curvature() const {
    if (!std::isfinite(config_.earth_radius_m) || !(config_.earth_radius_m > min_earth_radius_m)) {
      throw std::invalid_argument("earth radius must be finite and above 1e6 m");
    }
    if (config_.curvature == Curvature::none) {
      return;
    }
    if (scattering_on()) {
      throw std::invalid_argument("curved atmosphere cannot be combined with scattering");
    }
  }

  // The run key curvature_applied (optional: nullopt if the run lacks it) must match the mode.
  void check_run_curvature(std::optional<std::int64_t> applied) const {
    switch (config_.curvature) {
    case Curvature::none:
      if (applied == 1) {
        throw std::invalid_argument("input photons were produced with CORSIKA CURVED "
                                    "(curvature_applied = 1): select curvature mode 'corsika'");
      }
      break;
    case Curvature::corsika:
      if (applied != 1) {
        throw std::invalid_argument("curvature mode 'corsika' needs a run with curvature_applied = "
                                    "1 (CORSIKA CURVED option)");
      }
      break;
    case Curvature::spherical:
      if (!applied) {
        throw std::invalid_argument("curvature mode 'spherical' needs the run key "
                                    "curvature_applied (0: no curvature in the input)");
      }
      if (*applied == 1) {
        throw std::invalid_argument(
            "curvature is already applied in the input (curvature_applied = 1)");
      }
      break;
    }
  }

  static bool valid_factor_name(const std::string &name) {
    if (name.empty() || name == "atmosphere" || name == "wavelength" || name == "scattering") {
      return false;
    }
    return std::all_of(name.begin(), name.end(), [](char c) {
      return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
             c == '_' || c == '-';
    });
  }

  static double number_field(const obdeect::stream::Metadata &metadata, const char *key) {
    const auto it = metadata.find(key);
    if (it == metadata.end()) {
      throw std::invalid_argument(std::string("stream record lacks key '") + key + "'");
    }
    if (const auto *value = std::get_if<double>(&it->second)) {
      if (!std::isfinite(*value)) {
        throw std::invalid_argument(std::string("key '") + key + "' is not finite");
      }
      return *value;
    }
    if (const auto *value = std::get_if<std::int64_t>(&it->second)) {
      return static_cast<double>(*value);
    }
    throw std::invalid_argument(std::string("key '") + key + "' is not a number");
  }

  static std::int64_t checked_integer(const obdeect::stream::Metadata &metadata, const char *key,
                                      std::int64_t limit) {
    const auto it = metadata.find(key);
    if (it == metadata.end() || !std::holds_alternative<std::int64_t>(it->second)) {
      throw std::invalid_argument(std::string("stream record lacks integer key '") + key + "'");
    }
    const auto value = std::get<std::int64_t>(it->second);
    if (value < 0 || value > limit) {
      throw std::invalid_argument(std::string("key '") + key + "' out of range");
    }
    return value;
  }

  static const std::vector<double> &array_field(const obdeect::stream::Metadata &metadata,
                                                const char *key) {
    const auto it = metadata.find(key);
    if (it == metadata.end() || !std::holds_alternative<std::vector<double>>(it->second)) {
      throw std::invalid_argument(std::string("run record lacks array key '") + key + "'");
    }
    return std::get<std::vector<double>>(it->second);
  }

  // Efficiency status values that other stages wrote for factors this stage applies (never apply
  // twice).
  void check_incoming_efficiency_status(const obdeect::stream::Metadata &metadata) const {
    const auto state_of = [&metadata](const std::string &name) -> std::string {
      const auto it = metadata.find(name + ".state");
      if (it == metadata.end()) {
        return "pending";
      }
      if (const auto *text = std::get_if<std::string>(&it->second)) {
        return *text;
      }
      throw std::invalid_argument("efficiency_status state of " + name + " is not a string");
    };
    if (state_of("atmosphere") != "pending") {
      throw std::invalid_argument("atmosphere is already applied in the input efficiency status");
    }
    if (state_of("scattering") != "pending") {
      throw std::invalid_argument("scattering is already applied in the input efficiency status");
    }
    for (const auto &factor : config_.envelope) {
      if (state_of(factor.name) != "pending") {
        throw std::invalid_argument("envelope factor " + factor.name +
                                    " is already applied in the input efficiency status");
      }
    }
  }

  static int state_rank(const std::string &state) {
    if (state == "pending") {
      return 0;
    }
    return state == "envelope" ? 1 : 2;
  }

  static std::vector<std::string> factor_names(const obdeect::stream::Metadata &metadata) {
    std::vector<std::string> names;
    const auto it = metadata.find("factors");
    if (it == metadata.end()) {
      return names;
    }
    const auto *list = std::get_if<std::string>(&it->second);
    if (list == nullptr) {
      throw std::invalid_argument("efficiency status 'factors' is not a string");
    }
    std::size_t begin = 0;
    while (begin <= list->size()) {
      std::size_t end = list->find(',', begin);
      if (end == std::string::npos) {
        end = list->size();
      }
      const std::string token = list->substr(begin, end - begin);
      begin = end + 1;
      const auto first = token.find_first_not_of(" \t");
      if (first != std::string::npos) {
        names.push_back(token.substr(first, token.find_last_not_of(" \t") - first + 1));
      }
    }
    return names;
  }

  // Merges an efficiency status `add` into `base` (one run status record,
  // docs/propagation/propagation.md, "Efficiency status"): a factor that `base` has in the same
  // state is not written again; a factor that advances (pending -> envelope -> applied) takes the
  // keys of `add` (its own keys in `base` are replaced); a factor in a backward state is an error.
  // Factors of `base` that `add` does not name are kept.
  static obdeect::stream::Metadata merge_status(const obdeect::stream::Metadata &base,
                                                const obdeect::stream::Metadata &add) {
    obdeect::stream::Metadata merged = base;
    std::vector<std::string> names = factor_names(base);
    for (const auto &name : factor_names(add)) {
      const std::string state_key = name + ".state";
      const auto add_state_it = add.find(state_key);
      if (add_state_it == add.end() || !std::holds_alternative<std::string>(add_state_it->second)) {
        throw std::invalid_argument("efficiency status of " + name + " lacks a string state");
      }
      const std::string add_state = std::get<std::string>(add_state_it->second);
      const auto base_it = merged.find(state_key);
      if (base_it != merged.end()) {
        const auto *base_state = std::get_if<std::string>(&base_it->second);
        if (base_state == nullptr) {
          throw std::invalid_argument("efficiency status state of " + name + " is not a string");
        }
        if (*base_state == add_state) {
          continue; // already in the input in this state: not written again
        }
        if (state_rank(add_state) < state_rank(*base_state)) {
          throw std::invalid_argument("efficiency status of " + name + " goes back from " +
                                      *base_state + " to " + add_state);
        }
      } else {
        names.push_back(name);
      }
      const std::string prefix = name + ".";
      for (auto it = merged.begin(); it != merged.end();) {
        it = it->first.compare(0, prefix.size(), prefix) == 0 ? merged.erase(it) : std::next(it);
      }
      for (const auto &[key, value] : add) {
        if (key.compare(0, prefix.size(), prefix) == 0) {
          merged[key] = value;
        }
      }
    }
    std::string list;
    for (const auto &name : names) {
      list += (list.empty() ? "" : ",") + name;
    }
    merged["factors"] = list;
    return merged;
  }

  // Run keys and the efficiency status of this stage (docs/stream/efficiency-status.md).
  void begin_run(const obdeect::stream::Metadata &metadata, obdeect::stream::Metadata &out) {
    const auto &x_m = array_field(metadata, "telescope_x_m");
    const auto &y_m = array_field(metadata, "telescope_y_m");
    const auto &z_m = array_field(metadata, "telescope_z_m");
    if (y_m.size() != x_m.size() || z_m.size() != x_m.size()) {
      throw std::invalid_argument("telescope position arrays differ in length");
    }
    const double observation_m = number_field(metadata, "observation_level_m");
    if (std::abs(observation_m - config_.site_altitude_m) > altitude_tolerance_m) {
      throw std::invalid_argument("run observation level does not match the site altitude");
    }
    std::optional<std::int64_t> curvature_applied;
    if (metadata.count("curvature_applied") != 0) {
      curvature_applied = checked_integer(metadata, "curvature_applied", 1);
    }
    check_run_curvature(curvature_applied);
    if (config_.curvature == Curvature::spherical) {
      out["curvature_applied"] = std::int64_t{1}; // curvature is applied from here on
    }
    bool refraction_applied = false;
    if (metadata.count("refraction_applied") != 0) {
      refraction_applied = checked_integer(metadata, "refraction_applied", 1) == 1;
    } else if (config_.profile) {
      throw std::invalid_argument("run lacks refraction_applied (needed with a profile)");
    }
    const bool sphere_input = obdeect::stream::ground_representation(metadata) ==
                              obdeect::stream::GroundRepresentation::sphere;
    const bool sphere =
        config_.ground_representation == obdeect::stream::GroundRepresentation::sphere;
    config_.scattering.sphere_geometry = sphere;
    if (sphere_input && !sphere) {
      throw std::invalid_argument(
          "input run is in the sphere representation (ground_representation "
          "= sphere): select ground representation sphere");
    }
    std::vector<double> radius_m;
    if (sphere) {
      radius_m = array_field(metadata, "telescope_r_m");
      if (radius_m.size() != x_m.size()) {
        throw std::invalid_argument("sphere representation needs telescope_r_m, one per telescope");
      }
      for (const double radius : radius_m) {
        if (!(std::isfinite(radius) && radius > 0.0)) {
          throw std::invalid_argument("telescope_r_m must be finite and positive");
        }
      }
    }
    WavelengthRange wavelengths{};
    if (config_.wavelength_range) {
      wavelengths = *config_.wavelength_range;
    } else {
      wavelengths = {number_field(metadata, "wavelength_min_nm"),
                     number_field(metadata, "wavelength_max_nm")};
    }
    wavelengths.validate();
    if (config_.transmission) {
      config_.transmission->require_wavelength_range(wavelengths.lo_nm, wavelengths.hi_nm);
    }
    if (config_.dispersion_reference_nm) {
      if (refraction_applied) {
        throw std::invalid_argument(
            "dispersion needs photons that are not refracted yet (run key refraction_applied = 1: "
            "CORSIKA already refracted them); run without dispersion, then to-sky, or set "
            "refraction_applied = 0");
      }
      if (wavelengths.lo_nm < refractivity_min_nm || wavelengths.hi_nm > refractivity_max_nm) {
        throw std::invalid_argument("dispersion needs the wavelength range within 200 to 2000 nm");
      }
    }
    if (scattering_on()) {
      const auto &r_m = array_field(metadata, "telescope_r_m");
      if (r_m.size() != x_m.size()) {
        throw std::invalid_argument("single scattering needs telescope_r_m, one per telescope");
      }
      if (config_.scattering_primaries && config_.scattering_primaries->back() > x_m.size()) {
        throw std::invalid_argument("scattering_primaries names telescope " +
                                    std::to_string(config_.scattering_primaries->back()) +
                                    ", the run has " + std::to_string(x_m.size()));
      }
      sphere_.clear();
      for (std::size_t i = 0; i < x_m.size(); ++i) {
        sphere_.push_back({static_cast<std::uint32_t>(i + 1), x_m[i], y_m[i], z_m[i], r_m[i]});
      }
      scattered_.assign(x_m.size(), {});
    }
    if (config_.telescopes && config_.telescopes->back() > x_m.size()) {
      throw std::invalid_argument("telescope selection names telescope " +
                                  std::to_string(config_.telescopes->back()) + ", the run has " +
                                  std::to_string(x_m.size()));
    }
    std::vector<double> pathlen_m(x_m.size(), config_.default_pathlen_m);
    for (const auto &[id, pathlen] : config_.telescope_pathlen_m) {
      if (id > x_m.size()) {
        throw std::invalid_argument("telescope_pathlen_m names telescope " + std::to_string(id) +
                                    ", the run has " + std::to_string(x_m.size()));
      }
      pathlen_m[id - 1] = pathlen;
    }
    Run run{observation_m,
            refraction_applied,
            sphere,
            sphere_input,
            radius_m,
            x_m,
            y_m,
            z_m,
            pathlen_m,
            wavelengths,
            product_envelope(wavelengths),
            make_efficiency_status(wavelengths, pathlen_m)};
    if (config_.profile) {
      out["refraction_applied"] = std::int64_t{1}; // sky-frame photons are refracted from here on
    }
    if (sphere) {
      out["ground_representation"] = std::string("sphere");
    }
    run_ = std::move(run);
    telescope_.open = false;
    have_event_ = false;
  }

  // Product of the envelope factors on the integer wavelengths of the range.
  [[nodiscard]] Envelope product_envelope(const WavelengthRange &range) const {
    const double first = std::floor(range.lo_nm);
    const double last = std::ceil(range.hi_nm);
    std::vector<double> grid;
    for (double lambda = first; lambda <= last; lambda += 1.0) {
      grid.push_back(lambda);
    }
    std::vector<double> product(grid.size(), 1.0);
    std::vector<Envelope> factors;
    for (const auto &factor : config_.envelope) {
      if (factor.wavelength_nm.front() > first || factor.wavelength_nm.back() < last) {
        throw std::invalid_argument("envelope factor " + factor.name +
                                    " does not cover the wavelength range");
      }
      factors.emplace_back(factor.wavelength_nm, factor.value);
      for (std::size_t i = 0; i < grid.size(); ++i) {
        product[i] *= factors.back()(grid[i]);
      }
    }
    return Envelope(std::move(grid), std::move(product), std::move(factors));
  }

  [[nodiscard]] obdeect::stream::Metadata
  make_efficiency_status(const WavelengthRange &range, const std::vector<double> &pathlen_m) const {
    obdeect::stream::Metadata efficiency_status;
    std::string factors = "atmosphere,wavelength";
    if (config_.geometry_only) {
      efficiency_status["atmosphere.state"] = std::string("pending"); // applied by a later stage
    } else {
      efficiency_status["atmosphere.state"] = std::string("applied");
      efficiency_status["atmosphere.source"] = config_.transmission_source.empty()
                                                   ? std::string("transmission table")
                                                   : config_.transmission_source;
      efficiency_status["atmosphere.pathlen_m"] = pathlen_m;
    }
    if (config_.profile) {
      efficiency_status["atmosphere.profile_source"] = config_.profile_source.empty()
                                                           ? std::string("atmospheric profile")
                                                           : config_.profile_source;
    }
    if (config_.curvature != Curvature::none) {
      efficiency_status["atmosphere.geometry"] =
          std::string(config_.curvature == Curvature::corsika ? "corsika-curved" : "spherical");
      efficiency_status["atmosphere.earth_radius_m"] = config_.earth_radius_m;
    }
    if (!config_.geometry_only) {
      // The wavelength mode of the transmission table (WavelengthInterpolation; default linear).
      efficiency_status["atmosphere.wavelength_interpolation"] = std::string(
          config_.transmission->wavelength_interpolation() == WavelengthInterpolation::linear
              ? "linear"
              : "floor");
    }
    if (config_.dispersion_reference_nm) {
      efficiency_status["atmosphere.dispersion_reference_nm"] = *config_.dispersion_reference_nm;
    }
    efficiency_status["wavelength.state"] = std::string("applied");
    efficiency_status["wavelength.source"] = std::string("density 1/lambda^2, continuous");
    efficiency_status["wavelength.range_nm"] = std::vector<double>{range.lo_nm, range.hi_nm};
    if (scattering_on()) {
      factors += ",scattering";
      efficiency_status["scattering.state"] = std::string("applied");
      efficiency_status["scattering.source"] = config_.scattering_source.empty()
                                                   ? std::string("scattering table")
                                                   : config_.scattering_source;
      efficiency_status["scattering.points_per_photon"] =
          static_cast<std::int64_t>(config_.scattering.points_per_photon);
      if (config_.scattering_primaries) {
        const auto &ids = *config_.scattering_primaries;
        efficiency_status["scattering.primaries"] =
            std::vector<double>(ids.begin(), ids.end()); // telescope ids, 1-based
      }
    }
    for (const auto &factor : config_.envelope) {
      factors += "," + factor.name;
      efficiency_status[factor.name + ".state"] = std::string("envelope");
      efficiency_status[factor.name + ".source"] = factor.source;
      efficiency_status[factor.name + ".wavelength_nm"] = factor.wavelength_nm;
      efficiency_status[factor.name + ".value"] = factor.value;
    }
    efficiency_status["factors"] = factors;
    return efficiency_status;
  }

  template <class Sink> void flush_efficiency_status(Sink &sink) {
    if (!efficiency_status_pending_) {
      return;
    }
    efficiency_status_pending_ = false;
    obdeect::stream::Record out;
    out.type = obdeect::stream::RecordType::efficiency_status;
    out.metadata = merge_status(incoming_run_status_, run_->efficiency_status);
    incoming_run_status_.clear();
    sink(out);
  }

  void open_telescope(const obdeect::stream::Metadata &metadata, obdeect::stream::Metadata &out) {
    if (!run_) {
      throw std::invalid_argument("telescope record before run record");
    }
    if (!have_event_) {
      throw std::invalid_argument("telescope record before event record");
    }
    const auto id = checked_integer(metadata, "telescope_id", obdeect::stream::max_telescope);
    const auto frame = checked_integer(metadata, "frame", 3);
    if (frame != static_cast<std::int64_t>(obdeect::stream::Frame::sky) &&
        frame != static_cast<std::int64_t>(obdeect::stream::Frame::ground)) {
      throw std::invalid_argument("the atmosphere stage needs photons in frame sky or G");
    }
    if (config_.curvature == Curvature::corsika &&
        frame == static_cast<std::int64_t>(obdeect::stream::Frame::sky)) {
      throw std::invalid_argument("curvature mode 'corsika' needs photons in frame G (CORSIKA "
                                  "CURVED photons have no sky-frame geometry)");
    }
    if (id < 1 || static_cast<std::size_t>(id) > run_->x_m.size()) {
      throw std::invalid_argument("telescope_id outside the run's telescope positions (1-based)");
    }
    const auto index = static_cast<std::size_t>(id - 1);
    telescope_.open = true;
    telescope_.id = static_cast<std::uint32_t>(id);
    telescope_.x_m = run_->x_m[index];
    telescope_.y_m = run_->y_m[index];
    telescope_.z_m = run_->z_m[index];
    telescope_.radius_m = run_->sphere ? run_->radius_m[index] : 0.0;
    telescope_.altitude_m = run_->observation_m + telescope_.z_m; // the run's observation level
    telescope_.nx = 0.0;
    telescope_.ny = 0.0;
    telescope_.nz = 1.0;
    if (config_.curvature == Curvature::spherical) {
      // Position relative to C = (0, 0, -(R + observation level)): v = (x, y, R + a), a = altitude
      // of the flat plane. The altitude is |v| - R and n = v / |v|.
      const double radius = config_.earth_radius_m + telescope_.altitude_m;
      const double rho2 = telescope_.x_m * telescope_.x_m + telescope_.y_m * telescope_.y_m;
      const double norm = std::sqrt(rho2 + radius * radius);
      telescope_.altitude_m =
          spherical_altitude_m(rho2, telescope_.altitude_m, config_.earth_radius_m);
      telescope_.nx = telescope_.x_m / norm;
      telescope_.ny = telescope_.y_m / norm;
      telescope_.nz = radius / norm;
    }
    if (config_.curvature != Curvature::none) {
      validate(curved_path(telescope_));
    }
    telescope_.pathlen_m = run_->pathlen_m[index];
    if (config_.geometry_only) {
      telescope_.depths = nullptr;
    } else {
      const auto key = std::make_pair(telescope_.altitude_m, telescope_.pathlen_m);
      auto cached = depths_cache_.find(key);
      if (cached == depths_cache_.end()) {
        cached = depths_cache_
                     .emplace(key, config_.transmission->telescope_depths(telescope_.altitude_m,
                                                                          telescope_.pathlen_m))
                     .first;
      }
      telescope_.depths = &cached->second;
    }
    // Output rows continue across array reuses of the same telescope and event.
    telescope_.out_index = out_index_[telescope_.id];
    telescope_.weight_out = 0.0;
    out["frame"] = static_cast<std::int64_t>(obdeect::stream::Frame::ground);
  }

  // Keeps the bunch weights of the output; weighted output drops rows with w = 0.
  template <class Sink> void emit(obdeect::stream::PhotonChunk &chunk, Sink &sink) {
    using obdeect::stream::Column;
    if (config_.mode == OutputMode::weighted) {
      drop_zero_weight(chunk);
    }
    if (!config_.store_p) {
      chunk[Column::p].clear();
    }
    if (chunk.rows() == 0) {
      return;
    }
    chunk.first_index = telescope_.out_index;
    telescope_.out_index += chunk.rows();
    for (const float weight : chunk[Column::w]) {
      telescope_.weight_out += static_cast<double>(weight);
      summary_.weight_out += static_cast<double>(weight);
    }
    summary_.photons_out += chunk.rows();
    obdeect::stream::Record out;
    out.type = obdeect::stream::RecordType::photons;
    out.chunk = std::move(chunk);
    sink(out);
  }

  // Removes the rows with w = 0 of a chunk, in place (the rows keep their order; no copy unless
  // a row is removed).
  static void drop_zero_weight(obdeect::stream::PhotonChunk &chunk) {
    using obdeect::stream::Column;
    const std::size_t rows = chunk.rows();
    auto &w = chunk[Column::w];
    std::size_t kept = 0;
    for (const float weight : w) {
      kept += weight > 0.0F ? 1U : 0U;
    }
    if (kept == rows) {
      return;
    }
    std::size_t j = 0;
    for (std::size_t row = 0; row < rows; ++row) {
      if (!(w[row] > 0.0F)) {
        continue;
      }
      if (j != row) {
        for (auto &column : chunk.columns) {
          if (!column.empty()) {
            column[j] = column[row];
          }
        }
      }
      ++j;
    }
    chunk.truncate(kept);
  }

  // Spherical atmosphere of the open telescope (curvature modes other than none).
  [[nodiscard]] CurvedPath curved_path(const Telescope &telescope) const {
    return CurvedPath{config_.earth_radius_m, telescope.altitude_m, telescope.nx, telescope.ny,
                      telescope.nz};
  }

  // Curved refraction table of the spherical mode at the level a [m a.s.l.], built on first use.
  // The levels lie on the grid curved_altitude_step_m (or at the profile bottom, if the grid lies
  // below it). A telescope altitude h_t between two adjacent levels is served by linear
  // interpolation in h_t between the two tables (sky_plane): the error of that interpolation is
  // second order in the 10 m spacing. Each table is looked up at its own level.
  const CurvedRefractionTables &curved_tables_at(double level_m, bool accurate) {
    const auto key = std::make_pair(level_m, accurate);
    const auto found = curved_tables_.find(key);
    if (found != curved_tables_.end()) {
      return found->second;
    }
    return curved_tables_
        .try_emplace(key, *config_.profile, config_.earth_radius_m, level_m, accurate,
                     config_.threads)
        .first->second;
  }

  // Build the more accurate dispersion tables only for a chunk containing an inclined ray.
  // The scan is outside the photon kernel and is independent of chunk size and order.
  [[nodiscard]] bool
  needs_curved_dispersion_table(const obdeect::stream::PhotonChunk &chunk) const {
    using obdeect::stream::Column;
    if (!config_.dispersion_reference_nm || config_.curvature != Curvature::spherical ||
        !chunk.has(Column::z)) {
      return false;
    }
    // A lower cached altitude can have a slightly larger zenith angle.
    const double radius = config_.earth_radius_m + telescope_.altitude_m - curved_altitude_step_m;
    const double threshold =
        std::sin(CurvedRefractionTables::accurate_dispersion_blend_start_rad());
    for (std::size_t i = 0; i < chunk.rows(); ++i) {
      const double cx = chunk[Column::cx][i], cy = chunk[Column::cy][i];
      const double cz = std::sqrt(std::max(0.0, 1.0 - cx * cx - cy * cy));
      const double impact =
          impact_parameter_m(chunk[Column::x][i], chunk[Column::y][i], chunk[Column::z][i],
                             run_->observation_m, config_.earth_radius_m, cx, cy, cz);
      if (impact >= threshold * radius) {
        return true;
      }
    }
    return false;
  }

  // Geometry of the open telescope for sky-frame photons (refraction tables built on first use).
  SkyPlane sky_plane(const obdeect::stream::PhotonChunk &chunk) {
    SkyPlane plane;
    plane.dispersion_reference_nm = config_.dispersion_reference_nm.value_or(0.0);
    plane.telescope_x_m = telescope_.x_m;
    plane.telescope_y_m = telescope_.y_m;
    plane.telescope_z_m = telescope_.z_m;
    plane.observation_m = run_->observation_m;
    if (config_.curvature == Curvature::spherical) {
      plane.earth_radius_m = config_.earth_radius_m;
      plane.telescope_altitude_m = telescope_.altitude_m; // true altitude (r_t, nu_t)
      if (config_.profile) {
        plane.profile = &*config_.profile;
        const double h = telescope_.altitude_m;
        const double bottom = config_.profile->bottom_m();
        const double top = config_.profile->top_m();
        if (!(h >= bottom && h < top)) {
          throw std::invalid_argument("telescope altitude must lie inside the atmospheric profile");
        }
        if (run_->observation_m < bottom || run_->observation_m >= top) {
          throw std::invalid_argument("observation level must lie inside the atmospheric profile");
        }
        // Levels of the two tables around h: the grid level below h (not below the profile
        // bottom) and the next level above it.
        const double lower =
            std::max(std::floor(h / curved_altitude_step_m) * curved_altitude_step_m, bottom);
        const double upper = lower + curved_altitude_step_m;
        // No upper table if h is on the lower level (zero weight) or the upper level would not lie
        // below the profile top (the table needs its level inside the profile); then the lower
        // table alone is used, which is off by the first-order change over less than one step.
        const bool blend = h > lower && upper < top;
        const bool accurate = needs_curved_dispersion_table(chunk);
        plane.curved_tables = &curved_tables_at(lower, accurate);
        plane.curved_tables_upper = blend ? &curved_tables_at(upper, accurate) : nullptr;
        plane.curved_weight = blend ? (h - lower) / (upper - lower) : 0.0;
        plane.curved_bending = !run_->refraction_applied;
      }
      return plane;
    }
    if (config_.profile) {
      plane.profile = &*config_.profile;
      if (!run_->refraction_applied) {
        const double altitude = telescope_.altitude_m;
        auto cached = tables_cache_.find(altitude);
        if (cached == tables_cache_.end()) {
          cached = tables_cache_.try_emplace(altitude, *config_.profile, altitude).first;
        }
        plane.tables = &cached->second;
      }
    }
    return plane;
  }

  // Primary photons for the scattering kernel: the input with a wavelength column (assigned later)
  // and its cz column kept (signed: upward primaries are scattered upwards to the top of the
  // table, downward ones down to the observation level). Horizontal photons (cz = 0) are not
  // scattered by the kernel. The direction is snapped to cz, so the kernel's consistency check
  // holds.
  static obdeect::stream::PhotonChunk
  scattering_primaries(const obdeect::stream::PhotonChunk &input) {
    using obdeect::stream::Column;
    // Checked before the loop below writes out.wavelength; scatter_chunk checks the rest later.
    for (const Column required :
         {Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::t, Column::w}) {
      if (!input.has(required)) {
        throw std::invalid_argument("scattering needs the columns x y z cx cy t w");
      }
    }
    obdeect::stream::PhotonChunk out = input;
    if (!out.has(Column::wavelength)) {
      out[Column::wavelength].assign(out.rows(), 0.0F);
    }
    if (input.has(Column::cz)) {
      const auto &cz = input[Column::cz];
      for (std::size_t row = 0; row < out.rows(); ++row) {
        if (!std::isfinite(cz[row])) {
          throw std::invalid_argument("cz must be finite");
        }
        obdeect::stream::snap_direction_to_cz(out[Column::cx][row], out[Column::cy][row], cz[row]);
      }
    }
    return out;
  }

  // Single scattering of the primary photons of one event, in stream order over all telescope
  // blocks. The draws are keyed on the running index of the primary photon in the event, so the
  // result does not depend on the chunking. Fills scattered_ (one list per telescope).
  void scatter_event(const std::vector<obdeect::stream::Record> &records, std::uint32_t event) {
    if (!run_) {
      throw std::invalid_argument("event record before run record");
    }
    scatter_index_ = 0;
    for (auto &list : scattered_) {
      list.clear();
    }
    for (const auto &record : records) {
      if (record.type != obdeect::stream::RecordType::photons ||
          !is_scattering_primary(record.chunk.telescope_id)) {
        continue; // not a primary: its index does not advance either
      }
      const obdeect::stream::PhotonChunk &input = record.chunk;
      if (input.frame != obdeect::stream::Frame::sky) {
        throw std::invalid_argument("single scattering needs photons in frame sky");
      }
      obdeect::stream::PhotonChunk primaries = scattering_primaries(input);
      assign_wavelengths(primaries, run_->wavelengths, config_.seed, event, input.telescope_id);
      obdeect::propagation::scatter_chunk(config_.scattering, *config_.scattering_table,
                                          *config_.transmission, sphere_, primaries, config_.seed,
                                          event, scatter_index_, scattered_);
      scatter_index_ += input.rows();
    }
  }

  // Scattered photons of the open telescope: sky to G, envelope, output. Empties their list.
  // Telescopes that are not kept get no output (their photons are computed with the others, as
  // the draws of the kernel depend on the full sphere list).
  template <class Sink> void emit_scattered(Sink &sink) {
    auto &list = scattered_[telescope_.id - 1];
    if (!selected(telescope_.id)) {
      list.clear();
      return;
    }
    const bool weighted = config_.mode == OutputMode::weighted;
    for (auto &scattered : list) {
      obdeect::stream::PhotonChunk plane_input = std::move(scattered);
      // The stored cz of the scattering gives the direction (use_cz_column).
      const DropCounts counts =
          run_->sphere
              ? sky_to_sphere(plane_input, sky_plane(plane_input), telescope_.radius_m, true)
              : sky_to_ground(plane_input, sky_plane(plane_input), true);
      summary_.scattered_dropped_by_reason += counts;
      summary_.scattered_dropped = summary_.scattered_dropped_by_reason.total();
      const double weight_in = weight_sum(plane_input);
      auto outputs = apply_envelope_chunk(
          std::move(plane_input), run_->envelope, *config_.transmission, *telescope_.depths,
          weighted, config_.seed, event_number_, telescope_.id, telescope_.out_index);
      summary_.scattered_thinned_out += weight_in - weight_sum(outputs);
      for (auto &out : outputs) {
        emit(out, sink);
      }
    }
    list.clear();
  }

  // Telescopes with scattered photons but no block in the event get one (before event_end).
  template <class Sink>
  void write_scattered_blocks(const std::vector<bool> &has_block, Sink &sink) {
    for (std::size_t i = 0; i < scattered_.size(); ++i) {
      if (has_block[i] || scattered_[i].empty() || !selected(static_cast<std::uint32_t>(i + 1))) {
        continue;
      }
      obdeect::stream::Record open;
      open.type = obdeect::stream::RecordType::telescope;
      open.metadata["telescope_id"] = static_cast<std::int64_t>(i + 1);
      open.metadata["frame"] = static_cast<std::int64_t>(obdeect::stream::Frame::sky);
      process_record(open, sink);
      emit_scattered(sink);
      obdeect::stream::Record end;
      end.type = obdeect::stream::RecordType::telescope_end;
      end.metadata["telescope_id"] = static_cast<std::int64_t>(i + 1);
      end.metadata["photons"] = 0.0;
      process_record(end, sink);
    }
  }

  // Writes a held event: scattering of its primaries, then its records. The scattered photons
  // of a telescope are written in its block, before its telescope_end.
  template <class Sink> void flush_event(Sink &sink) {
    in_event_ = false;
    std::vector<obdeect::stream::Record> records;
    records.swap(event_records_);
    const auto event = static_cast<std::uint32_t>(
        checked_integer(records.front().metadata, "event_number", max_u32));
    if (scattering_on()) {
      check_single_array_index(records);
      check_one_block_per_telescope(records);
    }
    scatter_event(records, event);
    std::vector<bool> has_block(sphere_.size(), false);
    for (const auto &record : records) {
      if (record.type == obdeect::stream::RecordType::telescope) {
        const auto id =
            checked_integer(record.metadata, "telescope_id", obdeect::stream::max_telescope);
        if (id >= 1 && static_cast<std::size_t>(id) <= has_block.size()) {
          has_block[static_cast<std::size_t>(id - 1)] = true;
        }
      }
    }
    for (auto &record : records) {
      if (record.type == obdeect::stream::RecordType::telescope_end) {
        const auto id =
            checked_integer(record.metadata, "telescope_id", obdeect::stream::max_telescope);
        if (telescope_.open && static_cast<std::int64_t>(telescope_.id) == id) {
          emit_scattered(sink);
        }
      } else if (record.type == obdeect::stream::RecordType::event_end) {
        write_scattered_blocks(has_block, sink);
      }
      process_record(record, sink);
    }
    for (auto &list : scattered_) {
      list.clear();
    }
  }

  // Photons of the open telescope. The chunk is moved into chunk_ (the record is not needed after).
  template <class Sink> void photons(obdeect::stream::PhotonChunk &input, Sink &sink) {
    if (!run_ || !telescope_.open) {
      throw std::invalid_argument("photons outside a telescope record");
    }
    if (input.telescope_id != telescope_.id) {
      throw std::invalid_argument("photons of a telescope that is not open");
    }
    if (!selected(telescope_.id)) {
      return; // not in the selection: neither counted nor written
    }
    summary_.photons_in += input.rows();
    if (input.rows() == 0) {
      return;
    }
    check_first_index(input);
    chunk_ = std::move(input);
    if (chunk_.frame == obdeect::stream::Frame::sky) {
      if (config_.dispersion_reference_nm) {
        // The bending of a dispersive run needs each photon's wavelength. Assigning it here draws
        // the same numbers as the later assignment (keyed on the row, see wavelength.hpp).
        if (!chunk_.has(obdeect::stream::Column::wavelength)) {
          chunk_[obdeect::stream::Column::wavelength].assign(chunk_.rows(), 0.0F);
        }
        assign_wavelengths(chunk_, run_->wavelengths, config_.seed, event_number_, telescope_.id);
      }
      const DropCounts counts = run_->sphere
                                    ? sky_to_sphere(chunk_, sky_plane(chunk_), telescope_.radius_m)
                                    : sky_to_ground(chunk_, sky_plane(chunk_));
      summary_.dropped_by_reason += counts;
      summary_.dropped = summary_.dropped_by_reason.total();
    } else if (run_->sphere) {
      if (!run_->sphere_input) {
        throw std::invalid_argument(
            "sphere representation: frame G photons of a plane run cannot be "
            "converted to sphere rows (use sky-frame input)");
      }
      if (!chunk_.has(obdeect::stream::Column::path)) {
        chunk_[obdeect::stream::Column::path].assign(chunk_.rows(), 0.0F); // emitted at the row
      }
      obdeect::stream::validate_sphere_rows(chunk_, telescope_.radius_m);
    } else if (run_->sphere_input) {
      throw std::invalid_argument("input run is in the sphere representation: select ground "
                                  "representation sphere");
    }

    if (config_.geometry_only) {
      for (auto &out : assign_wavelength_chunk(std::move(chunk_), run_->wavelengths, config_.seed,
                                               event_number_, telescope_.id)) {
        emit(out, sink);
      }
      return;
    }

    // Ground steps (see the header comment). Wavelength assignment, extinction and
    // envelope thinning are one kernel call; sky steps came before.
    ExtinctionSettings settings;
    settings.table = &*config_.transmission;
    settings.depths = telescope_.depths;
    settings.envelope = &run_->envelope;
    settings.wavelengths = run_->wavelengths;
    settings.seed = config_.seed;
    settings.event = event_number_;
    settings.telescope = telescope_.id;
    settings.threads = config_.threads;
    settings.pool = &pool_;
    settings.weighted = config_.mode == OutputMode::weighted;
    if (run_->sphere) {
      settings.sphere =
          SphereGeometry{run_->observation_m,
                         config_.curvature == Curvature::spherical ? config_.earth_radius_m : 0.0,
                         telescope_.x_m, telescope_.y_m, telescope_.z_m};
    } else if (config_.curvature != Curvature::none) {
      settings.curved = curved_path(telescope_);
    }
    const double weight_in = weight_sum(chunk_);
    std::vector<obdeect::stream::PhotonChunk> outputs =
        apply_extinction_chunk(std::move(chunk_), settings, telescope_.out_index);
    summary_.thinned_out += weight_in - weight_sum(outputs);
    for (auto &out : outputs) {
      emit(out, sink);
    }
  }

  PropagatorConfig config_;
  PropagationSummary summary_;
  RowPool pool_; // threads of the chunks (config_.threads), kept for the run
  std::optional<Run> run_;
  Telescope telescope_;
  bool have_event_{};
  std::map<std::uint32_t, std::uint64_t> out_index_; // per telescope, reset at each event
  std::map<std::uint32_t, std::uint64_t> in_index_;  // per telescope: first_index of the next block
  std::uint32_t event_number_{};
  bool efficiency_status_pending_{};
  obdeect::stream::Metadata incoming_run_status_; // run status records of the input, before flush
  // Caches: telescope depths per (altitude, pathlen) (the table is fixed), refraction tables per
  // telescope altitude (flat geometry, sky-frame photons), curved tables per level of the 10 m grid
  // and per accuracy of the dispersion tables (spherical geometry).
  std::map<std::pair<double, double>, TelescopeDepths> depths_cache_;
  std::map<double, RefractionTables> tables_cache_;      // by altitude [m a.s.l.]
  static constexpr double curved_altitude_step_m = 10.0; // grid of the curved table levels [m]
  std::map<std::pair<double, bool>, CurvedRefractionTables> curved_tables_; // (level, accurate)
  obdeect::stream::PhotonChunk chunk_;                                      // reused by photons()
  // Single scattering: sphere of telescope i+1, scattered photons per telescope (sky frame),
  // the records of the held event, and the running index of its primary photons.
  std::vector<ScatterTelescope> sphere_;
  std::vector<std::vector<obdeect::stream::PhotonChunk>> scattered_;
  std::vector<obdeect::stream::Record> event_records_;
  bool in_event_{};
  std::uint64_t scatter_index_{};
};

// Propagates a whole stream: reads records from in and writes them to out, with the
// atmosphere stage applied. The output stream record is written by obdeect::stream::Writer, so the
// stream record of the input is not copied. Throws on malformed input or invalid settings.
inline PropagationSummary propagate_stream(std::istream &in, std::ostream &out,
                                           PropagatorConfig config,
                                           const std::string &producer = "obdeect propagation") {
  obdeect::stream::Writer writer(out, producer);
  // Blocks whose photons are decoded (nullopt: all). Without a selection every block is written,
  // so all are decoded. With a selection and no scattering: the kept blocks. With scattering, the
  // primaries are all blocks unless scattering_primaries names them: then the kept blocks and
  // those primaries are decoded, otherwise all.
  std::optional<std::vector<std::uint32_t>> decoded;
  if (config.telescopes && !config.scattering_table) {
    decoded = *config.telescopes;
  } else if (config.telescopes && config.scattering_primaries) {
    std::vector<std::uint32_t> ids = *config.telescopes;
    ids.insert(ids.end(), config.scattering_primaries->begin(), config.scattering_primaries->end());
    sort_unique_telescope_ids(ids);
    decoded = std::move(ids);
  }
  Propagator propagator(std::move(config));
  const auto sink = [&writer](const obdeect::stream::Record &record) {
    obdeect::stream::write_record(writer, record);
  };
  obdeect::stream::Reader reader(in);
  if (decoded) {
    reader.keep_photons_of(*decoded); // photons of other telescopes are not decoded
  }
  obdeect::stream::Record record;
  while (reader.next(record)) {
    propagator.process(record, sink);
  }
  if (!reader.ended()) {
    throw std::runtime_error("truncated photon stream (no stream_end record)");
  }
  propagator.finish(sink);
  writer.finish();
  return propagator.summary();
}

} // namespace obdeect::propagation
