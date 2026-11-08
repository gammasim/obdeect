#pragma once

// Atmosphere stage of a photon stream (docs/propagation.md). Records are read with
// gstrickd::Reader, passed to Propagator::process, and written with gstrickd::Writer.
//
// Per photon chunk, in this order (Propagator::photons):
//   1. sky photons only: propagation to the telescope plane, bending, air delay, frame G (sky.hpp)
//   2. wavelength assignment for rows without a wavelength (wavelength.hpp)
//   3. atmospheric extinction p_atm(lambda, zem, airmass) (extinction.hpp)
//   4. envelope E(lambda), the product of the envelope factors, then thinning or weighting
//      (thinning.hpp)
// Optional single scattering (scattering.hpp) runs on the primaries of an event before step 3.
// Its scattered photons go through step 1 (sky to G), then step 4 with the envelope and the
// in-telescope factor exp(-od_focus) only (Propagator::emit_scattered). An event is held until its
// end, see Propagator::flush_event. Each effect of the ledger is listed in make_ledger(); a new
// effect adds a step in Propagator::photons and an entry in make_ledger().

#include "dlufd/extinction.hpp"
#include "dlufd/pipeline.hpp"
#include "dlufd/profile.hpp"
#include "dlufd/refraction.hpp"
#include "dlufd/scattering.hpp"
#include "dlufd/sky.hpp"
#include "dlufd/streams.hpp"
#include "dlufd/thinning.hpp"
#include "dlufd/timing.hpp"
#include "dlufd/transmission.hpp"
#include "dlufd/wavelength.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <istream>
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

namespace dlufd {

enum class OutputMode { thinned, weighted };

// An envelope factor f^env(lambda) of a later efficiency (ledger state `envelope`).
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
  // extinction, envelope or thinning. All rows keep their weight; the ledger records
  // `atmosphere` as pending. Needs weighted mode, no transmission table, no envelope, no
  // scattering, no store_p and no in-telescope path.
  bool geometry_only{false};
  std::string transmission_source; // ledger provenance; default "transmission table"
  double site_altitude_m{};        // observation level [m a.s.l.]
  // Optical path inside a telescope [m] (sim_telarray pathlen_tel, the od2focus part of the
  // telescope depths). default_pathlen_m applies to every telescope; telescope_pathlen_m (key:
  // telescope_id, 1-based) overrides it for one telescope. 0: no path inside the telescope.
  double default_pathlen_m{0.0};
  std::map<std::uint32_t, double> telescope_pathlen_m;
  // Optional: air delay for sky photons, and refraction (bending) unless the run says that
  // refraction is already applied.
  std::optional<AtmosphereProfile> profile;
  std::string profile_source;
  std::vector<EnvelopeFactor> envelope;
  // Wavelengths drawn from 1/lambda^2 for rows without one. Default: the run keys
  // wavelength_min_nm and wavelength_max_nm.
  std::optional<WavelengthRange> wavelength_range;
  std::uint64_t seed{};
  OutputMode mode{OutputMode::thinned};
  bool store_p{false}; // keep the acceptance column p in the output
  // Optional single scattering of the sky photons (docs/scattering.md), off by default. When
  // enabled, the table is required and the run needs telescope_r_m (sphere radii).
  ScatteringConfig scattering;
  std::optional<ScatteringTable> scattering_table;
  std::string scattering_source; // ledger provenance, e.g. the table file name
};

struct PropagationSummary {
  std::uint64_t photons_in{};  // rows of the input photon records
  std::uint64_t photons_out{}; // rows written (thinned: surviving photons, weighted: w > 0)
  std::uint64_t dropped{};     // sky photons that do not reach the telescope plane
  double weight_out{};         // sum of w of the output
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
      throw std::invalid_argument("envelope table: each line needs wavelength and value");
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

class Propagator {
public:
  explicit Propagator(PropagatorConfig config) : config_(std::move(config)) {
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
    if (config_.wavelength_range) {
      config_.wavelength_range->validate();
    }
    check_pathlen(config_.default_pathlen_m);
    for (const auto &[id, pathlen] : config_.telescope_pathlen_m) {
      if (id < 1 || id > gstrickd::max_telescope) {
        throw std::invalid_argument("telescope_pathlen_m: telescope ids are 1-based");
      }
      check_pathlen(pathlen);
    }
    if (config_.scattering.enabled != config_.scattering_table.has_value()) {
      throw std::invalid_argument(config_.scattering.enabled
                                      ? "scattering needs a scattering table"
                                      : "scattering table given, but scattering is off");
    }
    if (config_.scattering.enabled &&
        (config_.scattering.points_per_photon < 1 ||
         config_.scattering.points_per_photon > max_points_per_photon)) {
      throw std::invalid_argument("scattering needs 1 to 1024 points per photon");
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

  // Processes one input record. Records of the output go to sink(const gstrickd::Record &).
  // The ledger of a run is written before the first record that is not a ledger after it.
  // Throws std::invalid_argument if an incoming ledger shows a factor applied already.
  template <class Sink> void process(const gstrickd::Record &record, Sink &&sink) {
    using gstrickd::RecordType;
    if (!config_.scattering.enabled) {
      process_record(record, sink);
      return;
    }
    // With scattering, an event is held until its end: the scattered photons of a telescope come
    // from the primaries of all telescopes of the event.
    if (in_event_) {
      if (record.type == RecordType::event_end) {
        event_records_.push_back(record);
        flush_event(sink);
        return;
      }
      const bool ends_event = record.type == RecordType::event || record.type == RecordType::run ||
                              record.type == RecordType::run_end ||
                              record.type == RecordType::stream ||
                              record.type == RecordType::stream_end;
      if (!ends_event) {
        event_records_.push_back(record);
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

  // Writes the held event (if any) and the ledger still pending.
  template <class Sink> void finish(Sink &&sink) {
    if (in_event_) {
      flush_event(sink);
    }
    flush_ledger(sink);
  }

  [[nodiscard]] const PropagationSummary &summary() const noexcept { return summary_; }

private:
  static constexpr std::int64_t max_u32 = 0xFFFFFFFFLL;

  // One record of the stream (process without the event hold).
  template <class Sink> void process_record(const gstrickd::Record &record, Sink &sink) {
    using gstrickd::RecordType;
    if (record.type == RecordType::ledger) {
      check_incoming_ledger(record.metadata);
      sink(record);
      return;
    }
    flush_ledger(sink);
    switch (record.type) {
    case RecordType::stream:
    case RecordType::event_end:
    case RecordType::run_end:
    case RecordType::blob:
      sink(record);
      break;
    case RecordType::run: {
      gstrickd::Record out = record;
      begin_run(record.metadata, out.metadata);
      sink(out);
      ledger_pending_ = true;
      break;
    }
    case RecordType::event:
      event_number_ =
          static_cast<std::uint32_t>(checked_integer(record.metadata, "event_number", max_u32));
      have_event_ = true;
      sink(record);
      break;
    case RecordType::telescope: {
      gstrickd::Record out = record;
      open_telescope(record.metadata, out.metadata);
      sink(out);
      break;
    }
    case RecordType::photons:
      photons(record.chunk, sink);
      break;
    case RecordType::telescope_end: {
      const auto id = checked_integer(record.metadata, "telescope_id", gstrickd::max_telescope);
      if (!telescope_.open || static_cast<std::int64_t>(telescope_.id) != id) {
        throw std::invalid_argument("telescope_end of a telescope that is not open");
      }
      gstrickd::Record out = record;
      out.metadata["photons"] = telescope_.weight_out;
      telescope_.open = false;
      sink(out);
      break;
    }
    case RecordType::ledger:
    case RecordType::stream_end:
      break;
    }
  }

  struct Run {
    double observation_m{};
    bool refraction_applied{};
    std::vector<double> x_m, y_m, z_m;
    std::vector<double> pathlen_m; // path inside each telescope [m], index = telescope_id - 1
    WavelengthRange wavelengths{};
    Envelope envelope;
    gstrickd::Metadata ledger;
  };

  struct Telescope {
    bool open{};
    std::uint32_t id{};
    double altitude_m{}; // a.s.l. of the telescope plane
    double pathlen_m{};  // path inside the telescope [m]
    double x_m{}, y_m{}, z_m{};
    const TelescopeDepths *depths{}; // cached per (altitude, pathlen) (depths_cache_)
    std::uint64_t out_index{};       // output rows so far (first_index of the next chunk)
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
    if (config_.scattering.enabled || config_.scattering_table) {
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

  static bool valid_factor_name(const std::string &name) {
    if (name.empty() || name == "atmosphere" || name == "wavelength" || name == "scattering") {
      return false;
    }
    return std::all_of(name.begin(), name.end(), [](char c) {
      return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
             c == '_' || c == '-';
    });
  }

  static double number_field(const gstrickd::Metadata &metadata, const char *key) {
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

  static std::int64_t checked_integer(const gstrickd::Metadata &metadata, const char *key,
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

  static const std::vector<double> &array_field(const gstrickd::Metadata &metadata,
                                                const char *key) {
    const auto it = metadata.find(key);
    if (it == metadata.end() || !std::holds_alternative<std::vector<double>>(it->second)) {
      throw std::invalid_argument(std::string("run record lacks array key '") + key + "'");
    }
    return std::get<std::vector<double>>(it->second);
  }

  // Ledger values that other stages wrote for factors this stage applies (never apply twice).
  void check_incoming_ledger(const gstrickd::Metadata &metadata) const {
    const auto state_of = [&metadata](const std::string &name) -> std::string {
      const auto it = metadata.find(name + ".state");
      if (it == metadata.end()) {
        return "pending";
      }
      if (const auto *text = std::get_if<std::string>(&it->second)) {
        return *text;
      }
      throw std::invalid_argument("ledger state of " + name + " is not a string");
    };
    if (state_of("atmosphere") != "pending") {
      throw std::invalid_argument("atmosphere is already applied in the input ledger");
    }
    if (state_of("scattering") != "pending") {
      throw std::invalid_argument("scattering is already applied in the input ledger");
    }
    for (const auto &factor : config_.envelope) {
      if (state_of(factor.name) != "pending") {
        throw std::invalid_argument("envelope factor " + factor.name +
                                    " is already applied in the input ledger");
      }
    }
  }

  // Run keys and the ledger of this stage (docs/format.md, "Efficiency ledger").
  void begin_run(const gstrickd::Metadata &metadata, gstrickd::Metadata &out) {
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
    bool refraction_applied = false;
    if (metadata.count("refraction_applied") != 0) {
      refraction_applied = checked_integer(metadata, "refraction_applied", 1) == 1;
    } else if (config_.profile) {
      throw std::invalid_argument("run lacks refraction_applied (needed with a profile)");
    }
    WavelengthRange wavelengths{};
    if (config_.wavelength_range) {
      wavelengths = *config_.wavelength_range;
    } else {
      wavelengths = {number_field(metadata, "wavelength_min_nm"),
                     number_field(metadata, "wavelength_max_nm")};
    }
    wavelengths.validate();
    if (config_.scattering.enabled) {
      const auto &r_m = array_field(metadata, "telescope_r_m");
      if (r_m.size() != x_m.size()) {
        throw std::invalid_argument("single scattering needs telescope_r_m, one per telescope");
      }
      sphere_.clear();
      for (std::size_t i = 0; i < x_m.size(); ++i) {
        sphere_.push_back({static_cast<std::uint32_t>(i + 1), x_m[i], y_m[i], z_m[i], r_m[i]});
      }
      scattered_.assign(x_m.size(), {});
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
            x_m,
            y_m,
            z_m,
            pathlen_m,
            wavelengths,
            product_envelope(wavelengths),
            make_ledger(wavelengths, pathlen_m)};
    if (config_.profile) {
      out["refraction_applied"] = std::int64_t{1}; // sky photons are refracted from here on
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
    for (const auto &factor : config_.envelope) {
      if (factor.wavelength_nm.front() > first || factor.wavelength_nm.back() < last) {
        throw std::invalid_argument("envelope factor " + factor.name +
                                    " does not cover the wavelength range");
      }
      const Envelope table(factor.wavelength_nm, factor.value);
      for (std::size_t i = 0; i < grid.size(); ++i) {
        product[i] *= table(grid[i]);
      }
    }
    return Envelope(std::move(grid), std::move(product));
  }

  [[nodiscard]] gstrickd::Metadata make_ledger(const WavelengthRange &range,
                                               const std::vector<double> &pathlen_m) const {
    gstrickd::Metadata ledger;
    std::string factors = "atmosphere,wavelength";
    if (config_.geometry_only) {
      ledger["atmosphere.state"] = std::string("pending"); // applied by a later stage
    } else {
      ledger["atmosphere.state"] = std::string("applied");
      ledger["atmosphere.source"] = config_.transmission_source.empty()
                                        ? std::string("transmission table")
                                        : config_.transmission_source;
      ledger["atmosphere.pathlen_m"] = pathlen_m;
    }
    if (config_.profile) {
      ledger["atmosphere.profile_source"] = config_.profile_source.empty()
                                                ? std::string("atmospheric profile")
                                                : config_.profile_source;
    }
    ledger["wavelength.state"] = std::string("applied");
    ledger["wavelength.source"] = std::string("density 1/lambda^2, integer nm");
    ledger["wavelength.range_nm"] = std::vector<double>{range.lo_nm, range.hi_nm};
    if (config_.scattering.enabled) {
      factors += ",scattering";
      ledger["scattering.state"] = std::string("applied");
      ledger["scattering.source"] = config_.scattering_source.empty()
                                        ? std::string("scattering table")
                                        : config_.scattering_source;
      ledger["scattering.points_per_photon"] =
          static_cast<std::int64_t>(config_.scattering.points_per_photon);
    }
    for (const auto &factor : config_.envelope) {
      factors += "," + factor.name;
      ledger[factor.name + ".state"] = std::string("envelope");
      ledger[factor.name + ".source"] = factor.source;
      ledger[factor.name + ".wavelength_nm"] = factor.wavelength_nm;
      ledger[factor.name + ".value"] = factor.value;
    }
    ledger["factors"] = factors;
    return ledger;
  }

  template <class Sink> void flush_ledger(Sink &sink) {
    if (!ledger_pending_) {
      return;
    }
    ledger_pending_ = false;
    gstrickd::Record out;
    out.type = gstrickd::RecordType::ledger;
    out.metadata = run_->ledger;
    sink(out);
  }

  void open_telescope(const gstrickd::Metadata &metadata, gstrickd::Metadata &out) {
    if (!run_) {
      throw std::invalid_argument("telescope record before run record");
    }
    if (!have_event_) {
      throw std::invalid_argument("telescope record before event record");
    }
    const auto id = checked_integer(metadata, "telescope_id", gstrickd::max_telescope);
    const auto frame = checked_integer(metadata, "frame", 3);
    if (frame != static_cast<std::int64_t>(gstrickd::Frame::sky) &&
        frame != static_cast<std::int64_t>(gstrickd::Frame::ground)) {
      throw std::invalid_argument("the atmosphere stage needs photons in frame sky or G");
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
    telescope_.altitude_m = run_->observation_m + telescope_.z_m; // the run's observation level
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
    telescope_.out_index = 0;
    telescope_.weight_out = 0.0;
    out["frame"] = static_cast<std::int64_t>(gstrickd::Frame::ground);
  }

  // Keeps the bunch weights of the output; weighted output drops rows with w = 0.
  template <class Sink> void emit(gstrickd::PhotonChunk &chunk, Sink &sink) {
    using gstrickd::Column;
    if (config_.mode == OutputMode::weighted) {
      chunk = without_zero_weight(chunk);
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
    gstrickd::Record out;
    out.type = gstrickd::RecordType::photons;
    out.chunk = std::move(chunk);
    sink(out);
  }

  static gstrickd::PhotonChunk without_zero_weight(const gstrickd::PhotonChunk &chunk) {
    using gstrickd::Column;
    std::size_t kept = 0;
    for (const float weight : chunk[Column::w]) {
      kept += weight > 0.0F ? 1U : 0U;
    }
    if (kept == chunk.rows()) {
      return chunk;
    }
    gstrickd::PhotonChunk out;
    out.telescope_id = chunk.telescope_id;
    out.frame = chunk.frame;
    out.first_index = chunk.first_index;
    out.resize(kept, chunk.mask());
    std::size_t j = 0;
    for (std::size_t row = 0; row < chunk.rows(); ++row) {
      if (!(chunk[Column::w][row] > 0.0F)) {
        continue;
      }
      for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
        if (chunk.has(static_cast<Column>(c))) {
          out.columns[c][j] = chunk.columns[c][row];
        }
      }
      ++j;
    }
    return out;
  }

  // Geometry of the open telescope for sky photons (refraction tables built on first use).
  SkyPlane sky_plane() {
    SkyPlane plane;
    plane.telescope_x_m = telescope_.x_m;
    plane.telescope_y_m = telescope_.y_m;
    plane.telescope_z_m = telescope_.z_m;
    plane.observation_m = run_->observation_m;
    if (config_.profile) {
      plane.profile = &*config_.profile;
      if (!run_->refraction_applied) {
        if (!tables_ || tables_altitude_m_ != telescope_.altitude_m) {
          tables_.emplace(*config_.profile, telescope_.altitude_m);
          tables_altitude_m_ = telescope_.altitude_m;
        }
        plane.tables = &*tables_;
      }
    }
    return plane;
  }

  // The chunk without column `dropped`; the other columns, telescope, frame and first_index are
  // kept.
  static gstrickd::PhotonChunk without_column(const gstrickd::PhotonChunk &input,
                                              gstrickd::Column dropped) {
    using gstrickd::Column;
    gstrickd::PhotonChunk out;
    out.telescope_id = input.telescope_id;
    out.frame = input.frame;
    out.first_index = input.first_index;
    out.resize(input.rows(), input.mask() & ~gstrickd::column_bit(dropped));
    for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
      const auto column = static_cast<Column>(c);
      if (column == dropped || !input.has(column)) {
        continue;
      }
      std::copy(input.columns[c].begin(), input.columns[c].end(), out.columns[c].begin());
    }
    return out;
  }

  // Primary photons for the scattering kernel: the input without cz (the kernel takes
  // cz = -sqrt(1 - cx^2 - cy^2)), with a wavelength column. Rows of upward or horizontal photons
  // get w = 0 (they are not scattered, as sky_to_ground drops them).
  static gstrickd::PhotonChunk scattering_primaries(const gstrickd::PhotonChunk &input) {
    using gstrickd::Column;
    // Checked before the loop below writes out.w; scatter_chunk checks the rest later.
    for (const Column required :
         {Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::t, Column::w}) {
      if (!input.has(required)) {
        throw std::invalid_argument("scattering needs the columns x y z cx cy t w");
      }
    }
    gstrickd::PhotonChunk out = without_column(input, Column::cz);
    if (!out.has(Column::wavelength)) {
      out[Column::wavelength].assign(out.rows(), 0.0F);
    }
    if (input.has(Column::cz)) {
      const auto &cz = input[Column::cz];
      for (std::size_t row = 0; row < out.rows(); ++row) {
        if (!std::isfinite(cz[row])) {
          throw std::invalid_argument("cz must be finite");
        }
        if (!(cz[row] < 0.0F)) {
          out[Column::w][row] = 0.0F;
        }
      }
    }
    return out;
  }

  // Single scattering of the primary photons of one event, in stream order over all telescope
  // blocks. The draws are keyed on the running index of the primary photon in the event, so the
  // result does not depend on the chunking. Fills scattered_ (one list per telescope).
  void scatter_event(const std::vector<gstrickd::Record> &records, std::uint32_t event) {
    if (!run_) {
      throw std::invalid_argument("event record before run record");
    }
    scatter_index_ = 0;
    for (auto &list : scattered_) {
      list.clear();
    }
    for (const auto &record : records) {
      if (record.type != gstrickd::RecordType::photons) {
        continue;
      }
      const gstrickd::PhotonChunk &input = record.chunk;
      if (input.frame != gstrickd::Frame::sky) {
        throw std::invalid_argument("single scattering needs photons in frame sky");
      }
      gstrickd::PhotonChunk primaries = scattering_primaries(input);
      assign_wavelengths(primaries, run_->wavelengths, config_.seed, event, input.telescope_id);
      dlufd::scatter_chunk(config_.scattering, *config_.scattering_table, *config_.transmission,
                           sphere_, primaries, config_.seed, event, scatter_index_, scattered_);
      scatter_index_ += input.rows();
    }
  }

  // Scattered photons of the open telescope: sky to G, envelope, output. Empties their list.
  template <class Sink> void emit_scattered(Sink &sink) {
    auto &list = scattered_[telescope_.id - 1];
    const bool weighted = config_.mode == OutputMode::weighted;
    for (const auto &scattered : list) {
      // Scattered photons point downwards by construction. Their cz column (float32) can disagree
      // with cx, cy by more than the tolerance of sky_to_ground for nearly horizontal photons, so
      // their direction is taken from cx, cy, scaled inside the unit disc (float32 rounding).
      gstrickd::PhotonChunk plane_input = without_column(scattered, gstrickd::Column::cz);
      auto &cx = plane_input[gstrickd::Column::cx];
      auto &cy = plane_input[gstrickd::Column::cy];
      for (std::size_t row = 0; row < plane_input.rows(); ++row) {
        const double r2 =
            static_cast<double>(cx[row]) * cx[row] + static_cast<double>(cy[row]) * cy[row];
        if (r2 > 1.0) {
          const double scale = (1.0 - 1e-6) / std::sqrt(r2);
          cx[row] = static_cast<float>(cx[row] * scale);
          cy[row] = static_cast<float>(cy[row] * scale);
        }
      }
      sky_to_ground(plane_input, sky_plane());
      for (auto &out : apply_envelope_chunk(std::move(plane_input), run_->envelope,
                                            *telescope_.depths, weighted, config_.seed,
                                            event_number_, telescope_.id, telescope_.out_index)) {
        emit(out, sink);
      }
    }
    list.clear();
  }

  // Telescopes with scattered photons but no block in the event get one (before event_end).
  template <class Sink>
  void write_scattered_blocks(const std::vector<bool> &has_block, Sink &sink) {
    for (std::size_t i = 0; i < scattered_.size(); ++i) {
      if (has_block[i] || scattered_[i].empty()) {
        continue;
      }
      gstrickd::Record open;
      open.type = gstrickd::RecordType::telescope;
      open.metadata["telescope_id"] = static_cast<std::int64_t>(i + 1);
      open.metadata["frame"] = static_cast<std::int64_t>(gstrickd::Frame::sky);
      process_record(open, sink);
      emit_scattered(sink);
      gstrickd::Record end;
      end.type = gstrickd::RecordType::telescope_end;
      end.metadata["telescope_id"] = static_cast<std::int64_t>(i + 1);
      end.metadata["photons"] = 0.0;
      process_record(end, sink);
    }
  }

  // Writes a held event: scattering of its primaries, then its records. The scattered photons
  // of a telescope are written in its block, before its telescope_end.
  template <class Sink> void flush_event(Sink &sink) {
    in_event_ = false;
    std::vector<gstrickd::Record> records;
    records.swap(event_records_);
    const auto event = static_cast<std::uint32_t>(
        checked_integer(records.front().metadata, "event_number", max_u32));
    scatter_event(records, event);
    std::vector<bool> has_block(sphere_.size(), false);
    for (const auto &record : records) {
      if (record.type == gstrickd::RecordType::telescope) {
        const auto id = checked_integer(record.metadata, "telescope_id", gstrickd::max_telescope);
        if (id >= 1 && static_cast<std::size_t>(id) <= has_block.size()) {
          has_block[static_cast<std::size_t>(id - 1)] = true;
        }
      }
    }
    for (const auto &record : records) {
      if (record.type == gstrickd::RecordType::telescope_end) {
        const auto id = checked_integer(record.metadata, "telescope_id", gstrickd::max_telescope);
        if (telescope_.open && static_cast<std::int64_t>(telescope_.id) == id) {
          emit_scattered(sink);
        }
      } else if (record.type == gstrickd::RecordType::event_end) {
        write_scattered_blocks(has_block, sink);
      }
      process_record(record, sink);
    }
    for (auto &list : scattered_) {
      list.clear();
    }
  }

  template <class Sink> void photons(const gstrickd::PhotonChunk &input, Sink &sink) {
    if (!run_ || !telescope_.open) {
      throw std::invalid_argument("photons outside a telescope record");
    }
    if (input.telescope_id != telescope_.id) {
      throw std::invalid_argument("photons of a telescope that is not open");
    }
    summary_.photons_in += input.rows();
    if (input.rows() == 0) {
      return;
    }
    chunk_ = input;
    if (input.frame == gstrickd::Frame::sky) {
      summary_.dropped += sky_to_ground(chunk_, sky_plane());
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
    settings.weighted = config_.mode == OutputMode::weighted;
    std::vector<gstrickd::PhotonChunk> outputs =
        apply_extinction_chunk(std::move(chunk_), settings, telescope_.out_index);
    for (auto &out : outputs) {
      emit(out, sink);
    }
  }

  PropagatorConfig config_;
  PropagationSummary summary_;
  std::optional<Run> run_;
  Telescope telescope_;
  bool have_event_{};
  std::uint32_t event_number_{};
  bool ledger_pending_{};
  // Caches: telescope depths per (altitude, pathlen) (the table is fixed), refraction tables of
  // the last altitude used for sky photons.
  std::map<std::pair<double, double>, TelescopeDepths> depths_cache_;
  std::optional<RefractionTables> tables_;
  double tables_altitude_m_{};
  gstrickd::PhotonChunk chunk_; // reused by photons()
  // Single scattering: sphere of telescope i+1, scattered photons per telescope (sky frame),
  // the records of the held event, and the running index of its primary photons.
  std::vector<ScatterTelescope> sphere_;
  std::vector<std::vector<gstrickd::PhotonChunk>> scattered_;
  std::vector<gstrickd::Record> event_records_;
  bool in_event_{};
  std::uint64_t scatter_index_{};
};

// Propagates a whole stream: reads records from in and writes them to out, with the
// atmosphere stage applied. The output stream record is written by gstrickd::Writer, so the
// stream record of the input is not copied. Throws on malformed input or invalid settings.
inline PropagationSummary propagate_stream(std::istream &in, std::ostream &out,
                                           PropagatorConfig config,
                                           const std::string &producer = "dlufd") {
  gstrickd::Writer writer(out, producer);
  Propagator propagator(std::move(config));
  const auto sink = [&writer](const gstrickd::Record &record) {
    switch (record.type) {
    case gstrickd::RecordType::stream:
      break;
    case gstrickd::RecordType::photons:
      writer.write_chunk(record.chunk);
      break;
    case gstrickd::RecordType::blob:
      writer.write_blob(record.blob_eventio_type, record.blob_eventio_id, record.blob);
      break;
    default:
      writer.write_metadata(record.type, record.metadata);
      break;
    }
  };
  gstrickd::Reader reader(in);
  gstrickd::Record record;
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

} // namespace dlufd
