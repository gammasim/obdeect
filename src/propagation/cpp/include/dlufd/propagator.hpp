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
// Each effect of the ledger is listed in make_ledger(); a new effect adds a step in
// Propagator::photons and an entry in make_ledger().

#include "dlufd/extinction.hpp"
#include "dlufd/pipeline.hpp"
#include "dlufd/profile.hpp"
#include "dlufd/refraction.hpp"
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
  // Required. Its H2 must match site_altitude_m (within 25 m).
  std::optional<TransmissionTable> transmission;
  std::string transmission_source; // ledger provenance; default "transmission table"
  double site_altitude_m{};        // observation level [m a.s.l.]
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
    if (!config_.transmission) {
      throw std::invalid_argument("propagator needs a transmission table");
    }
    if (!std::isfinite(config_.site_altitude_m)) {
      throw std::invalid_argument("site altitude must be finite");
    }
    if (std::abs(config_.transmission->observation_level_m() - config_.site_altitude_m) >
        altitude_tolerance_m) {
      throw std::invalid_argument("transmission table H2 does not match the site altitude");
    }
    if (config_.wavelength_range) {
      config_.wavelength_range->validate();
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

  // Writes the ledger still pending after the last record.
  template <class Sink> void finish(Sink &&sink) { flush_ledger(sink); }

  [[nodiscard]] const PropagationSummary &summary() const noexcept { return summary_; }

private:
  static constexpr std::int64_t max_u32 = 0xFFFFFFFFLL;

  struct Run {
    double observation_m{};
    bool refraction_applied{};
    std::vector<double> x_m, y_m, z_m;
    WavelengthRange wavelengths{};
    Envelope envelope;
    gstrickd::Metadata ledger;
  };

  struct Telescope {
    bool open{};
    std::uint32_t id{};
    double altitude_m{}; // a.s.l. of the telescope plane
    double x_m{}, y_m{}, z_m{};
    const TelescopeDepths *depths{}; // cached per altitude (depths_cache_)
    std::uint64_t out_index{};       // output rows so far (first_index of the next chunk)
    double weight_out{};
  };

  static bool valid_factor_name(const std::string &name) {
    if (name.empty() || name == "atmosphere" || name == "wavelength") {
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
    Run run{observation_m,
            refraction_applied,
            x_m,
            y_m,
            z_m,
            wavelengths,
            product_envelope(wavelengths),
            make_ledger(wavelengths)};
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

  [[nodiscard]] gstrickd::Metadata make_ledger(const WavelengthRange &range) const {
    gstrickd::Metadata ledger;
    std::string factors = "atmosphere,wavelength";
    ledger["atmosphere.state"] = std::string("applied");
    ledger["atmosphere.source"] = config_.transmission_source.empty()
                                      ? std::string("transmission table")
                                      : config_.transmission_source;
    if (config_.profile) {
      ledger["atmosphere.profile_source"] = config_.profile_source.empty()
                                                ? std::string("atmospheric profile")
                                                : config_.profile_source;
    }
    ledger["wavelength.state"] = std::string("applied");
    ledger["wavelength.source"] = std::string("density 1/lambda^2, integer nm");
    ledger["wavelength.range_nm"] = std::vector<double>{range.lo_nm, range.hi_nm};
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
    auto cached = depths_cache_.find(telescope_.altitude_m);
    if (cached == depths_cache_.end()) {
      cached = depths_cache_
                   .emplace(telescope_.altitude_m,
                            config_.transmission->telescope_depths(telescope_.altitude_m, 0.0))
                   .first;
    }
    telescope_.depths = &cached->second;
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
      summary_.dropped += sky_to_ground(chunk_, plane);
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
  // Caches: telescope depths per altitude (the table is fixed), refraction tables of the
  // last altitude used for sky photons.
  std::map<double, TelescopeDepths> depths_cache_;
  std::optional<RefractionTables> tables_;
  double tables_altitude_m_{};
  gstrickd::PhotonChunk chunk_; // reused by photons()
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
