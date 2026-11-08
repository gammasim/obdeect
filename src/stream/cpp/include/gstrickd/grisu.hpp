#pragma once

// Photon stream to GrISU/GrOptics photon text input (ps2grisu). Header-only,
// standard library only.
//
// Input: ground-frame photons (frame G, CORSIKA axes: x north, y west) with
// the zem and wavelength columns, atmosphere applied upstream (efficiency status factor
// `atmosphere` in state `applied`). Output: a header block, then per event
// `S` lines and per photon `P` lines (docs/usage.md, "ps2grisu").
//
// Frames: x_g = -y_c, y_g = -x_c for positions. Azimuths: az_g = (3 pi / 2 -
// az_c) mod 2 pi, with az_c = atan2(cy, cx) for photons. The direction cosines
// written are sin(ze) cos(az_g) and sin(ze) sin(az_g) with sin(ze) =
// sqrt(cx^2 + cy^2). Values below 1e-8 in magnitude are written as 0.

#include "gstrickd/simtel.hpp"
#include "gstrickd/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <istream>
#include <numbers>
#include <optional>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace gstrickd {

// A GrISU telescope location (TLLOC line): id and position [m], x east, y north.
struct GrisuTelescope {
  std::uint32_t id{};
  double x{};
  double y{};
};

struct GrisuOptions {
  double qe{1.0};                   // global QE factor written on the `R` line, not applied
  bool expand_bunches{};            // write round(w) copies of a bunch with w > 1
  bool event_details{};             // write a `C` line after each `S` line
  std::string producer{"ps2grisu"}; // name and version for the header
  std::optional<std::vector<GrisuTelescope>> telescopes; // TLLOC map; none: id = stream id
  std::ostream *warnings{};
};

namespace grisu_detail {

inline constexpr double two_pi = 2.0 * std::numbers::pi;
inline constexpr double tiny = 1e-8;
inline constexpr double match_distance_m = 0.5;
inline constexpr double unit_weight_tolerance = 1e-6;

inline double zero_tiny(double value) { return std::abs(value) < tiny ? 0.0 : value; }

inline double wrap_two_pi(double angle) {
  double wrapped = std::fmod(angle, two_pi);
  if (wrapped < 0.0) {
    wrapped += two_pi;
  }
  return wrapped;
}

} // namespace grisu_detail

// CORSIKA ground position (x north, y west) to GrISU (x east, y north).
[[nodiscard]] inline std::pair<double, double> corsika_to_grisu_position(double x, double y) {
  return {-y, -x};
}

// GrISU position to the CORSIKA ground frame (the same mapping is its own inverse up to
// the exchange of the roles of the two axes).
[[nodiscard]] inline std::pair<double, double> grisu_to_corsika_position(double x, double y) {
  return {-y, -x};
}

// GrISU direction cosines (dcx, dcy) from a CORSIKA zenith and azimuth [rad].
[[nodiscard]] inline std::pair<double, double> grisu_direction_from_angles(double zenith,
                                                                           double azimuth) {
  const double azimuth_g = grisu_detail::wrap_two_pi(1.5 * std::numbers::pi - azimuth);
  const double sin_ze = std::sin(zenith);
  return {sin_ze * std::cos(azimuth_g), sin_ze * std::sin(azimuth_g)};
}

// GrISU direction cosines of a photon from its CORSIKA ground-frame (cx, cy).
[[nodiscard]] inline std::pair<double, double> corsika_to_grisu_direction(double cx, double cy) {
  const double azimuth = std::atan2(cy, cx);
  const double azimuth_g = grisu_detail::wrap_two_pi(1.5 * std::numbers::pi - azimuth);
  const double sin_ze = std::sqrt(cx * cx + cy * cy);
  return {sin_ze * std::cos(azimuth_g), sin_ze * std::sin(azimuth_g)};
}

// Inverse of corsika_to_grisu_direction: CORSIKA (cx, cy) from GrISU (dcx, dcy).
[[nodiscard]] inline std::pair<double, double> grisu_to_corsika_direction(double dcx, double dcy) {
  const double azimuth_g = std::atan2(dcy, dcx);
  const double azimuth = grisu_detail::wrap_two_pi(1.5 * std::numbers::pi - azimuth_g);
  const double sin_ze = std::sqrt(dcx * dcx + dcy * dcy);
  return {sin_ze * std::cos(azimuth), sin_ze * std::sin(azimuth)};
}

// Reads `* TLLOC id x y z offset px py` lines (x east, y north [m]). Other lines are
// comments. Throws on a malformed TLLOC line or a repeated id.
inline std::vector<GrisuTelescope> read_tlloc(std::istream &in) {
  std::vector<GrisuTelescope> telescopes;
  std::string line;
  std::size_t number = 0;
  while (std::getline(in, line)) {
    ++number;
    std::istringstream words(line);
    std::string first;
    words >> first;
    std::string name;
    if (first == "*") {
      words >> name;
    } else if (first.size() > 1 && first[0] == '*') {
      name = first.substr(1);
    }
    if (name != "TLLOC") {
      continue;
    }
    long long id = 0;
    double values[6] = {};
    words >> id;
    for (double &value : values) {
      words >> value;
    }
    if (words.fail() || id < 1 || id > 0xFFFFFFFFLL || !std::isfinite(values[0]) ||
        !std::isfinite(values[1])) {
      throw std::runtime_error("TLLOC line " + std::to_string(number) +
                               " is not `* TLLOC id x y z offset px py`");
    }
    for (const GrisuTelescope &known : telescopes) {
      if (known.id == id) {
        throw std::runtime_error("TLLOC telescope " + std::to_string(id) + " listed twice");
      }
    }
    telescopes.push_back({static_cast<std::uint32_t>(id), values[0], values[1]});
  }
  return telescopes;
}

// Stream records to GrISU photon text.
class PsToGrisu {
public:
  PsToGrisu(std::ostream &out, GrisuOptions options) : out_{out}, options_{std::move(options)} {
    if (!std::isfinite(options_.qe)) {
      throw std::invalid_argument("the QE factor must be finite");
    }
  }

  void write(const Record &record) {
    switch (record.type) {
    case RecordType::stream:
      if (const auto it = record.metadata.find("producer"); it != record.metadata.end()) {
        if (const auto *name = std::get_if<std::string>(&it->second)) {
          stream_producer_ = *name;
        }
      }
      return;
    case RecordType::run:
      end_run();
      run_metadata_ = record.metadata;
      run_open_ = true;
      return;
    case RecordType::efficiency_status:
      efficiency_status(record.metadata);
      return;
    case RecordType::event:
      begin_event(record.metadata);
      return;
    case RecordType::telescope:
      begin_telescope(record.metadata);
      return;
    case RecordType::photons:
      photons(record.chunk);
      return;
    case RecordType::telescope_end:
      telescope_open_ = false;
      return;
    case RecordType::event_end:
      event_open_ = false;
      return;
    case RecordType::run_end:
      end_run();
      return;
    case RecordType::blob:
    case RecordType::stream_end:
      return; // CORSIKA blocks are not part of the GrISU input
    }
  }

  void finish() {
    end_run();
    if (options_.warnings != nullptr && !skipped_.empty()) {
      *options_.warnings << "ps2grisu: skipped " << skipped_.size()
                         << " telescope(s) without a TLLOC match within "
                         << grisu_detail::match_distance_m << " m\n";
    }
    out_.flush();
  }

private:
  static double number(const Metadata &metadata, const char *key, const char *what) {
    const auto it = metadata.find(key);
    if (it == metadata.end()) {
      throw std::runtime_error(std::string("ps2grisu: stream lacks '") + key + "' (" + what + ")");
    }
    if (const auto *integer = std::get_if<std::int64_t>(&it->second)) {
      return static_cast<double>(*integer);
    }
    if (const auto *real = std::get_if<double>(&it->second)) {
      return *real;
    }
    throw std::runtime_error(std::string("ps2grisu: stream key '") + key + "' is not a number");
  }

  static std::optional<double> optional_number(const Metadata &metadata, const char *key) {
    const auto it = metadata.find(key);
    if (it == metadata.end()) {
      return std::nullopt;
    }
    return number(metadata, key, "");
  }

  static const std::vector<double> *array(const Metadata &metadata, const char *key) {
    const auto it = metadata.find(key);
    return it == metadata.end() ? nullptr : std::get_if<std::vector<double>>(&it->second);
  }

  void put(const char *format, auto... values) {
    char buffer[512];
    const int n = std::snprintf(buffer, sizeof buffer, format, values...);
    if (n < 0 || static_cast<std::size_t>(n) >= sizeof buffer) {
      throw std::runtime_error("ps2grisu: output line too long");
    }
    out_.write(buffer, n);
  }

  void efficiency_status(const Metadata &metadata) {
    if (!run_open_) {
      throw std::runtime_error("ps2grisu: efficiency status record outside a run");
    }
    if (telescope_open_) {
      if (efficiency_status_checked_) {
        throw std::runtime_error("ps2grisu: efficiency status for telescope " +
                                 std::to_string(telescope_id_) + " comes after its photons");
      }
      telescope_efficiency_statuses_.push_back(metadata);
    } else {
      run_efficiency_statuses_.push_back(metadata);
    }
  }

  // Header block; the primary lines are taken from the first event of the run.
  void write_header(const Metadata *event) {
    if (header_done_) {
      return;
    }
    header_done_ = true;
    const double height = number(run_metadata_, "observation_level_m", "observation height");
    out_ << "* HEADF  <-- Start of header flag\n";
    out_ << "* Produced by " << options_.producer;
    if (!stream_producer_.empty()) {
      out_ << " from a photon stream written by " << stream_producer_;
    }
    out_ << "\n";
    const auto line = [&](const char *label, std::optional<double> value, double scale = 1.0) {
      out_ << "* " << label << ": ";
      if (value) {
        char buffer[64];
        std::snprintf(buffer, sizeof buffer, "%.9g", *value * scale);
        out_ << buffer;
      } else {
        out_ << "n/a";
      }
      out_ << "\n";
    };
    line("CORSIKA run number", optional_number(run_metadata_, "run_number"));
    {
      const auto low = optional_number(run_metadata_, "energy_min_gev");
      const auto high = optional_number(run_metadata_, "energy_max_gev");
      out_ << "* Energy range [TeV]: ";
      if (low && high) {
        char buffer[96];
        std::snprintf(buffer, sizeof buffer, "%.9g %.9g", *low * 1e-3, *high * 1e-3);
        out_ << buffer;
      } else {
        out_ << "n/a";
      }
      out_ << "\n";
    }
    line("Spectral slope", optional_number(run_metadata_, "spectral_index"));
    line("Primary id", event != nullptr ? optional_number(*event, "primary_id") : std::nullopt);
    const double degrees = 180.0 / std::numbers::pi;
    line("Primary zenith [deg]",
         event != nullptr ? optional_number(*event, "zenith_rad") : std::nullopt, degrees);
    line("Primary azimuth [deg, CORSIKA]",
         event != nullptr ? optional_number(*event, "azimuth_rad") : std::nullopt, degrees);
    line("Observation height [m]", height);
    out_ << "* DATAF  <-- end of header flag\n";
    put("R %.9g\n", options_.qe);
    put("H %.9g\n", height);
  }

  void begin_event(const Metadata &metadata) {
    if (!run_open_) {
      throw std::runtime_error("ps2grisu: event record outside a run");
    }
    write_header(&metadata);
    event_ = metadata;
    event_open_ = true;
    last_array_ = -1;
    event_energy_tev_ = number(metadata, "energy_gev", "primary energy") * 1e-3;
    event_zenith_ = number(metadata, "zenith_rad", "primary zenith");
    event_azimuth_ = number(metadata, "azimuth_rad", "primary azimuth");
    event_first_interaction_ = number(metadata, "first_interaction_m", "first interaction height");
    event_number_ = static_cast<std::int64_t>(number(metadata, "event_number", "event number"));
  }

  void begin_telescope(const Metadata &metadata) {
    if (!event_open_) {
      throw std::runtime_error("ps2grisu: telescope record outside an event");
    }
    telescope_open_ = true;
    telescope_efficiency_statuses_.clear();
    efficiency_status_checked_ = false;
    mapped_ = false;
    const auto id_it = metadata.find("telescope_id");
    const auto *id = id_it == metadata.end() ? nullptr : std::get_if<std::int64_t>(&id_it->second);
    if (id == nullptr || *id < 1 || *id > 0xFFFFFFFFLL) {
      throw std::runtime_error("ps2grisu: telescope record needs an integer telescope_id >= 1");
    }
    telescope_id_ = static_cast<std::uint32_t>(*id);
    // Position in the CORSIKA ground frame [m].
    double x = 0.0;
    double y = 0.0;
    if (const auto px = optional_number(metadata, "x_m"), py = optional_number(metadata, "y_m");
        px && py) {
      x = *px;
      y = *py;
    } else {
      const auto *xs = array(run_metadata_, "telescope_x_m");
      const auto *ys = array(run_metadata_, "telescope_y_m");
      if (xs == nullptr || ys == nullptr || telescope_id_ == 0 || telescope_id_ > xs->size() ||
          telescope_id_ > ys->size()) {
        throw std::runtime_error("ps2grisu: no position for telescope " +
                                 std::to_string(telescope_id_));
      }
      x = (*xs)[telescope_id_ - 1];
      y = (*ys)[telescope_id_ - 1];
    }
    telescope_x_ = x;
    telescope_y_ = y;
    grisu_id_ = telescope_id_;
    if (options_.telescopes) {
      const auto [gx, gy] = corsika_to_grisu_position(x, y);
      double best = grisu_detail::match_distance_m;
      std::optional<std::uint32_t> match;
      for (const GrisuTelescope &candidate : *options_.telescopes) {
        const double distance = std::hypot(candidate.x - gx, candidate.y - gy);
        if (distance <= best) {
          best = distance;
          match = candidate.id;
        }
      }
      if (!match) {
        skipped_.insert(telescope_id_);
        return;
      }
      grisu_id_ = *match;
    }
    mapped_ = true;
    const auto array_it = metadata.find("array_index");
    const std::int64_t array_index =
        array_it != metadata.end() && std::holds_alternative<std::int64_t>(array_it->second)
            ? std::get<std::int64_t>(array_it->second)
            : 0;
    if (array_index != last_array_) {
      last_array_ = array_index;
      write_shower(array_index);
    }
  }

  // `S` line (and optional `C` line) of the shower on one array position.
  void write_shower(std::int64_t array_index) {
    double core_x = 0.0;
    double core_y = 0.0;
    const auto *xs = array(event_, "core_x_m");
    const auto *ys = array(event_, "core_y_m");
    if (xs != nullptr && ys != nullptr) {
      if (array_index < 0 || static_cast<std::size_t>(array_index) >= xs->size() ||
          static_cast<std::size_t>(array_index) >= ys->size()) {
        throw std::runtime_error("ps2grisu: no core position for array index " +
                                 std::to_string(array_index));
      }
      core_x = (*xs)[static_cast<std::size_t>(array_index)];
      core_y = (*ys)[static_cast<std::size_t>(array_index)];
    }
    const auto [gx, gy] = corsika_to_grisu_position(core_x, core_y);
    const auto [dcos, dsin] = grisu_direction_from_angles(event_zenith_, event_azimuth_);
    put("S %.9g %.9g %.9g %.9g %.9g %.9g -1 -1 -1\n", event_energy_tev_,
        grisu_detail::zero_tiny(gx), grisu_detail::zero_tiny(gy), grisu_detail::zero_tiny(dcos),
        grisu_detail::zero_tiny(dsin), event_first_interaction_);
    if (options_.event_details) {
      put("C %.9g -1 %lld\n", event_first_interaction_, static_cast<long long>(event_number_));
    }
  }

  // The atmosphere must be applied upstream: state of the factor from the run efficiency status
  // records, overridden per factor by the efficiency status records of the telescope.
  void check_atmosphere() {
    std::string state;
    for (const auto *efficiency_statuses :
         {&run_efficiency_statuses_, &telescope_efficiency_statuses_}) {
      for (const Metadata &efficiency_status : *efficiency_statuses) {
        const auto states = simtel_detail::efficiency_status_states(efficiency_status);
        if (const auto it = states.find("atmosphere"); it != states.end()) {
          state = it->second;
        }
      }
    }
    if (state != "applied") {
      throw std::runtime_error(
          "ps2grisu: efficiency status of telescope " + std::to_string(telescope_id_) +
          (state.empty() ? " lacks the factor 'atmosphere'" : ": factor 'atmosphere' is " + state) +
          "; the atmosphere must be applied upstream (state applied)");
    }
  }

  void photons(const PhotonChunk &chunk) {
    if (!telescope_open_) {
      throw std::runtime_error("ps2grisu: photons outside a telescope record");
    }
    if (!efficiency_status_checked_) {
      check_atmosphere();
      efficiency_status_checked_ = true;
    }
    if (!mapped_ || chunk.rows() == 0) {
      return;
    }
    if (chunk.frame != Frame::ground) {
      throw std::runtime_error("ps2grisu: photons of telescope " + std::to_string(telescope_id_) +
                               " are not in the ground frame");
    }
    if (!chunk.has(Column::zem)) {
      throw std::runtime_error("ps2grisu: photons need the zem column");
    }
    if (!chunk.has(Column::wavelength)) {
      throw std::runtime_error("ps2grisu: wavelength is not assigned (no wavelength column)");
    }
    const auto &x = chunk[Column::x];
    const auto &y = chunk[Column::y];
    const auto &cx = chunk[Column::cx];
    const auto &cy = chunk[Column::cy];
    const auto &t = chunk[Column::t];
    const auto &w = chunk[Column::w];
    const auto &zem = chunk[Column::zem];
    const auto &lambda = chunk[Column::wavelength];
    std::string text;
    char buffer[256];
    for (std::size_t i = 0; i < chunk.rows(); ++i) {
      if (!(lambda[i] >= 0.5F) || !std::isfinite(lambda[i])) {
        throw std::runtime_error("ps2grisu: wavelength is not assigned (value " +
                                 std::to_string(lambda[i]) + " nm)");
      }
      long long copies = 1;
      if (!std::isfinite(w[i])) {
        throw std::runtime_error("ps2grisu: photon weight is not finite");
      }
      if (!(std::abs(w[i] - 1.0F) <= grisu_detail::unit_weight_tolerance)) {
        if (!(w[i] > 1.0F)) {
          throw std::runtime_error("ps2grisu: photon weight " + std::to_string(w[i]) +
                                   " is below 1; thin the photons upstream");
        }
        if (!options_.expand_bunches) {
          throw std::runtime_error("ps2grisu: photon weight " + std::to_string(w[i]) +
                                   " is above 1; use --expand-bunches to write copies");
        }
        if (w[i] > 1e6F) {
          throw std::runtime_error("ps2grisu: photon weight " + std::to_string(w[i]) +
                                   " is too large to expand");
        }
        copies = std::llround(static_cast<double>(w[i]));
      }
      const auto [gx, gy] = corsika_to_grisu_position(telescope_x_ + static_cast<double>(x[i]),
                                                      telescope_y_ + static_cast<double>(y[i]));
      const auto [dcx, dcy] = corsika_to_grisu_direction(cx[i], cy[i]);
      const int n = std::snprintf(
          buffer, sizeof buffer, "P %+.9g %+.9g %+.9g %+.9g %+.9g %+.9g %lld 3 %u\n",
          grisu_detail::zero_tiny(gx), grisu_detail::zero_tiny(gy), grisu_detail::zero_tiny(dcx),
          grisu_detail::zero_tiny(dcy), static_cast<double>(zem[i]), static_cast<double>(t[i]),
          std::llround(static_cast<double>(lambda[i])), grisu_id_);
      if (n < 0 || static_cast<std::size_t>(n) >= sizeof buffer) {
        throw std::runtime_error("ps2grisu: output line too long");
      }
      for (long long k = 0; k < copies; ++k) {
        text.append(buffer, static_cast<std::size_t>(n));
      }
    }
    out_.write(text.data(), static_cast<std::streamsize>(text.size()));
  }

  void end_run() {
    if (!run_open_) {
      return;
    }
    write_header(nullptr);
    run_open_ = false;
    header_done_ = false;
    run_efficiency_statuses_.clear();
    run_metadata_.clear();
  }

  std::ostream &out_;
  GrisuOptions options_;
  std::string stream_producer_;
  Metadata run_metadata_;
  Metadata event_;
  std::vector<Metadata> run_efficiency_statuses_;
  std::vector<Metadata> telescope_efficiency_statuses_;
  std::set<std::uint32_t> skipped_;
  bool run_open_{};
  bool header_done_{};
  bool event_open_{};
  bool telescope_open_{};
  bool efficiency_status_checked_{};
  bool mapped_{};
  std::int64_t last_array_{-1};
  std::int64_t event_number_{};
  double event_energy_tev_{};
  double event_zenith_{};
  double event_azimuth_{};
  double event_first_interaction_{};
  std::uint32_t telescope_id_{};
  std::uint32_t grisu_id_{};
  double telescope_x_{};
  double telescope_y_{};
};

// Photon stream from `in` to GrISU photon text on `out`. Throws on a truncated stream or
// on input the format cannot represent.
inline void stream_to_grisu(std::istream &in, std::ostream &out, GrisuOptions options) {
  PsToGrisu converter{out, std::move(options)};
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
