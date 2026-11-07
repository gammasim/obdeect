#pragma once

#include "dlamb/source.hpp"
#include "dlamb/vec3.hpp"

#include "gstrickd/stream.hpp"

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dlamb {

// Telescope (fiducial sphere) of a run. Positions are relative to the observation level.
struct TelescopeSphere {
  std::uint32_t id{1};
  Vec3 centre{};
  double radius_m{};
};

struct RunInfo {
  int run_number{1};
  double observation_level_m{};
  double wavelength_min_nm{};
  double wavelength_max_nm{};
  std::vector<TelescopeSphere> telescopes{};
};

struct EventInfo {
  std::uint32_t event_number{1};
  int primary_id{};
  double energy_gev{};
  double core_x_m{};
  double core_y_m{};
};

// Totals of a written stream, for diagnostics.
struct StreamSummary {
  std::uint64_t photons{}; // Rows written.
  double weight{};         // Sum of photon weights.
  double expected{};       // Expected sum of weights from the program's own formula.
};

// Seed of source `index` of a run with `seed` (splitmix64 mix): every source gets its own
// Philox key, so two sources never share draws.
[[nodiscard]] inline std::uint64_t source_seed(std::uint64_t seed, std::uint64_t index) {
  std::uint64_t z = seed + 0x9E3779B97F4A7C15ULL * (index + 1);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

// Writes one run of photons in the sky frame: stream and run records, then per event the
// event record, one block per telescope (telescope record, photon chunks, telescope_end)
// and event_end. Photons are buffered and written in chunks of at most 2^16 rows, with
// first_index counting the photons of each telescope and event.
class RunWriter {
public:
  RunWriter(std::ostream &out, const std::string &producer, const RunInfo &run)
      : writer_(out, producer) {
    gstrickd::Metadata meta;
    meta["run_number"] = static_cast<std::int64_t>(run.run_number);
    meta["observation_level_m"] = run.observation_level_m;
    meta["wavelength_min_nm"] = run.wavelength_min_nm;
    meta["wavelength_max_nm"] = run.wavelength_max_nm;
    meta["refraction_applied"] = std::int64_t{0};
    meta["n_telescopes"] = static_cast<std::int64_t>(run.telescopes.size());
    std::vector<double> x, y, z, r;
    for (const auto &t : run.telescopes) {
      x.push_back(t.centre.x);
      y.push_back(t.centre.y);
      z.push_back(t.centre.z);
      r.push_back(t.radius_m);
    }
    meta["telescope_x_m"] = x;
    meta["telescope_y_m"] = y;
    meta["telescope_z_m"] = z;
    meta["telescope_r_m"] = r;
    writer_.write_metadata(gstrickd::RecordType::run, meta);
  }

  void begin_event(const EventInfo &event) {
    gstrickd::Metadata meta;
    meta["event_number"] = static_cast<std::int64_t>(event.event_number);
    meta["primary_id"] = static_cast<std::int64_t>(event.primary_id);
    meta["energy_gev"] = event.energy_gev;
    meta["core_x_m"] = std::vector<double>{event.core_x_m};
    meta["core_y_m"] = std::vector<double>{event.core_y_m};
    writer_.write_metadata(gstrickd::RecordType::event, meta);
  }

  void begin_telescope(const TelescopeSphere &telescope, gstrickd::Frame frame) {
    if (in_telescope_)
      throw std::logic_error("RunWriter: telescope block already open");
    in_telescope_ = true;
    telescope_id_ = telescope.id;
    frame_ = frame;
    next_index_ = 0;
    weight_sum_ = 0.0;
    pending_.clear();
    gstrickd::Metadata meta;
    meta["telescope_id"] = static_cast<std::int64_t>(telescope.id);
    meta["frame"] = static_cast<std::int64_t>(frame);
    meta["x_m"] = telescope.centre.x;
    meta["y_m"] = telescope.centre.y;
    meta["z_m"] = telescope.centre.z;
    meta["r_m"] = telescope.radius_m;
    writer_.write_metadata(gstrickd::RecordType::telescope, meta);
  }

  // Buffers one photon of the open telescope block.
  void push(const Photon &photon) {
    if (!in_telescope_)
      throw std::logic_error("RunWriter: push outside a telescope block");
    pending_.push_back(photon);
    weight_sum_ += photon.weight;
    weight_total_ += photon.weight;
    ++photons_;
    if (pending_.size() == gstrickd::max_chunk_rows)
      flush();
  }

  void end_telescope() {
    if (!in_telescope_)
      throw std::logic_error("RunWriter: no telescope block open");
    flush();
    gstrickd::Metadata meta;
    meta["telescope_id"] = static_cast<std::int64_t>(telescope_id_);
    meta["photons"] = weight_sum_;
    writer_.write_metadata(gstrickd::RecordType::telescope_end, meta);
    in_telescope_ = false;
  }

  void end_event() { writer_.write_metadata(gstrickd::RecordType::event_end, {}); }

  // Writes run_end and stream_end.
  void finish() {
    writer_.write_metadata(gstrickd::RecordType::run_end, {});
    writer_.finish();
  }

  [[nodiscard]] std::uint64_t photons() const { return photons_; }
  [[nodiscard]] double weight() const { return weight_total_; }

private:
  void flush() {
    if (pending_.empty())
      return;
    gstrickd::PhotonChunk chunk;
    chunk.telescope_id = telescope_id_;
    chunk.frame = frame_;
    chunk.first_index = next_index_;
    chunk.resize(pending_.size(), source_columns);
    auto &x = chunk[gstrickd::Column::x], &y = chunk[gstrickd::Column::y],
         &z = chunk[gstrickd::Column::z], &cx = chunk[gstrickd::Column::cx],
         &cy = chunk[gstrickd::Column::cy], &cz = chunk[gstrickd::Column::cz],
         &t = chunk[gstrickd::Column::t], &w = chunk[gstrickd::Column::w],
         &wl = chunk[gstrickd::Column::wavelength];
    for (std::size_t i = 0; i < pending_.size(); ++i) {
      const Photon &p = pending_[i];
      x[i] = static_cast<float>(p.position.x);
      y[i] = static_cast<float>(p.position.y);
      z[i] = static_cast<float>(p.position.z);
      cx[i] = static_cast<float>(p.direction.x);
      cy[i] = static_cast<float>(p.direction.y);
      cz[i] = static_cast<float>(p.direction.z);
      t[i] = static_cast<float>(p.time_ns);
      w[i] = static_cast<float>(p.weight);
      wl[i] = static_cast<float>(p.wavelength_nm);
    }
    writer_.write_chunk(chunk);
    next_index_ += pending_.size();
    pending_.clear();
  }

  gstrickd::Writer writer_;
  std::vector<Photon> pending_;
  std::uint32_t telescope_id_{};
  gstrickd::Frame frame_{gstrickd::Frame::sky};
  std::uint64_t next_index_{};
  std::uint64_t photons_{};
  double weight_sum_{};
  double weight_total_{};
  bool in_telescope_{false};
};

} // namespace dlamb
