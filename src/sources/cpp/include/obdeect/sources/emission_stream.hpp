#pragma once

// Run writer for sources that push photons one by one (muon track): gstrickd records with
// buffered chunks (docs/spec/stream-writer.md, writer C).

#include "obdeect/sources/source.hpp"
#include "obdeect/sources/vec3.hpp"

#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::sources {

// Seed of emission i of a base seed (splitmix64 of base + golden * (i + 1)).
inline std::uint64_t source_seed(std::uint64_t base, std::uint64_t i) {
  std::uint64_t z = base + 0x9E3779B97F4A7C15ULL * (i + 1);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

// Sphere around a telescope: the region photons are generated for.
struct TelescopeSphere {
  std::uint32_t id{1};
  Vec3 centre{};
  double radius_m{};
};

struct RunInfo {
  std::int64_t run_number{1};
  double observation_level_m{2150.0};
  double wavelength_min_nm{};
  double wavelength_max_nm{};
  std::vector<TelescopeSphere> telescopes;
};

struct EventInfo {
  std::int64_t event_number{1};
  std::int64_t primary_id{};
  double energy_gev{};
  double core_x_m{};
  double core_y_m{};
};

struct StreamSummary {
  std::uint64_t photons{}; // rows written
  double weight{};         // sum of the weights
  double expected{};       // expected sum of weights (set by the caller)
};

class RunWriter {
public:
  RunWriter(std::ostream &out, const std::string &producer, const RunInfo &run)
      : writer_(out, producer), observation_level_m_(run.observation_level_m) {
    obdeect::stream::Metadata m;
    m["run_number"] = run.run_number;
    m["observation_level_m"] = run.observation_level_m;
    m["wavelength_min_nm"] = run.wavelength_min_nm;
    m["wavelength_max_nm"] = run.wavelength_max_nm;
    m["refraction_applied"] = std::int64_t{0};
    m["n_telescopes"] = static_cast<std::int64_t>(run.telescopes.size());
    std::vector<double> x, y, z, r;
    for (const auto &t : run.telescopes) {
      x.push_back(t.centre.x);
      y.push_back(t.centre.y);
      z.push_back(t.centre.z);
      r.push_back(t.radius_m);
    }
    m["telescope_x_m"] = x;
    m["telescope_y_m"] = y;
    m["telescope_z_m"] = z;
    m["telescope_r_m"] = r;
    writer_.write_metadata(obdeect::stream::RecordType::run, m);
  }

  void begin_event(const EventInfo &event) {
    if (in_event_)
      throw std::logic_error("RunWriter: nested event");
    obdeect::stream::Metadata m;
    m["event_number"] = event.event_number;
    m["primary_id"] = event.primary_id;
    m["energy_gev"] = event.energy_gev;
    m["core_x_m"] = std::vector<double>{event.core_x_m};
    m["core_y_m"] = std::vector<double>{event.core_y_m};
    writer_.write_metadata(obdeect::stream::RecordType::event, m);
    in_event_ = true;
  }

  void begin_telescope(const TelescopeSphere &t) {
    if (!in_event_ || in_telescope_)
      throw std::logic_error("RunWriter: telescope block outside an event or nested");
    obdeect::stream::Metadata m;
    m["telescope_id"] = static_cast<std::int64_t>(t.id);
    m["frame"] = static_cast<std::int64_t>(obdeect::stream::Frame::sky);
    m["x_m"] = t.centre.x;
    m["y_m"] = t.centre.y;
    m["z_m"] = t.centre.z;
    m["r_m"] = t.radius_m;
    writer_.write_metadata(obdeect::stream::RecordType::telescope, m);
    telescope_id_ = t.id;
    in_telescope_ = true;
    flushed_ = 0;
    block_weight_ = 0.0;
    buffer_.clear();
  }

  void push(const Photon &p) {
    if (!in_telescope_)
      throw std::logic_error("RunWriter: photon outside a telescope block");
    buffer_.push_back(p);
    block_weight_ += p.weight;
    ++summary_.photons;
    summary_.weight += p.weight;
    if (buffer_.size() == obdeect::stream::max_chunk_rows)
      flush();
  }

  void end_telescope() {
    if (!in_telescope_)
      throw std::logic_error("RunWriter: no telescope block to end");
    flush();
    obdeect::stream::Metadata m;
    m["telescope_id"] = static_cast<std::int64_t>(telescope_id_);
    m["photons"] = block_weight_;
    writer_.write_metadata(obdeect::stream::RecordType::telescope_end, m);
    in_telescope_ = false;
  }

  void end_event() {
    if (!in_event_ || in_telescope_)
      throw std::logic_error("RunWriter: no event to end or telescope block still open");
    writer_.write_metadata(obdeect::stream::RecordType::event_end, {});
    in_event_ = false;
  }

  void finish() {
    if (in_event_)
      throw std::logic_error("RunWriter: event still open");
    writer_.write_metadata(obdeect::stream::RecordType::run_end, {});
    writer_.finish();
  }

  [[nodiscard]] const StreamSummary &summary() const { return summary_; }

private:
  void flush() {
    if (buffer_.empty())
      return;
    using obdeect::stream::Column;
    const std::size_t n = buffer_.size();
    chunk_.resize(n, source_columns | obdeect::stream::column_bit(Column::zem));
    chunk_.telescope_id = telescope_id_;
    chunk_.frame = obdeect::stream::Frame::sky;
    chunk_.first_index = flushed_;
    for (std::size_t i = 0; i < n; ++i) {
      const Photon &p = buffer_[i];
      chunk_[Column::x][i] = static_cast<float>(p.position.x);
      chunk_[Column::y][i] = static_cast<float>(p.position.y);
      chunk_[Column::z][i] = static_cast<float>(p.position.z);
      chunk_[Column::cx][i] = static_cast<float>(p.direction.x);
      chunk_[Column::cy][i] = static_cast<float>(p.direction.y);
      chunk_[Column::cz][i] = static_cast<float>(p.direction.z);
      chunk_[Column::t][i] = static_cast<float>(p.time_ns);
      chunk_[Column::w][i] = static_cast<float>(p.weight);
      chunk_[Column::wavelength][i] = static_cast<float>(p.wavelength_nm);
      chunk_[Column::zem][i] =
          static_cast<float>(observation_level_m_ + static_cast<double>(chunk_[Column::z][i]));
    }
    writer_.write_chunk(chunk_);
    flushed_ += n;
    buffer_.clear();
  }

  obdeect::stream::Writer writer_;
  double observation_level_m_;
  bool in_event_{false}, in_telescope_{false};
  std::uint32_t telescope_id_{};
  std::uint64_t flushed_{};
  double block_weight_{};
  std::vector<Photon> buffer_;
  obdeect::stream::PhotonChunk chunk_;
  StreamSummary summary_;
};

} // namespace obdeect::sources
