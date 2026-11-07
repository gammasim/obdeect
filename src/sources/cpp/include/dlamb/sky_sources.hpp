#pragma once

// Shared parts of the sky-level light-source programs (ports of sim_telarray
// LightEmission): run and event layout of the photon stream, emissions of core
// sources, and vector helpers with the rotation sign convention of LightEmission
// (SpaceVect::Rotate_y, Rotate_z).

#include "dlamb/source.hpp"
#include "dlamb/streams.hpp"
#include "dlamb/vec3.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace dlamb {

// Telescope key of all sky-frame photons: one emission, identical for every telescope
// of the run (LightEmission bunches are written once, telescopes only select them).
inline constexpr std::uint32_t shared_telescope_key = 0;

// A telescope fiducial sphere in the sky frame (metres; z above the observation level).
struct SkyTelescope {
  std::uint32_t id{1};
  Vec3 centre{};
  double radius_m{};
};

// Run keys (docs/format.md, "Run keys") and program-specific keys in `extra`.
struct SkyRun {
  std::int64_t run_number{1};
  double observation_level_m{2150.0};
  std::int64_t atmosphere_id{-1};
  double wavelength_min_nm{0.0};
  double wavelength_max_nm{0.0};
  std::vector<SkyTelescope> telescopes{};
  gstrickd::Metadata extra{};
};

// Rotations with the sign convention of LightEmission (SpaceVect::Rotate_y/Rotate_z).
[[nodiscard]] inline Vec3 rotate_y(Vec3 v, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return {v.x * c + v.z * s, v.y, v.z * c - v.x * s};
}
[[nodiscard]] inline Vec3 rotate_z(Vec3 v, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return {v.x * c - v.y * s, v.x * s + v.y * c, v.z};
}

// Unit vector of DirectionVect(altitude, azimuth) in LightEmission.
[[nodiscard]] inline Vec3 direction_from_altitude(double altitude_rad, double azimuth_rad) {
  return {std::cos(altitude_rad) * std::cos(azimuth_rad),
          std::cos(altitude_rad) * std::sin(azimuth_rad), std::sin(altitude_rad)};
}

// Largest reach of the telescopes from the origin: max(|centre| + radius), as Run::Rmax().
[[nodiscard]] inline double fiducial_super_radius(const std::vector<SkyTelescope> &telescopes) {
  double r = 0.0;
  for (const auto &t : telescopes)
    r = std::max(r, norm(t.centre) + t.radius_m);
  return r;
}

// Wavelength range of a spectrum for the run keys (both 0 when not known).
[[nodiscard]] inline std::pair<double, double> spectrum_range(const Spectrum &spectrum) {
  return std::visit(
      [](const auto &s) -> std::pair<double, double> {
        using S = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<S, Fixed>) {
          return {s.value, s.value};
        } else if constexpr (std::is_same_v<S, UniformRange> || std::is_same_v<S, InverseSquare>) {
          return {s.low, s.high};
        } else {
          return {s.density.lower(), s.density.upper()};
        }
      },
      spectrum);
}

// Emission of a core source for one event. Each row is one bunch of `bunch` photons
// (LightEmission bunch; gstrickd column w = photons in the bunch). `time_offset_ns`
// is added to the emission time.
struct SourceEmission {
  Source source{};
  std::uint64_t seed{};
  double bunch{1.0};
  double time_offset_ns{0.0};

  [[nodiscard]] std::uint64_t count(std::uint32_t event) const {
    return photon_count(source, seed, event, shared_telescope_key);
  }
  void fill(std::uint32_t event, std::uint64_t first, std::size_t rows,
            gstrickd::PhotonChunk &chunk) const {
    emit(source, seed, event, shared_telescope_key, first, rows, chunk);
    auto &w = chunk[gstrickd::Column::w];
    auto &t = chunk[gstrickd::Column::t];
    for (std::size_t i = 0; i < rows; ++i) {
      w[i] *= static_cast<float>(bunch);
      t[i] += static_cast<float>(time_offset_ns);
    }
  }
};

// Writes the records of one stream: stream, run (with telescope keys), events with
// telescope blocks, run_end and stream_end (docs/format.md, nesting).
class SkyStream {
public:
  SkyStream(std::ostream &out, const std::string &producer) : writer_(out, producer) {}

  void begin_run(const SkyRun &run) {
    telescopes_ = run.telescopes;
    gstrickd::Metadata m;
    m["run_number"] = run.run_number;
    m["observation_level_m"] = run.observation_level_m;
    m["atmosphere_id"] = run.atmosphere_id;
    m["refraction_applied"] = std::int64_t{0};
    m["n_telescopes"] = static_cast<std::int64_t>(telescopes_.size());
    std::vector<double> x, y, z, r;
    for (const auto &t : telescopes_) {
      x.push_back(t.centre.x);
      y.push_back(t.centre.y);
      z.push_back(t.centre.z);
      r.push_back(t.radius_m);
    }
    m["telescope_x_m"] = x;
    m["telescope_y_m"] = y;
    m["telescope_z_m"] = z;
    m["telescope_r_m"] = r;
    if (run.wavelength_max_nm > 0.0) {
      m["wavelength_min_nm"] = run.wavelength_min_nm;
      m["wavelength_max_nm"] = run.wavelength_max_nm;
    }
    m.insert(run.extra.begin(), run.extra.end());
    writer_.write_metadata(gstrickd::RecordType::run, m);
  }

  // One event: `emissions` are generated for every telescope; photon indices run over
  // all emissions of the telescope and event.
  template <class Emission>
  void write_event(std::uint32_t event, const std::vector<Emission> &emissions) {
    gstrickd::Metadata e;
    e["event_number"] = static_cast<std::int64_t>(event);
    writer_.write_metadata(gstrickd::RecordType::event, e);
    std::vector<std::uint64_t> counts;
    counts.reserve(emissions.size());
    for (const auto &em : emissions)
      counts.push_back(em.count(event));
    for (const auto &tel : telescopes_) {
      gstrickd::Metadata t;
      t["telescope_id"] = static_cast<std::int64_t>(tel.id);
      t["frame"] = static_cast<std::int64_t>(gstrickd::Frame::sky);
      t["x_m"] = tel.centre.x;
      t["y_m"] = tel.centre.y;
      t["z_m"] = tel.centre.z;
      t["r_m"] = tel.radius_m;
      writer_.write_metadata(gstrickd::RecordType::telescope, t);
      double weight_sum = 0.0;
      std::uint64_t index = 0;
      gstrickd::PhotonChunk chunk;
      for (std::size_t k = 0; k < emissions.size(); ++k) {
        for (std::uint64_t first = 0; first < counts[k]; first += gstrickd::max_chunk_rows) {
          const auto rows = static_cast<std::size_t>(
              std::min<std::uint64_t>(gstrickd::max_chunk_rows, counts[k] - first));
          emissions[k].fill(event, first, rows, chunk);
          chunk.telescope_id = tel.id;
          chunk.frame = gstrickd::Frame::sky;
          chunk.first_index = index + first;
          for (const float w : chunk[gstrickd::Column::w])
            weight_sum += w;
          writer_.write_chunk(chunk);
        }
        index += counts[k];
      }
      gstrickd::Metadata end;
      end["telescope_id"] = static_cast<std::int64_t>(tel.id);
      end["photons"] = weight_sum;
      writer_.write_metadata(gstrickd::RecordType::telescope_end, end);
    }
    writer_.write_metadata(gstrickd::RecordType::event_end, {});
  }

  void end_run() {
    writer_.write_metadata(gstrickd::RecordType::run_end, {});
    writer_.finish();
  }

private:
  gstrickd::Writer writer_;
  std::vector<SkyTelescope> telescopes_;
};

// A complete program output: run keys, and the emissions of each case. Each case is
// written `events_per_case` times (event numbers continue across cases).
template <class Emission> struct SkyPlan {
  SkyRun run{};
  std::uint32_t events_per_case{1};
  std::vector<std::vector<Emission>> cases{};
};

template <class Emission>
void write_plan(std::ostream &out, const std::string &producer, const SkyPlan<Emission> &plan) {
  SkyStream stream(out, producer);
  stream.begin_run(plan.run);
  std::uint32_t event = 0;
  for (const auto &emissions : plan.cases)
    for (std::uint32_t j = 0; j < plan.events_per_case; ++j)
      stream.write_event(++event, emissions);
  stream.end_run();
}

} // namespace dlamb
