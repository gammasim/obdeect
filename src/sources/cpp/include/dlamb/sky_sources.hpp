#pragma once

// Sky-frame run/plan layer shared by the point-source, beam-scan, laser-beam and night-sky
// tools (docs/spec/sky-sources.md).

#include "dlamb/background.hpp"
#include "dlamb/source.hpp"
#include "dlamb/vec3.hpp"

#include "gstrickd/photons.hpp"
#include "gstrickd/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

namespace dlamb {

// Every sky emission is generated once with this telescope key; all telescopes get the same
// photons.
inline constexpr std::uint32_t shared_telescope_key = 0;

inline Vec3 direction_from_altitude(double altitude_rad, double azimuth_rad) {
  return {std::cos(altitude_rad) * std::cos(azimuth_rad),
          std::cos(altitude_rad) * std::sin(azimuth_rad), std::sin(altitude_rad)};
}

struct SkyTelescope {
  std::uint32_t id{1};
  Vec3 centre{};
  double radius_m{};
};

inline double fiducial_super_radius(const std::vector<SkyTelescope> &telescopes) {
  double r = 0.0;
  for (const auto &t : telescopes)
    r = std::max(r, norm(t.centre) + t.radius_m);
  return r;
}

struct SkyRun {
  std::int64_t run_number{1};
  double observation_level_m{2150.0};
  std::int64_t atmosphere_id{-1};
  double wavelength_min_nm{};
  double wavelength_max_nm{};
  std::vector<SkyTelescope> telescopes;
  gstrickd::Metadata extra; // added to the run record without overwriting
};

// A Source emitted with a bunch size and a time offset, for all telescopes alike.
struct SourceEmission {
  Source source;
  std::uint64_t seed{1};
  double bunch{1.0};
  double time_offset_ns{};

  [[nodiscard]] std::uint64_t count(std::uint32_t event) const {
    return photon_count(source, seed, event, shared_telescope_key);
  }

  void fill(std::uint32_t event, std::uint64_t first, std::size_t rows,
            gstrickd::PhotonChunk &chunk) const {
    emit(source, seed, event, shared_telescope_key, first, rows, chunk);
    const float bunch_f = static_cast<float>(bunch);
    const float offset_f = static_cast<float>(time_offset_ns);
    auto &w = chunk[gstrickd::Column::w];
    auto &t = chunk[gstrickd::Column::t];
    for (std::size_t i = 0; i < rows; ++i) {
      w[i] = w[i] * bunch_f;
      t[i] = t[i] + offset_f;
    }
  }
};

// Background photons with weight = bunch size.
struct BackgroundEmission {
  BackgroundSource source;
  std::uint64_t seed{1};
  double bunch{1.0};
  double mean_bunches{};

  [[nodiscard]] std::uint64_t count(std::uint32_t event) const {
    return background_count(mean_bunches, seed, event, shared_telescope_key);
  }

  void fill(std::uint32_t event, std::uint64_t first, std::size_t rows,
            gstrickd::PhotonChunk &chunk) const {
    emit_background(source, seed, event, shared_telescope_key, first, rows, bunch, chunk);
  }
};

template <class Emission> struct SkyPlan {
  SkyRun run;
  std::uint64_t events_per_case{1};
  std::vector<std::vector<Emission>> cases;
};

inline gstrickd::Metadata sky_run_metadata(const SkyRun &run) {
  gstrickd::Metadata m;
  m["run_number"] = run.run_number;
  m["observation_level_m"] = run.observation_level_m;
  m["atmosphere_id"] = run.atmosphere_id;
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
  if (run.wavelength_max_nm > 0.0) {
    m["wavelength_min_nm"] = run.wavelength_min_nm;
    m["wavelength_max_nm"] = run.wavelength_max_nm;
  }
  for (const auto &[key, value] : run.extra)
    m.emplace(key, value);
  return m;
}

// Writes the stream of a plan: each case is written events_per_case times, events numbered 1, 2...
template <class Emission>
void write_plan(std::ostream &out, const std::string &producer, const SkyPlan<Emission> &plan) {
  using gstrickd::Column;
  gstrickd::Writer writer(out, producer);
  writer.write_metadata(gstrickd::RecordType::run, sky_run_metadata(plan.run));
  gstrickd::PhotonChunk chunk;
  std::uint64_t event_number = 0;
  std::vector<std::uint64_t> counts;
  for (const auto &emissions : plan.cases) {
    for (std::uint64_t repeat = 0; repeat < plan.events_per_case; ++repeat) {
      ++event_number;
      const auto event = static_cast<std::uint32_t>(event_number);
      writer.write_metadata(gstrickd::RecordType::event,
                            {{"event_number", static_cast<std::int64_t>(event_number)}});
      counts.clear();
      for (const auto &emission : emissions)
        counts.push_back(emission.count(event));
      for (const auto &tel : plan.run.telescopes) {
        gstrickd::Metadata tm;
        tm["telescope_id"] = static_cast<std::int64_t>(tel.id);
        tm["frame"] = static_cast<std::int64_t>(gstrickd::Frame::sky);
        tm["x_m"] = tel.centre.x;
        tm["y_m"] = tel.centre.y;
        tm["z_m"] = tel.centre.z;
        tm["r_m"] = tel.radius_m;
        writer.write_metadata(gstrickd::RecordType::telescope, tm);
        std::uint64_t written = 0;
        double weight = 0.0;
        for (std::size_t k = 0; k < emissions.size(); ++k) {
          for (std::uint64_t done = 0; done < counts[k];) {
            const auto rows = static_cast<std::size_t>(
                std::min<std::uint64_t>(gstrickd::max_chunk_rows, counts[k] - done));
            emissions[k].fill(event, done, rows, chunk);
            chunk.telescope_id = tel.id;
            chunk.first_index = written;
            auto &zem = chunk[Column::zem];
            zem.resize(rows);
            for (std::size_t i = 0; i < rows; ++i) {
              zem[i] = static_cast<float>(plan.run.observation_level_m +
                                          static_cast<double>(chunk[Column::z][i]));
              weight += static_cast<double>(chunk[Column::w][i]);
            }
            writer.write_chunk(chunk);
            written += rows;
            done += rows;
          }
        }
        writer.write_metadata(
            gstrickd::RecordType::telescope_end,
            {{"telescope_id", static_cast<std::int64_t>(tel.id)}, {"photons", weight}});
      }
      writer.write_metadata(gstrickd::RecordType::event_end, {});
    }
  }
  writer.write_metadata(gstrickd::RecordType::run_end, {});
  writer.finish();
}

} // namespace dlamb
