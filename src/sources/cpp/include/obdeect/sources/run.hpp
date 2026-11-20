#pragma once

// A run of light-source events written as an obdeect photon stream (docs/spec/run.md): one run
// record, then per event one telescope block per telescope holding the photons of all emissions.

#include "obdeect/sources/emit.hpp"
#include "obdeect/sources/source.hpp"

#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace obdeect::sources {

struct Telescope {
  std::uint32_t id{1};
  Vec3 position{};
  double radius_m{};
};

// A source with a bunch size (photons per row: candidates = photons / bunch, weight x bunch),
// a time offset, and optionally a straight-line move of every photon to the plane z = 0.
struct Emission {
  Source source;
  double bunch_size{1.0};
  double time_offset_ns{};
  bool to_plane_z0{};
};

// Emissions written `repeat` times as events numbered on from the previous ones.
struct Event {
  std::uint64_t repeat{1};
  std::vector<Emission> emissions;
  obdeect::stream::Metadata metadata; // added to each event record
};

struct Run {
  std::int64_t number{1};
  std::uint64_t seed{1};
  double observation_level_m{};
  std::int64_t atmosphere_id{-1};
  obdeect::stream::Frame frame{obdeect::stream::Frame::sky};
  std::vector<Telescope> telescopes;
  std::vector<Event> events;
  obdeect::stream::Metadata metadata; // added to the run record (writer keys take precedence)
  std::string producer{"obdeect sources"};
};

struct RunSummary {
  std::uint64_t events{};
  std::uint64_t rows{}; // rows written, all telescopes
  double weight{};      // sum of weights, all telescopes
  double expected{};    // expected sum of weights, all telescopes
};

// Sphere that holds every telescope sphere: centred on the mean telescope position.
inline SphereTarget fiducial_sphere(const std::vector<Telescope> &telescopes) {
  SphereTarget sphere;
  if (telescopes.empty())
    return sphere;
  for (const auto &t : telescopes)
    sphere.centre = sphere.centre + t.position;
  sphere.centre = sphere.centre * (1.0 / static_cast<double>(telescopes.size()));
  for (const auto &t : telescopes)
    sphere.radius_m = std::max(sphere.radius_m, norm(t.position - sphere.centre) + t.radius_m);
  return sphere;
}

// Target "auto" of a source, from the fiducial sphere (centre C, radius R): no target for a point
// within 1.02 R of C, the sphere for an isotropic point, none for a fixed angle, otherwise a disc
// of radius R about C in the plane z = C.z. A line with a fixed angle or isotropic emission gets
// the sphere; other lines and geometries no target.
inline Target resolve_auto_target(const Source &s, const std::vector<Telescope> &telescopes) {
  const SphereTarget sphere = fiducial_sphere(telescopes);
  if (!(sphere.radius_m > 0.0))
    throw std::invalid_argument("an automatic target needs at least one telescope");
  if (const auto *l = std::get_if<Line>(&s.geometry)) {
    if (std::holds_alternative<FixedAngle>(l->angular) ||
        std::holds_alternative<Isotropic>(l->angular))
      return sphere;
    return std::monostate{};
  }
  const auto *p = std::get_if<Point>(&s.geometry);
  if (!p || norm(p->position - sphere.centre) <= 1.02 * sphere.radius_m)
    return std::monostate{};
  if (std::holds_alternative<Isotropic>(p->angular))
    return sphere;
  if (std::holds_alternative<FixedAngle>(p->angular))
    return std::monostate{};
  if (p->position.z == sphere.centre.z)
    throw std::invalid_argument("automatic target: a non-isotropic point in the plane of the "
                                "telescopes outside them is not supported");
  return DiscTarget{sphere.centre, {0.0, 0.0, 1.0}, sphere.radius_m};
}

namespace detail {

inline void validate_run(const Run &run) {
  require(run.frame == obdeect::stream::Frame::sky ||
              run.frame == obdeect::stream::Frame::telescope ||
              run.frame == obdeect::stream::Frame::camera,
          "run: frame must be sky, telescope or camera");
  require(std::isfinite(run.observation_level_m), "run: observation level must be finite");
  require(!run.telescopes.empty(), "run: needs at least one telescope");
  std::uint64_t event_count = 0;
  for (const auto &event : run.events) {
    require(event.repeat <= std::numeric_limits<std::uint32_t>::max() - event_count,
            "run: total event count exceeds the Philox event key range");
    event_count += event.repeat;
  }
  for (std::size_t i = 0; i < run.telescopes.size(); ++i) {
    const auto &t = run.telescopes[i];
    require(within_limits(t.position) && positive_within_limits(t.radius_m),
            "run: telescope position or radius out of range");
    for (std::size_t k = 0; k < i; ++k)
      require(run.telescopes[k].id != t.id, "run: telescope ids must be distinct");
  }
  for (const auto &event : run.events)
    for (const auto &e : event.emissions) {
      require(std::isfinite(e.bunch_size) && e.bunch_size >= 1.0,
              "emission: bunch size must be finite and >= 1");
      require(within_limits(e.time_offset_ns), "emission: time offset out of range");
    }
}

// Rows of one emission chunk: bunch weight, time offset, move to z = 0 and zem (sky frame).
inline void finish_rows(const Emission &e, const Run &run, obdeect::stream::PhotonChunk &chunk) {
  using obdeect::stream::Column;
  const std::size_t rows = chunk.rows();
  auto &x = chunk[Column::x];
  auto &y = chunk[Column::y];
  auto &z = chunk[Column::z];
  auto &t = chunk[Column::t];
  auto &w = chunk[Column::w];
  for (std::size_t i = 0; i < rows; ++i) {
    double time = static_cast<double>(t[i]) + e.time_offset_ns;
    if (e.to_plane_z0) {
      const double cz = chunk[Column::cz][i];
      if (!(cz < 0.0) || !(z[i] >= 0.0F))
        throw std::invalid_argument("a photon does not travel down to the plane z = 0");
      const double path = static_cast<double>(z[i]) / -cz;
      x[i] = static_cast<float>(x[i] + path * chunk[Column::cx][i]);
      y[i] = static_cast<float>(y[i] + path * chunk[Column::cy][i]);
      z[i] = 0.0F;
      time += path / speed_of_light_m_per_ns;
    }
    t[i] = static_cast<float>(time);
    w[i] = static_cast<float>(static_cast<double>(w[i]) * e.bunch_size);
  }
  if (run.frame == obdeect::stream::Frame::sky) {
    auto &zem = chunk[Column::zem];
    zem.resize(rows);
    for (std::size_t i = 0; i < rows; ++i)
      zem[i] = static_cast<float>(run.observation_level_m + static_cast<double>(z[i]));
  }
}

} // namespace detail

// Emissions of a run with automatic targets resolved, bunch sizes applied and seeds assigned:
// emission k of the run (in event order) uses source_seed(run.seed, k).
struct PlannedEmission {
  Emission emission;
  std::uint64_t seed{};
};

inline std::vector<std::vector<PlannedEmission>> plan(const Run &run) {
  detail::validate_run(run);
  std::vector<std::vector<PlannedEmission>> planned;
  std::uint64_t k = 0;
  for (const auto &event : run.events) {
    auto &list = planned.emplace_back();
    for (const auto &e : event.emissions) {
      PlannedEmission p{e, source_seed(run.seed, k++)};
      if (std::holds_alternative<AutoTarget>(p.emission.source.target))
        p.emission.source.target = resolve_auto_target(p.emission.source, run.telescopes);
      p.emission.source.photons = e.source.photons / e.bunch_size;
      validate(p.emission.source);
      list.push_back(std::move(p));
    }
  }
  return planned;
}

// Bytes per generated photon in a cached chunk (ten float32 columns including zem).
inline constexpr std::uint64_t cached_bytes_per_photon = 10 * sizeof(float);
// The photons of one event are the same for every telescope: they are cached up to this size and
// generated again per telescope above it. The output is the same either way.
inline constexpr std::uint64_t max_cached_event_bytes = std::uint64_t{1} << 30;

// Expected sum of weights of the whole run (all events and telescopes).
inline double expected_weight(const Run &run) {
  const auto planned = plan(run);
  double sum = 0.0;
  for (std::size_t i = 0; i < planned.size(); ++i)
    for (const auto &p : planned[i])
      sum += expected_weight(p.emission.source) * p.emission.bunch_size *
             static_cast<double>(run.events[i].repeat);
  return sum * static_cast<double>(run.telescopes.size());
}

// Writes the run. Photons are keyed on (emission seed, event number, telescope 0), so all
// telescopes get the same photons.
inline RunSummary write_run(std::ostream &out, const Run &run,
                            std::uint64_t max_cached_bytes = max_cached_event_bytes) {
  using obdeect::stream::Column;
  using obdeect::stream::RecordType;
  const auto planned = plan(run);
  constexpr std::uint32_t key = 0;

  obdeect::stream::Metadata m;
  m["run_number"] = run.number;
  m["seed"] = static_cast<std::int64_t>(run.seed);
  m["observation_level_m"] = run.observation_level_m;
  m["atmosphere_id"] = run.atmosphere_id;
  m["refraction_applied"] = std::int64_t{0};
  m["curvature_applied"] = std::int64_t{0};
  m["n_telescopes"] = static_cast<std::int64_t>(run.telescopes.size());
  std::vector<double> tx, ty, tz, tr;
  for (const auto &t : run.telescopes) {
    tx.push_back(t.position.x);
    ty.push_back(t.position.y);
    tz.push_back(t.position.z);
    tr.push_back(t.radius_m);
  }
  m["telescope_x_m"] = tx;
  m["telescope_y_m"] = ty;
  m["telescope_z_m"] = tz;
  m["telescope_r_m"] = tr;
  bool any = false;
  double wl_min = 0.0, wl_max = 0.0;
  for (const auto &list : planned)
    for (const auto &p : list) {
      const auto [lo, hi] = spectrum_range(p.emission.source.spectrum);
      wl_min = any ? std::min(wl_min, lo) : lo;
      wl_max = any ? std::max(wl_max, hi) : hi;
      any = true;
    }
  if (any) {
    m["wavelength_min_nm"] = wl_min;
    m["wavelength_max_nm"] = wl_max;
  }
  for (const auto &[name, value] : run.metadata)
    m.emplace(name, value);

  obdeect::stream::Writer writer(out, run.producer);
  writer.write_metadata(RecordType::run, m);
  RunSummary summary;
  summary.expected = expected_weight(run);
  obdeect::stream::PhotonChunk chunk;
  std::vector<obdeect::stream::PhotonChunk> cached;
  std::vector<std::uint64_t> counts;
  std::uint64_t event_number = 0;
  for (std::size_t i = 0; i < planned.size(); ++i) {
    const auto &emissions = planned[i];
    for (std::uint64_t repeat = 0; repeat < run.events[i].repeat; ++repeat) {
      ++event_number;
      ++summary.events;
      const auto event = static_cast<std::uint32_t>(event_number);
      obdeect::stream::Metadata em{{"event_number", static_cast<std::int64_t>(event_number)}};
      for (const auto &[name, value] : run.events[i].metadata)
        em.emplace(name, value);
      writer.write_metadata(RecordType::event, em);
      counts.clear();
      std::uint64_t total = 0;
      for (const auto &p : emissions) {
        counts.push_back(photon_count(p.emission.source, p.seed, event, key));
        total += counts.back();
      }
      // Calls f(chunk) for every chunk of the event, in emission and photon order.
      const auto generate = [&](auto &&f) {
        for (std::size_t k = 0; k < emissions.size(); ++k) {
          const auto &p = emissions[k];
          for (std::uint64_t done = 0; done < counts[k];) {
            const auto rows = static_cast<std::size_t>(
                std::min<std::uint64_t>(obdeect::stream::max_chunk_rows, counts[k] - done));
            emit(p.emission.source, p.seed, event, key, done, rows, chunk);
            detail::finish_rows(p.emission, run, chunk);
            done += rows;
            if (chunk.rows() > 0)
              f(chunk);
          }
        }
      };
      const bool use_cache =
          run.telescopes.size() > 1 && total <= max_cached_bytes / cached_bytes_per_photon;
      cached.clear();
      if (use_cache)
        generate([&](const obdeect::stream::PhotonChunk &c) { cached.push_back(c); });
      for (const auto &tel : run.telescopes) {
        obdeect::stream::Metadata tm;
        tm["telescope_id"] = static_cast<std::int64_t>(tel.id);
        tm["frame"] = static_cast<std::int64_t>(run.frame);
        tm["x_m"] = tel.position.x;
        tm["y_m"] = tel.position.y;
        tm["z_m"] = tel.position.z;
        tm["r_m"] = tel.radius_m;
        writer.write_metadata(RecordType::telescope, tm);
        std::uint64_t written = 0;
        double weight = 0.0;
        const auto write_block = [&](obdeect::stream::PhotonChunk &block) {
          block.telescope_id = tel.id;
          block.frame = run.frame;
          block.first_index = written;
          for (const float w : block[Column::w])
            weight += static_cast<double>(w);
          writer.write_chunk(block);
          written += block.rows();
        };
        if (use_cache)
          for (auto &block : cached)
            write_block(block);
        else
          generate(write_block);
        writer.write_metadata(
            RecordType::telescope_end,
            {{"telescope_id", static_cast<std::int64_t>(tel.id)}, {"photons", weight}});
        summary.rows += written;
        summary.weight += weight;
      }
      writer.write_metadata(RecordType::event_end, {});
    }
  }
  writer.write_metadata(RecordType::run_end, {});
  writer.finish();
  return summary;
}

} // namespace obdeect::sources
