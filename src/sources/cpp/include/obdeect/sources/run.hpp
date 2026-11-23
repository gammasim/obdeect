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

// Emissions written `repeat` times as events numbered on from the previous ones. In a ground run
// (frame G, sphere representation) the positions and directions of an emission are relative to
// the centre of the telescope sphere, and to_plane_z0 is not allowed.
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
  double expected{};    // expected sum of weights of all emitted photons, all telescopes
  // Photons of to_plane_z0 emissions that move away from the plane z = 0 and are dropped. Counted
  // once per event (the photons are the same for every telescope); not in rows or weight. The
  // run_end record holds the same count.
  std::uint64_t not_on_plane{};
  // Photons of a ground run whose path misses the telescope sphere (dropped, like not_on_plane).
  std::uint64_t misses_sphere{};
};

// Photons dropped by finish_rows (counted in the run summary).
struct Drops {
  std::uint64_t not_on_plane{};
  std::uint64_t misses_sphere{};
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

// The sphere that the automatic target of a run aims at: the fiducial sphere of the telescopes in
// sky, T and C runs; in a ground run the telescope sphere at the origin of the G frame.
inline SphereTarget target_sphere(const Run &run) {
  if (run.frame == obdeect::stream::Frame::ground && !run.telescopes.empty())
    return SphereTarget{Vec3{}, run.telescopes.front().radius_m};
  return fiducial_sphere(run.telescopes);
}

// Target "auto" of a source, from a sphere (centre C, radius R): no target for a point within
// 1.02 R of C, the sphere for an isotropic point, none for a fixed angle, otherwise a disc of
// radius R about C in the plane z = C.z. A line with a fixed angle or isotropic emission gets the
// sphere; other lines and geometries no target.
inline Target resolve_auto_target(const Source &s, const SphereTarget &sphere) {
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

// Automatic target from the telescopes of a sky, T or C run (the fiducial sphere).
inline Target resolve_auto_target(const Source &s, const std::vector<Telescope> &telescopes) {
  const SphereTarget sphere = fiducial_sphere(telescopes);
  if (!(sphere.radius_m > 0.0))
    throw std::invalid_argument("an automatic target needs at least one telescope");
  return resolve_auto_target(s, sphere);
}

namespace detail {

inline void validate_run(const Run &run) {
  require(run.frame == obdeect::stream::Frame::sky ||
              run.frame == obdeect::stream::Frame::telescope ||
              run.frame == obdeect::stream::Frame::camera ||
              run.frame == obdeect::stream::Frame::ground,
          "run: frame must be sky, telescope, camera or ground");
  const bool ground = run.frame == obdeect::stream::Frame::ground;
  require(std::isfinite(run.observation_level_m), "run: observation level must be finite");
  require(!run.telescopes.empty(), "run: needs at least one telescope");
  std::uint64_t event_count = 0;
  for (const auto &event : run.events) {
    require(event.repeat <= std::numeric_limits<std::uint32_t>::max() - event_count,
            "run: total event count exceeds the Philox event key range");
    event_count += event.repeat;
  }
  require(!ground || run.telescopes.size() == 1,
          "run: a ground run (frame ground) needs exactly one telescope: its sphere is the origin "
          "of the frame");
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
      require(!(ground && e.to_plane_z0),
              "emission: to_plane_z0 is not used in a ground run (its rows are on the telescope "
              "sphere)");
    }
}

// Rows of one emission chunk in a ground run: the sphere representation (docs/stream/format.md).
// A photon at P (relative to the sphere centre, radius R) with direction d gives the row P when
// |P| <= R, else the entry point P + s d with the smallest s > 0 where the line meets the sphere
// (path s, time + s / c). A line that misses the sphere is dropped and counted in misses_sphere.
// zem is the emission altitude: observation level + telescope z + z of the emission point P.
inline void finish_sphere_rows(const Emission &e, const Run &run,
                               obdeect::stream::PhotonChunk &chunk, Drops &drops) {
  using obdeect::stream::Column;
  const Telescope &telescope = run.telescopes.front();
  const double radius2 = telescope.radius_m * telescope.radius_m;
  const double emission_zem_base = run.observation_level_m + telescope.position.z;
  const std::size_t rows = chunk.rows();
  auto &x = chunk[Column::x];
  auto &y = chunk[Column::y];
  auto &z = chunk[Column::z];
  auto &cx = chunk[Column::cx];
  auto &cy = chunk[Column::cy];
  auto &cz = chunk[Column::cz];
  auto &t = chunk[Column::t];
  auto &w = chunk[Column::w];
  chunk[Column::zem].assign(rows, 0.0F);
  chunk[Column::path].assign(rows, 0.0F);
  auto &zem = chunk[Column::zem];
  auto &path = chunk[Column::path];
  std::size_t kept = 0;
  for (std::size_t i = 0; i < rows; ++i) {
    // Read the row before it is written to the position kept (kept <= i).
    const Vec3 start{x[i], y[i], z[i]};
    const Vec3 direction{cx[i], cy[i], cz[i]};
    double time = static_cast<double>(t[i]) + e.time_offset_ns;
    const double weight = static_cast<double>(w[i]) * e.bunch_size;
    const double emission_zem = emission_zem_base + start.z;
    Vec3 row = start;
    double s = 0.0;
    const double r2 = dot(start, start);
    if (r2 > radius2) {
      // Solve |start + s d|^2 = R^2 (d is a unit vector): s^2 + 2 b s + (r2 - R^2) = 0.
      const double b = dot(start, direction);
      const double disc = b * b - (r2 - radius2);
      s = -b - std::sqrt(std::max(0.0, disc));
      // Both roots have the same sign (r2 > R^2): the smaller one must lie ahead of the photon.
      if (!(disc >= 0.0) || !(s > 0.0)) {
        ++drops.misses_sphere;
        continue;
      }
      row = start + direction * s;
      time += s / speed_of_light_m_per_ns;
    }
    if (kept != i)
      for (auto &column : chunk.columns)
        if (!column.empty())
          column[kept] = column[i];
    x[kept] = static_cast<float>(row.x);
    y[kept] = static_cast<float>(row.y);
    z[kept] = static_cast<float>(row.z);
    t[kept] = static_cast<float>(time);
    w[kept] = static_cast<float>(weight);
    zem[kept] = static_cast<float>(emission_zem);
    path[kept] = static_cast<float>(s);
    ++kept;
  }
  chunk.truncate(kept);
}

// Rows of one emission chunk: bunch weight, time offset, move to z = 0 and zem (sky frame).
// With to_plane_z0 a row reaches the plane when z == 0 (no move) or z * cz < 0 (towards the
// plane from either side); other rows are dropped and counted in drops.not_on_plane. Kept rows
// are compacted to the front in their order, so the rows of the chunk stay in emission order.
// Ground runs write sphere rows instead (finish_sphere_rows).
inline void finish_rows(const Emission &e, const Run &run, obdeect::stream::PhotonChunk &chunk,
                        Drops &drops) {
  if (run.frame == obdeect::stream::Frame::ground) {
    finish_sphere_rows(e, run, chunk, drops);
    return;
  }
  using obdeect::stream::Column;
  const std::size_t rows = chunk.rows();
  auto &x = chunk[Column::x];
  auto &y = chunk[Column::y];
  auto &z = chunk[Column::z];
  auto &t = chunk[Column::t];
  auto &w = chunk[Column::w];
  std::size_t kept = 0;
  for (std::size_t i = 0; i < rows; ++i) {
    // Read the row before it is written to the position kept (kept <= i).
    float xi = x[i];
    float yi = y[i];
    float zi = z[i];
    double time = static_cast<double>(t[i]) + e.time_offset_ns;
    if (e.to_plane_z0 && zi != 0.0F) {
      const double zd = static_cast<double>(zi);
      const double cz = static_cast<double>(chunk[Column::cz][i]);
      if (!(zd * cz < 0.0)) {
        ++drops.not_on_plane;
        continue;
      }
      const double path = -zd / cz;
      xi = static_cast<float>(x[i] + path * chunk[Column::cx][i]);
      yi = static_cast<float>(y[i] + path * chunk[Column::cy][i]);
      zi = 0.0F;
      time += path / speed_of_light_m_per_ns;
    }
    const double weight = static_cast<double>(w[i]) * e.bunch_size;
    if (kept != i)
      for (auto &column : chunk.columns)
        if (!column.empty())
          column[kept] = column[i];
    x[kept] = xi;
    y[kept] = yi;
    z[kept] = zi;
    t[kept] = static_cast<float>(time);
    w[kept] = static_cast<float>(weight);
    ++kept;
  }
  chunk.truncate(kept);
  if (run.frame == obdeect::stream::Frame::sky) {
    auto &zem = chunk[Column::zem];
    zem.resize(kept);
    for (std::size_t i = 0; i < kept; ++i)
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
  const SphereTarget sphere = target_sphere(run);
  std::vector<std::vector<PlannedEmission>> planned;
  std::uint64_t k = 0;
  for (const auto &event : run.events) {
    auto &list = planned.emplace_back();
    for (const auto &e : event.emissions) {
      PlannedEmission p{e, source_seed(run.seed, k++)};
      if (std::holds_alternative<AutoTarget>(p.emission.source.target))
        p.emission.source.target = resolve_auto_target(p.emission.source, sphere);
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
  if (run.frame == obdeect::stream::Frame::ground)
    m["ground_representation"] = std::string("sphere");
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
      // Calls f(chunk) for every chunk of the event, in emission and photon order; dropped photons
      // are added to dropped_sink.
      Drops drops;
      Drops drops_again; // drops of the telescopes after the first (not counted)
      const auto generate = [&](auto &&f, Drops &drops_sink) {
        for (std::size_t k = 0; k < emissions.size(); ++k) {
          const auto &p = emissions[k];
          for (std::uint64_t done = 0; done < counts[k];) {
            const auto rows = static_cast<std::size_t>(
                std::min<std::uint64_t>(obdeect::stream::max_chunk_rows, counts[k] - done));
            emit(p.emission.source, p.seed, event, key, done, rows, chunk);
            detail::finish_rows(p.emission, run, chunk, drops_sink);
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
        generate([&](const obdeect::stream::PhotonChunk &c) { cached.push_back(c); }, drops);
      for (const auto &tel : run.telescopes) {
        // The photons of the event are the same for every telescope: count the drops once.
        Drops &drops_sink = &tel == &run.telescopes.front() ? drops : drops_again;
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
          generate(write_block, drops_sink);
        writer.write_metadata(
            RecordType::telescope_end,
            {{"telescope_id", static_cast<std::int64_t>(tel.id)}, {"photons", weight}});
        summary.rows += written;
        summary.weight += weight;
      }
      summary.not_on_plane += drops.not_on_plane;
      summary.misses_sphere += drops.misses_sphere;
      writer.write_metadata(RecordType::event_end, {});
    }
  }
  writer.write_metadata(RecordType::run_end,
                        {{"not_on_plane", static_cast<std::int64_t>(summary.not_on_plane)},
                         {"misses_sphere", static_cast<std::int64_t>(summary.misses_sphere)}});
  writer.finish();
  return summary;
}

} // namespace obdeect::sources
