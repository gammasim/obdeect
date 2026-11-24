#pragma once

// Per-photon sampling of every source geometry and emission of photon chunks
// (docs/sources/spec/source-core.md).

#include "obdeect/sources/source.hpp"

#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <variant>

namespace obdeect::sources {

// Rejection attempts of one photon towards a disc target before the sampler gives up. A photon
// needs about 1 / (fraction of light reaching the disc) attempts, so this is a bug guard only.
inline constexpr std::uint32_t max_disc_attempts = 1U << 24;

// Per-emission constants of the sampler, computed once per emission (not per photon).
struct Prepared {
  double reference{1.0}; // line_reference_fraction of a line with a sphere target, else 1
  Basis basis{};         // orthonormal basis of the emission axis (point, line or sky background)
};

namespace detail {

// Direction at angle theta (cosine cos_theta) and azimuth about `axis`, with its basis cached.
inline Vec3 cone_direction_about(const Vec3 &axis, const Basis &basis, double cos_theta,
                                 double azimuth) {
  const double sin_theta = std::sqrt(std::max(0.0, (1.0 - cos_theta) * (1.0 + cos_theta)));
  return direction_about(axis, basis, cos_theta, sin_theta, azimuth);
}

inline Vec3 cone_direction(const Vec3 &axis, double cos_theta, double azimuth) {
  return cone_direction_about(axis, orthonormal_basis(axis), cos_theta, azimuth);
}

inline double line_position(const Line &l, double u) {
  return l.profile.empty() ? l.length_m * u : l.profile.sample(u);
}

// True when the ray from origin along the unit direction meets the disc in front of the origin.
inline bool hits_disc(const Vec3 &origin, const Vec3 &direction, const DiscTarget &disc) {
  const double denominator = dot(direction, disc.normal);
  if (denominator == 0.0)
    return false;
  const double t = dot(disc.centre - origin, disc.normal) / denominator;
  if (!(t > 0.0))
    return false;
  const Vec3 offset = origin + direction * t - disc.centre;
  return dot(offset, offset) <= disc.radius_m * disc.radius_m;
}

// A point with a disc target: the photons are the directions of the profile that reach the disc,
// each with weight 1, so the candidate count is the expected number of hits (expected_count).
// The directions are drawn in the smaller of two regions (disc_region): in the profile's cone
// (from the profile; kept when they hit the disc), or in the cone enclosing the disc (uniform in
// solid angle; kept when they hit the disc and u[2] < I / I_max, a rejection on the profile). Both
// give the profile's distribution over the disc. Attempt 0 uses the direction draws ud; later
// attempts their own draws.
inline void sample_point_disc(const Point &e, const Basis &basis, const DiscTarget &disc,
                              const Uniforms &ud, std::uint64_t seed, std::uint32_t event,
                              std::uint32_t telescope, std::uint64_t index, Photon &p) {
  const DiscRegion region = disc_region(e, disc);
  const Basis cone_basis = region.disc_cone ? orthonormal_basis(region.axis) : basis;
  const double maximum = intensity_maximum(e.angular);
  for (std::uint32_t attempt = 0; attempt < max_disc_attempts; ++attempt) {
    const Uniforms u = attempt == 0 ? ud
                                    : philox_draw_attempt(seed, event, telescope, index,
                                                          Stream::direction, attempt);
    Vec3 direction;
    if (region.disc_cone) {
      direction =
          cone_direction_about(region.axis, cone_basis, 1.0 - u[0] * region.u_cone, two_pi * u[1]);
      const double cosine = std::clamp(dot(direction, e.axis), -1.0, 1.0);
      if (!(u[2] * maximum < angle_intensity(e.angular, 1.0 - cosine)))
        continue;
    } else {
      direction =
          cone_direction_about(e.axis, basis, 1.0 - angle_quantile(e.angular, u[0]), two_pi * u[1]);
    }
    if (hits_disc(e.position, direction, disc)) {
      p.direction = direction;
      return;
    }
  }
  throw std::runtime_error("disc target: no photon reached the disc within " +
                           std::to_string(max_disc_attempts) +
                           " attempts (the disc subtends a negligible fraction of the light)");
}

inline void sample_point(const Point &e, const Basis &basis, const Target &target,
                         const Uniforms &ud, std::uint64_t seed, std::uint32_t event,
                         std::uint32_t telescope, std::uint64_t index, Photon &p) {
  p.position = e.position;
  if (const auto *disc = std::get_if<DiscTarget>(&target)) {
    sample_point_disc(e, basis, *disc, ud, seed, event, telescope, index, p);
  } else if (const auto *sphere = std::get_if<SphereTarget>(&target)) {
    const Vec3 to_centre = sphere->centre - e.position;
    const double d = norm(to_centre);
    const double ratio = sphere->radius_m / d;
    const double cos_max = std::sqrt(1.0 - ratio * ratio);
    p.direction =
        cone_direction(to_centre * (1.0 / d), 1.0 - ud[0] * (1.0 - cos_max), two_pi * ud[1]);
  } else {
    p.direction =
        cone_direction_about(e.axis, basis, 1.0 - angle_quantile(e.angular, ud[0]), two_pi * ud[1]);
  }
}

// False when the photon misses the target (the row is dropped).
inline bool sample_line(const Line &l, const Basis &basis, const Target &target, double reference,
                        const Uniforms &ud, std::uint64_t seed, std::uint32_t event,
                        std::uint32_t telescope, std::uint64_t index, Photon &p) {
  const double s = line_position(l, philox_draw(seed, event, telescope, index, Stream::target)[0]);
  p.position = l.start + l.direction * s;
  if (l.beta > 0.0)
    p.time_ns += s / (l.beta * speed_of_light_m_per_ns);
  const auto *sphere = std::get_if<SphereTarget>(&target);
  if (sphere && !(reference > 0.0))
    return false; // no photon of this line reaches the sphere (the expected count is 0)
  if (sphere && std::holds_alternative<Isotropic>(l.angular)) {
    // Cone about the sphere seen from the point; weight = fraction of 4 pi / reference.
    const Vec3 to_centre = sphere->centre - p.position;
    const double d = norm(to_centre);
    const double fraction = sphere_fraction(p.position, *sphere);
    p.direction =
        d > 0.0
            ? cone_direction(to_centre * (1.0 / d), 1.0 - ud[0] * 2.0 * fraction, two_pi * ud[1])
            : cone_direction_about(l.direction, basis, 1.0 - 2.0 * ud[0], two_pi * ud[1]);
    p.weight = fraction / reference;
  } else if (sphere) {
    const double cos_theta = 1.0 - std::get<FixedAngle>(l.angular).u0;
    const Arc arc = arc_towards(p.position, l.direction, cos_theta, *sphere);
    if (!arc.hit)
      return false;
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    p.direction = direction_about(l.direction, basis, cos_theta, sin_theta,
                                  arc.centre + arc.half * (2.0 * ud[0] - 1.0));
    p.weight = arc.half / pi / reference;
  } else {
    p.direction = cone_direction_about(l.direction, basis, 1.0 - angle_quantile(l.angular, ud[0]),
                                       two_pi * ud[1]);
  }
  return true;
}

inline void sample_plane_wave(const PlaneWave &w, std::uint64_t seed, std::uint32_t event,
                              std::uint32_t telescope, std::uint64_t index, Photon &p) {
  const Uniforms ut = philox_draw(seed, event, telescope, index, Stream::target);
  const double r = w.disc_radius_m * std::sqrt(ut[0]);
  const double phi = two_pi * ut[1];
  const Vec3 target{r * std::cos(phi), r * std::sin(phi), w.disc_z_m};
  const Vec3 centre{0.0, 0.0, w.disc_z_m};
  p.direction = -w.towards_source;
  if (w.start_distance_m) {
    // Start plane: (X - centre) . towards_source = start_distance.
    const double start_distance = *w.start_distance_m;
    if (std::isinf(w.distance_m)) {
      p.position =
          target + w.towards_source * (start_distance - dot(target - centre, w.towards_source));
      p.time_ns += dot(p.position - centre, p.direction) / speed_of_light_m_per_ns;
      return;
    }
    const Vec3 source = centre + w.towards_source * w.distance_m;
    p.direction = normalised(target - source);
    // The source lies at depth distance_m along towards_source: the photons start there when the
    // start plane is beyond it, else where the ray crosses the plane.
    p.position = w.distance_m <= start_distance
                     ? source
                     : source + p.direction * ((start_distance - w.distance_m) /
                                               dot(p.direction, w.towards_source));
    p.time_ns += (norm(p.position - source) - w.distance_m) / speed_of_light_m_per_ns;
    return;
  }
  const double start_z = plane_wave_start_z(w);
  if (std::isinf(w.distance_m)) {
    p.position = target + p.direction * ((start_z - target.z) / p.direction.z);
    p.time_ns += dot(p.position - centre, p.direction) / speed_of_light_m_per_ns;
    return;
  }
  const Vec3 source = centre + w.towards_source * w.distance_m;
  p.direction = normalised(target - source);
  // The photons start at the source when it is already on the disc side of the start plane,
  // otherwise where the ray crosses the plane.
  const bool source_is_on_path = p.direction.z * (source.z - start_z) >= 0.0;
  p.position =
      source_is_on_path ? source : target + p.direction * ((start_z - target.z) / p.direction.z);
  p.time_ns += (norm(p.position - source) - w.distance_m) / speed_of_light_m_per_ns;
}

inline void sample_sky_background(const SkyBackground &b, const Basis &basis, const Uniforms &ud,
                                  std::uint64_t seed, std::uint32_t event, std::uint32_t telescope,
                                  std::uint64_t index, Photon &p) {
  const double cos_cone = std::cos(b.viewcone_deg * pi / 180.0);
  p.direction = cone_direction_about(b.axis, basis, 1.0 - ud[0] * (1.0 - cos_cone), two_pi * ud[1]);
  // Start point: uniform over the disc seen along the direction, lifted onto the sphere.
  const Uniforms ut = philox_draw(seed, event, telescope, index, Stream::target);
  const double r = std::sqrt(ut[0]);
  const double phi = two_pi * ut[1];
  const double h = std::sqrt(std::max(0.0, 1.0 - r * r));
  const Basis start_basis = orthonormal_basis(p.direction);
  const Vec3 unit =
      start_basis.e1 * (r * std::cos(phi)) + start_basis.e2 * (r * std::sin(phi)) - p.direction * h;
  p.position = b.centre + unit * b.radius_m;
  p.time_ns -= h * b.radius_m / speed_of_light_m_per_ns;
}

} // namespace detail

// Reference fraction of a source for sample_photon (1 unless a line has a sphere target).
inline double line_reference(const Source &s) {
  const auto *l = std::get_if<Line>(&s.geometry);
  const auto *sphere = std::get_if<SphereTarget>(&s.target);
  return l && sphere ? line_reference_fraction(*l, *sphere) : 1.0;
}

// The constants of a source for the sampler (its validation is the caller's business).
inline Prepared prepare(const Source &s) {
  Prepared prepared;
  prepared.reference = line_reference(s);
  if (const auto *e = std::get_if<Point>(&s.geometry))
    prepared.basis = orthonormal_basis(e->axis);
  else if (const auto *l = std::get_if<Line>(&s.geometry))
    prepared.basis = orthonormal_basis(l->direction);
  else if (const auto *b = std::get_if<SkyBackground>(&s.geometry))
    prepared.basis = orthonormal_basis(b->axis);
  return prepared;
}

// Candidate photon `index` (no validation: use emit() or validate() first). False when a line
// photon misses its target.
inline bool sample_photon(const Source &s, const Prepared &prepared, std::uint64_t seed,
                          std::uint32_t event, std::uint32_t telescope, std::uint64_t index,
                          Photon &p) {
  // Draws that a profile ignores are not computed (same result, fewer Philox calls).
  const Uniforms none{};
  p = Photon{};
  p.time_ns = sample(s.pulse, ignores_draw(s.pulse)
                                  ? none
                                  : philox_draw(seed, event, telescope, index, Stream::pulse));
  p.wavelength_nm =
      sample(s.spectrum, ignores_draw(s.spectrum)
                             ? none
                             : philox_draw(seed, event, telescope, index, Stream::wavelength));
  if (const auto *w = std::get_if<PlaneWave>(&s.geometry)) {
    detail::sample_plane_wave(*w, seed, event, telescope, index, p);
    return true;
  }
  const Uniforms ud = philox_draw(seed, event, telescope, index, Stream::direction);
  if (const auto *e = std::get_if<Point>(&s.geometry)) {
    detail::sample_point(*e, prepared.basis, s.target, ud, seed, event, telescope, index, p);
    return true;
  }
  if (const auto *l = std::get_if<Line>(&s.geometry))
    return detail::sample_line(*l, prepared.basis, s.target, prepared.reference, ud, seed, event,
                               telescope, index, p);
  detail::sample_sky_background(std::get<SkyBackground>(s.geometry), prepared.basis, ud, seed,
                                event, telescope, index, p);
  return true;
}

// Candidates first .. first + candidates - 1 as float32 rows (source_columns) of a validated
// source with its constants prepare(s). Line photons that miss their target are dropped; returns
// the rows written. chunk.first_index is left to the caller.
inline std::size_t emit_prepared(const Source &s, const Prepared &prepared, std::uint64_t seed,
                                 std::uint32_t event, std::uint32_t telescope, std::uint64_t first,
                                 std::size_t candidates, obdeect::stream::PhotonChunk &chunk) {
  if (candidates > obdeect::stream::max_chunk_rows)
    throw std::invalid_argument("a chunk holds at most 65536 rows");
  chunk.resize(candidates, source_columns);
  chunk.telescope_id = telescope;
  using obdeect::stream::Column;
  std::size_t row = 0;
  Photon p;
  for (std::size_t i = 0; i < candidates; ++i) {
    if (!sample_photon(s, prepared, seed, event, telescope, first + i, p))
      continue;
    chunk[Column::x][row] = static_cast<float>(p.position.x);
    chunk[Column::y][row] = static_cast<float>(p.position.y);
    chunk[Column::z][row] = static_cast<float>(p.position.z);
    chunk[Column::cx][row] = static_cast<float>(p.direction.x);
    chunk[Column::cy][row] = static_cast<float>(p.direction.y);
    chunk[Column::cz][row] = static_cast<float>(p.direction.z);
    chunk[Column::t][row] = static_cast<float>(p.time_ns);
    chunk[Column::w][row] = static_cast<float>(p.weight);
    chunk[Column::wavelength][row] = static_cast<float>(p.wavelength_nm);
    ++row;
  }
  chunk.truncate(row);
  return row;
}

// Candidates first .. first + candidates - 1 of a source (validates it and prepares it per call;
// the run writer uses emit_prepared with the constants of each emission).
inline std::size_t emit(const Source &s, std::uint64_t seed, std::uint32_t event,
                        std::uint32_t telescope, std::uint64_t first, std::size_t candidates,
                        obdeect::stream::PhotonChunk &chunk) {
  validate(s);
  return emit_prepared(s, prepare(s), seed, event, telescope, first, candidates, chunk);
}

// Expected sum of weights of the emitted photons: the candidate count times the mean weight
// (midpoint rule along a line with a sphere target; exact otherwise: a disc target gives weight 1
// per candidate).
inline double expected_weight(const Source &s, int steps = 2000) {
  const double count = expected_count(s);
  if (!(count > 0.0))
    return 0.0;
  const auto *l = std::get_if<Line>(&s.geometry);
  const auto *sphere = std::get_if<SphereTarget>(&s.target);
  if (!l || !sphere)
    return count;
  const bool isotropic = std::holds_alternative<Isotropic>(l->angular);
  const double cos_theta = isotropic ? 0.0 : 1.0 - std::get<FixedAngle>(l->angular).u0;
  const double reference = line_reference_fraction(*l, *sphere);
  if (!(reference > 0.0))
    return 0.0;
  double sum = 0.0;
  const int line_steps = steps * 10;
  for (int j = 0; j < line_steps; ++j) {
    const Vec3 point = l->start + l->direction * detail::line_position(*l, (j + 0.5) / line_steps);
    sum += (isotropic ? sphere_fraction(point, *sphere)
                      : detail::arc_towards(point, l->direction, cos_theta, *sphere).half / pi) /
           reference;
  }
  return count * sum / line_steps;
}

} // namespace obdeect::sources
