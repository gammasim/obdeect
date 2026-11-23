#pragma once

// Per-photon sampling of every source geometry and emission of photon chunks
// (docs/spec/source-core.md).

#include "obdeect/sources/source.hpp"

#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <variant>

namespace obdeect::sources {

namespace detail {

inline Vec3 cone_direction(const Vec3 &axis, double cos_theta, double azimuth) {
  const Basis basis = orthonormal_basis(axis);
  const double sin_theta = std::sqrt(std::max(0.0, (1.0 - cos_theta) * (1.0 + cos_theta)));
  return direction_about(axis, basis, cos_theta, sin_theta, azimuth);
}

inline double line_position(const Line &l, double u) {
  return l.profile.empty() ? l.length_m * u : l.profile.sample(u);
}

inline void sample_point(const Point &e, const Target &target, const Uniforms &ud,
                         std::uint64_t seed, std::uint32_t event, std::uint32_t telescope,
                         std::uint64_t index, Photon &p) {
  p.position = e.position;
  if (const auto *disc = std::get_if<DiscTarget>(&target)) {
    const Uniforms ut = philox_draw(seed, event, telescope, index, Stream::target);
    const Basis basis = orthonormal_basis(disc->normal);
    const double r = disc->radius_m * std::sqrt(ut[0]);
    const double phi = two_pi * ut[1];
    const Vec3 point =
        disc->centre + basis.e1 * (r * std::cos(phi)) + basis.e2 * (r * std::sin(phi));
    const Vec3 to_point = point - e.position;
    const double distance = norm(to_point);
    p.direction = to_point * (1.0 / distance);
    const double h = std::abs(dot(disc->centre - e.position, disc->normal));
    const double cosine = std::clamp(dot(p.direction, e.axis), -1.0, 1.0);
    const double ratio = h / distance;
    p.weight = std::visit([&](const auto &a) { return a.intensity(1.0 - cosine); }, e.angular) /
               intensity_maximum(e.angular) * ratio * ratio * ratio;
  } else if (const auto *sphere = std::get_if<SphereTarget>(&target)) {
    const Vec3 to_centre = sphere->centre - e.position;
    const double d = norm(to_centre);
    const double ratio = sphere->radius_m / d;
    const double cos_max = std::sqrt(1.0 - ratio * ratio);
    p.direction =
        cone_direction(to_centre * (1.0 / d), 1.0 - ud[0] * (1.0 - cos_max), two_pi * ud[1]);
  } else {
    p.direction = cone_direction(e.axis, 1.0 - angle_quantile(e.angular, ud[0]), two_pi * ud[1]);
  }
}

// False when the photon misses the target (the row is dropped).
inline bool sample_line(const Line &l, const Target &target, double reference, const Uniforms &ud,
                        std::uint64_t seed, std::uint32_t event, std::uint32_t telescope,
                        std::uint64_t index, Photon &p) {
  const double s = line_position(l, philox_draw(seed, event, telescope, index, Stream::target)[0]);
  p.position = l.start + l.direction * s;
  if (l.beta > 0.0)
    p.time_ns += s / (l.beta * speed_of_light_m_per_ns);
  const auto *sphere = std::get_if<SphereTarget>(&target);
  if (sphere && std::holds_alternative<Isotropic>(l.angular)) {
    // Cone about the sphere seen from the point; weight = fraction of 4 pi / reference.
    const Vec3 to_centre = sphere->centre - p.position;
    const double d = norm(to_centre);
    const double fraction = sphere_fraction(p.position, *sphere);
    p.direction = d > 0.0 ? cone_direction(to_centre * (1.0 / d), 1.0 - ud[0] * 2.0 * fraction,
                                           two_pi * ud[1])
                          : cone_direction(l.direction, 1.0 - 2.0 * ud[0], two_pi * ud[1]);
    p.weight = fraction / reference;
  } else if (sphere) {
    const double cos_theta = 1.0 - std::get<FixedAngle>(l.angular).u0;
    const Arc arc = arc_towards(p.position, l.direction, cos_theta, *sphere);
    if (!arc.hit)
      return false;
    const Basis basis = orthonormal_basis(l.direction);
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    p.direction = direction_about(l.direction, basis, cos_theta, sin_theta,
                                  arc.centre + arc.half * (2.0 * ud[0] - 1.0));
    p.weight = arc.half / pi / reference;
  } else {
    p.direction =
        cone_direction(l.direction, 1.0 - angle_quantile(l.angular, ud[0]), two_pi * ud[1]);
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
  if (std::isinf(w.distance_m)) {
    p.position = target + p.direction * ((w.start_z_m - target.z) / p.direction.z);
    p.time_ns += dot(p.position - centre, p.direction) / speed_of_light_m_per_ns;
    return;
  }
  const Vec3 source = centre + w.towards_source * w.distance_m;
  p.direction = normalised(target - source);
  // The photons start at the source when it is already on the disc side of the start plane,
  // otherwise where the ray crosses the plane.
  const bool source_is_on_path = p.direction.z * (source.z - w.start_z_m) >= 0.0;
  p.position = source_is_on_path
                   ? source
                   : target + p.direction * ((w.start_z_m - target.z) / p.direction.z);
  p.time_ns += (norm(p.position - source) - w.distance_m) / speed_of_light_m_per_ns;
}

inline void sample_sky_background(const SkyBackground &b, const Uniforms &ud, std::uint64_t seed,
                                  std::uint32_t event, std::uint32_t telescope, std::uint64_t index,
                                  Photon &p) {
  const double cos_cone = std::cos(b.viewcone_deg * pi / 180.0);
  p.direction = cone_direction(b.axis, 1.0 - ud[0] * (1.0 - cos_cone), two_pi * ud[1]);
  // Start point: uniform over the disc seen along the direction, lifted onto the sphere.
  const Uniforms ut = philox_draw(seed, event, telescope, index, Stream::target);
  const double r = std::sqrt(ut[0]);
  const double phi = two_pi * ut[1];
  const double h = std::sqrt(std::max(0.0, 1.0 - r * r));
  const Basis basis = orthonormal_basis(p.direction);
  const Vec3 unit =
      basis.e1 * (r * std::cos(phi)) + basis.e2 * (r * std::sin(phi)) - p.direction * h;
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

// Candidate photon `index` (no validation: use emit() or validate() first); `reference` is
// line_reference(s). False when a line photon misses its target.
inline bool sample_photon(const Source &s, std::uint64_t seed, std::uint32_t event,
                          std::uint32_t telescope, std::uint64_t index, Photon &p,
                          double reference) {
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
    detail::sample_point(*e, s.target, ud, seed, event, telescope, index, p);
    return true;
  }
  if (const auto *l = std::get_if<Line>(&s.geometry))
    return detail::sample_line(*l, s.target, reference, ud, seed, event, telescope, index, p);
  detail::sample_sky_background(std::get<SkyBackground>(s.geometry), ud, seed, event, telescope,
                                index, p);
  return true;
}

// Candidates first .. first + candidates - 1 as float32 rows (source_columns). Line photons that
// miss their target are dropped; returns the rows written. chunk.first_index is left to the
// caller.
inline std::size_t emit(const Source &s, std::uint64_t seed, std::uint32_t event,
                        std::uint32_t telescope, std::uint64_t first, std::size_t candidates,
                        obdeect::stream::PhotonChunk &chunk) {
  validate(s);
  if (candidates > obdeect::stream::max_chunk_rows)
    throw std::invalid_argument("a chunk holds at most 65536 rows");
  chunk.resize(candidates, source_columns);
  chunk.telescope_id = telescope;
  using obdeect::stream::Column;
  std::size_t row = 0;
  const double reference = line_reference(s);
  Photon p;
  for (std::size_t i = 0; i < candidates; ++i) {
    if (!sample_photon(s, seed, event, telescope, first + i, p, reference))
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

// Expected sum of weights of the emitted photons: the candidate count times the mean weight
// (midpoint rule for a disc target and for a line with a sphere target; exact otherwise).
inline double expected_weight(const Source &s, int steps = 2000) {
  const double count = expected_count(s);
  if (!(count > 0.0))
    return 0.0;
  double sum = 0.0;
  if (const auto *p = std::get_if<Point>(&s.geometry)) {
    const auto *disc = std::get_if<DiscTarget>(&s.target);
    if (!disc)
      return count;
    // Mean of the disc weight over uniform points of the disc (midpoint in r^2 and azimuth).
    const Basis basis = orthonormal_basis(disc->normal);
    const double h = std::abs(dot(disc->centre - p->position, disc->normal));
    const double imax = intensity_maximum(p->angular);
    for (int i = 0; i < steps; ++i) {
      const double r = disc->radius_m * std::sqrt((i + 0.5) / steps);
      for (int k = 0; k < steps; ++k) {
        const double phi = two_pi * (k + 0.5) / steps;
        const Vec3 to_point = disc->centre + basis.e1 * (r * std::cos(phi)) +
                              basis.e2 * (r * std::sin(phi)) - p->position;
        const double distance = norm(to_point);
        const double cosine = std::clamp(dot(to_point, p->axis) / distance, -1.0, 1.0);
        const double ratio = h / distance;
        sum += std::visit([&](const auto &a) { return a.intensity(1.0 - cosine); }, p->angular) /
               imax * ratio * ratio * ratio;
      }
    }
    return count * sum / (static_cast<double>(steps) * steps);
  }
  const auto *l = std::get_if<Line>(&s.geometry);
  const auto *sphere = std::get_if<SphereTarget>(&s.target);
  if (!l || !sphere)
    return count;
  const bool isotropic = std::holds_alternative<Isotropic>(l->angular);
  const double cos_theta = isotropic ? 0.0 : 1.0 - std::get<FixedAngle>(l->angular).u0;
  const double reference = line_reference_fraction(*l, *sphere);
  if (!(reference > 0.0))
    return 0.0;
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
