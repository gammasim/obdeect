#pragma once

// Cherenkov-like track: a line segment emitting at a fixed angle to its direction, with
// importance sampling towards a sphere (docs/spec/track.md).

#include "dlamb/density.hpp"
#include "dlamb/source.hpp"
#include "dlamb/vec3.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <variant>
#include <vector>

namespace dlamb {

struct Track {
  Vec3 start{};                   // upper end of the segment
  Vec3 direction{0.0, 0.0, -1.0}; // unit direction of motion
  double length_m{};
  double cos_theta{1.0}; // cosine of the emission angle, 1 / (n beta)
  double beta{1.0};
  double time_ns{}; // emission time at the start
  Spectrum spectrum{Fixed{400.0}};
  double photons{}; // mean number emitted around the whole azimuth
  SphereTarget target{};
};

// Azimuths (relative to the basis of the track direction) that reach the sphere:
// centre +- half; fraction half / pi.
struct TrackArc {
  bool hit{};
  double centre{};
  double half{};
};

inline void validate(const Track &t) {
  using detail::require;
  require(std::isfinite(t.photons) && t.photons >= 0.0 && t.photons <= max_photon_count,
          "track: photons must be finite, >= 0 and at most 2^53");
  require(within_limits(t.start), "track: start out of range");
  require(is_unit(t.direction), "track: direction must be a unit vector");
  require(within_limits(t.length_m) && t.length_m >= 0.0, "track: length out of range");
  require(t.cos_theta > 0.0 && t.cos_theta <= 1.0, "track: cos_theta must lie in (0, 1]");
  require(t.beta >= 0.0 && t.beta <= 1.0, "track: beta must lie in [0, 1]");
  require(within_limits(t.time_ns), "track: time out of range");
  require(within_limits(t.target.centre) && positive_within_limits(t.target.radius_m),
          "track: sphere centre or radius out of range");
  require(std::holds_alternative<Fixed>(t.spectrum) ||
              std::holds_alternative<InverseSquare>(t.spectrum),
          "track: the spectrum must be fixed or inverse-square");
  validate_spectrum(t.spectrum);
}

namespace detail {

// Arc without validation (for loops over a validated track).
inline TrackArc arc_at(const Track &t, double s) {
  const Vec3 point = t.start + t.direction * s;
  const Vec3 to_centre = t.target.centre - point;
  const double distance = norm(to_centre);
  const double radius = t.target.radius_m;
  if (distance <= radius)
    return {true, 0.0, pi};
  const double ratio = radius / distance;
  const double cos_alpha = std::sqrt(1.0 - ratio * ratio);
  const Vec3 unit = to_centre * (1.0 / distance);
  const double cos_b = dot(t.direction, unit);
  const Vec3 perp = unit - t.direction * cos_b;
  const double sin_b = norm(perp);
  const Basis basis = orthonormal_basis(t.direction);
  const double phi0 = std::atan2(dot(perp, basis.e2), dot(perp, basis.e1));
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - t.cos_theta * t.cos_theta));
  const double a = cos_alpha - t.cos_theta * cos_b;
  const double b = sin_theta * sin_b;
  if (b <= 0.0)
    return a < 0.0 ? TrackArc{true, 0.0, pi} : TrackArc{};
  const double x = a / b;
  if (x >= 1.0)
    return {};
  if (x <= -1.0)
    return {true, 0.0, pi};
  return {true, phi0, std::acos(x)};
}

} // namespace detail

// Arc of azimuths from the point at distance s along a (validated) track.
inline TrackArc track_arc(const Track &t, double s) {
  validate(t);
  return detail::arc_at(t, s);
}

inline std::uint64_t track_count(const Track &t, std::uint64_t seed, std::uint32_t event,
                                 std::uint32_t telescope) {
  validate(t);
  return round_count(t.photons, seed, event, telescope);
}

// Appends the photons of a track that reach the sphere (misses consume their index).
inline void emit_track(const Track &t, std::uint64_t seed, std::uint32_t event,
                       std::uint32_t telescope, std::vector<Photon> &out) {
  const std::uint64_t n = track_count(t, seed, event, telescope);
  const Basis basis = orthonormal_basis(t.direction);
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - t.cos_theta * t.cos_theta));
  for (std::uint64_t k = 0; k < n; ++k) {
    const double s = t.length_m * philox_draw(seed, event, telescope, k, Stream::target)[0];
    const TrackArc arc = detail::arc_at(t, s);
    if (!arc.hit)
      continue;
    const double u = philox_draw(seed, event, telescope, k, Stream::direction)[0];
    const double phi = arc.centre + arc.half * (2.0 * u - 1.0);
    Photon p;
    p.position = t.start + t.direction * s;
    p.direction = direction_about(t.direction, basis, t.cos_theta, sin_theta, phi);
    p.time_ns = t.beta > 0.0 ? t.time_ns + s / (t.beta * speed_of_light_m_per_ns) : t.time_ns;
    p.weight = arc.half / pi;
    p.wavelength_nm =
        sample(t.spectrum, philox_draw(seed, event, telescope, k, Stream::wavelength));
    out.push_back(p);
  }
}

// Expected sum of weights (midpoint rule along the track).
inline double expected_track_weight(const Track &t, int steps = 20000) {
  validate(t);
  if (!(t.length_m > 0.0) || !(t.photons > 0.0))
    return 0.0;
  double sum = 0.0;
  for (int j = 0; j < steps; ++j)
    sum += detail::arc_at(t, t.length_m * (j + 0.5) / steps).half / pi;
  return t.photons * sum / steps;
}

} // namespace dlamb
