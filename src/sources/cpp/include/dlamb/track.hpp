#pragma once

#include "dlamb/density.hpp"
#include "dlamb/source.hpp"
#include "dlamb/streams.hpp"
#include "dlamb/vec3.hpp"

#include "gstrickd/random.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <variant>
#include <vector>

namespace dlamb {

// Track: a moving line source (LightEmission LightSource with length and beta) that emits
// at a fixed angle to its direction of motion, as the Cherenkov light of a charged particle
// (fake-muon). Photons are uniform along the segment and have a uniform azimuth. The
// sphere target keeps only the photons that reach it, with weight = the fraction of the
// azimuths that do (importance sampling; LightEmission keeps the same photons as bunches).
struct Track {
  Vec3 start{};                   // Upper end of the segment, altitude above observation level [m].
  Vec3 direction{0.0, 0.0, -1.0}; // Unit direction of motion, the emission axis.
  double length_m{};              // Segment length along the direction [m].
  double cos_theta{1.0};          // Cosine of the emission angle to the direction, 1/(n beta).
  double beta{1.0};               // Speed in units of c; 0 = instantaneous along the segment.
  double time_ns{};               // Emission time at the start of the segment [ns].
  Spectrum spectrum{Fixed{400.0}};
  double photons{};      // Mean photons emitted into 4 pi from the segment.
  SphereTarget target{}; // Telescope sphere.
};

// Throws std::invalid_argument for an inconsistent track.
inline void validate(const Track &track) {
  const auto fail = [](const char *message) { throw std::invalid_argument(message); };
  if (!(track.photons >= 0.0) || track.photons > max_photon_count)
    fail("track photons must be finite, >= 0 and at most 2^53");
  if (!magnitude_ok(track.start) || !is_finite(track.direction) ||
      std::abs(norm(track.direction) - 1.0) > 1e-9)
    fail("a track needs a finite start within 1e30 and a unit direction");
  if (!magnitude_ok(track.length_m) || track.length_m < 0.0)
    fail("track length must be finite, >= 0 and at most 1e30");
  if (!(track.cos_theta > 0.0 && track.cos_theta <= 1.0))
    fail("track cos(theta) must lie in (0, 1]");
  if (!(track.beta >= 0.0 && track.beta <= 1.0))
    fail("track beta must lie in [0, 1]");
  if (!magnitude_ok(track.time_ns))
    fail("track time must be finite and at most 1e30 in magnitude");
  if (!magnitude_ok(track.target.centre) || !positive_magnitude_ok(track.target.radius_m))
    fail("a track needs a sphere target with a centre and a positive radius, within 1e30");
  if (const auto *fixed = std::get_if<Fixed>(&track.spectrum)) {
    if (!positive_magnitude_ok(fixed->value))
      fail("a fixed wavelength must be positive and within 1e30");
  } else if (const auto *inv = std::get_if<InverseSquare>(&track.spectrum)) {
    if (!positive_magnitude_ok(inv->low) || !positive_magnitude_ok(inv->high) ||
        !(inv->high > inv->low))
      fail("an inverse-square wavelength range needs 0 < low < high, within 1e30");
  } else {
    fail("a track takes a fixed or inverse-square spectrum");
  }
}

// Photons of one track in one (event, telescope): the mean rounded up with probability
// equal to its fractional part (as Source, stream Stream::count).
[[nodiscard]] inline std::uint64_t track_count(const Track &track, std::uint64_t seed,
                                               std::uint32_t event, std::uint32_t telescope) {
  const double mean = track.photons;
  const double whole = std::floor(mean);
  const auto u = gstrickd::uniform4(
      {seed, stage_sources, event, telescope, 0, static_cast<std::uint32_t>(Stream::count)});
  return static_cast<std::uint64_t>(whole) + (u[0] < mean - whole ? 1U : 0U);
}

// Azimuths (about the track direction) of the photons that reach the sphere from a point
// at distance s along the track: those with azimuth in [centre - half, centre + half].
struct TrackArc {
  bool hit{};
  double centre{};
  double half{};
};

// Photon from a point P in direction e (angle theta to the track axis) reaches the sphere
// if its angle to the direction of the centre is below alpha = asin(R / |C - P|) (the ray
// points away from the sphere in its forward direction). With cos(angle) = cos(theta)
// cos(beta) + sin(theta) sin(beta) cos(phi - phi0) this gives an arc in phi. A point inside
// the sphere is hit in every direction.
[[nodiscard]] inline TrackArc track_arc(const Track &track, double s) {
  const Vec3 point = track.start + track.direction * s;
  const Vec3 to_centre = track.target.centre - point;
  const double distance = norm(to_centre);
  const double radius = track.target.radius_m;
  if (distance <= radius)
    return {true, 0.0, pi};
  const double ratio = radius / distance;
  const double cos_alpha = std::sqrt(1.0 - ratio * ratio);
  const Vec3 unit = to_centre * (1.0 / distance);
  const double cos_beta = dot(track.direction, unit);
  const Vec3 perpendicular = unit - track.direction * cos_beta;
  const Basis basis = orthonormal_basis(track.direction);
  const double sin_beta = norm(perpendicular);
  const double phi0 = std::atan2(dot(perpendicular, basis.e2), dot(perpendicular, basis.e1));
  const double cos_t = track.cos_theta;
  const double a = cos_alpha - cos_t * cos_beta;
  const double b = std::sqrt(std::max(0.0, 1.0 - cos_t * cos_t)) * sin_beta;
  if (b <= 0.0)
    return a < 0.0 ? TrackArc{true, 0.0, pi} : TrackArc{};
  const double x = a / b;
  if (x >= 1.0)
    return {};
  return {true, phi0, x <= -1.0 ? pi : std::acos(x)};
}

// Appends the photons of one track in one (event, telescope) that reach its sphere. Photon
// k draws its position, azimuth and wavelength from its own keys (Philox), so the result
// does not depend on the order of the photons. Times are emission times at the source.
inline void emit_track(const Track &track, std::uint64_t seed, std::uint32_t event,
                       std::uint32_t telescope, std::vector<Photon> &out) {
  validate(track);
  const std::uint64_t n = track_count(track, seed, event, telescope);
  const Basis basis = orthonormal_basis(track.direction);
  const double cos_t = track.cos_theta;
  const double sin_t = std::sqrt(std::max(0.0, 1.0 - cos_t * cos_t));
  for (std::uint64_t k = 0; k < n; ++k) {
    const auto key = [&](Stream stream) {
      return gstrickd::DrawKey{
          seed, stage_sources, event, telescope, k, static_cast<std::uint32_t>(stream)};
    };
    const double s = track.length_m * gstrickd::uniform4(key(Stream::target))[0];
    const TrackArc arc = track_arc(track, s);
    if (!arc.hit)
      continue;
    const double u = gstrickd::uniform4(key(Stream::direction))[0];
    const double phi = arc.centre + arc.half * (2.0 * u - 1.0);
    Photon photon;
    photon.position = track.start + track.direction * s;
    photon.direction =
        track.direction * cos_t + (basis.e1 * std::cos(phi) + basis.e2 * std::sin(phi)) * sin_t;
    photon.time_ns =
        track.time_ns + (track.beta > 0.0 ? s / (track.beta * speed_of_light_m_per_ns) : 0.0);
    photon.weight = arc.half / pi;
    const Uniforms wavelength_u = gstrickd::uniform4(key(Stream::wavelength));
    photon.wavelength_nm = std::visit(
        [&](const auto &spectrum) { return spectrum.sample(wavelength_u); }, track.spectrum);
    out.push_back(photon);
  }
}

// Expected sum of photon weights of one track (midpoint rule over the segment). Equals
// photons times the mean of the arc fractions along the segment.
[[nodiscard]] inline double expected_track_weight(const Track &track, std::size_t steps = 20000) {
  validate(track);
  if (track.length_m <= 0.0 || track.photons <= 0.0)
    return 0.0;
  double sum = 0.0;
  for (std::size_t k = 0; k < steps; ++k) {
    const double s = track.length_m * (static_cast<double>(k) + 0.5) / steps;
    const TrackArc arc = track_arc(track, s);
    if (arc.hit)
      sum += arc.half / pi;
  }
  return track.photons * sum / static_cast<double>(steps);
}

} // namespace dlamb
