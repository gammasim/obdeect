#pragma once

// Light sources: a geometry (point, line, plane wave, sky background) with angular, spectral and
// pulse profiles, an optional sampling target, photon counts and validation
// (docs/spec/source-core.md). Per-photon sampling is in emit.hpp.

#include "obdeect/sources/density.hpp"
#include "obdeect/sources/streams.hpp"
#include "obdeect/sources/vec3.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/random.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>

namespace obdeect::sources {

inline constexpr double max_magnitude = 1e30;
inline constexpr double max_photon_count = 9007199254740992.0; // 2^53
inline constexpr double min_positive = 1e-30;

// Columns written for the photons of a source (no zem; the run writer adds it in the sky frame).
inline constexpr std::uint32_t source_columns =
    obdeect::stream::column_bit(obdeect::stream::Column::x) |
    obdeect::stream::column_bit(obdeect::stream::Column::y) |
    obdeect::stream::column_bit(obdeect::stream::Column::z) |
    obdeect::stream::column_bit(obdeect::stream::Column::cx) |
    obdeect::stream::column_bit(obdeect::stream::Column::cy) |
    obdeect::stream::column_bit(obdeect::stream::Column::cz) |
    obdeect::stream::column_bit(obdeect::stream::Column::t) |
    obdeect::stream::column_bit(obdeect::stream::Column::w) |
    obdeect::stream::column_bit(obdeect::stream::Column::wavelength);

struct Photon {
  Vec3 position{};
  Vec3 direction{};
  double time_ns{};
  double weight{1.0};
  double wavelength_nm{};
};

// Emission from one point; the angular profile is about `axis`.
struct Point {
  Vec3 position{};
  Vec3 axis{0.0, 0.0, 1.0};
  Angular angular{Isotropic{}};
};

// Emission along a segment from `start` in `direction`; the angular profile is about the
// direction (a fixed angle is a Cherenkov cone). Points along the line follow `profile`
// (relative photons per metre against the distance from the start, empty = uniform). The
// emission time grows by s / (beta c) along the line (beta = 0: all at the pulse time).
struct Line {
  Vec3 start{};
  Vec3 direction{0.0, 0.0, -1.0};
  double length_m{};
  double beta{1.0};
  Angular angular{FixedAngle{0.0}};
  LinearDensity profile{};
};

// Parallel beam (infinite distance) or point source at a finite distance, through a disc at
// height disc_z_m; rays start on the plane z = start_z_m. Photons travel along -towards_source:
// towards_source.z > 0 is a beam from above (downward photons), z < 0 a beam from below (upward
// photons). The start plane lies on the source side of the disc. Horizontal beams (z = 0) need
// the sphere representation (not yet supported).
struct PlaneWave {
  Vec3 towards_source{0.0, 0.0, 1.0};
  double distance_m{std::numeric_limits<double>::infinity()};
  double disc_radius_m{};
  double disc_z_m{};
  double start_z_m{};
};

// Diffuse light: photons start on a sphere and travel inwards, directions uniform in solid angle
// within viewcone_deg of `axis`.
struct SkyBackground {
  Vec3 centre{};
  double radius_m{};
  Vec3 axis{0.0, 0.0, -1.0};
  double viewcone_deg{90.0};
};

using Geometry = std::variant<Point, Line, PlaneWave, SkyBackground>;

struct DiscTarget {
  Vec3 centre{};
  Vec3 normal{0.0, 0.0, 1.0};
  double radius_m{};
};

struct SphereTarget {
  Vec3 centre{};
  double radius_m{};
};

// Resolved by the run writer from the telescopes (resolve_auto_target in run.hpp).
struct AutoTarget {};

using Target = std::variant<std::monostate, DiscTarget, SphereTarget, AutoTarget>;

struct Source {
  Geometry geometry{Point{}};
  Spectrum spectrum{Fixed{400.0}};
  Pulse pulse{Fixed{0.0}};
  double photons{}; // mean photons: into 4 pi (point), along the line, on the disc, on the sphere
  Target target{};
};

// Draws of one photon for one stream.
inline Uniforms philox_draw(std::uint64_t seed, std::uint32_t event, std::uint32_t telescope,
                            std::uint64_t index, Stream stream) {
  return obdeect::stream::uniform4(
      {seed, stage_sources, event, telescope, index, static_cast<std::uint32_t>(stream)});
}

// Seed of emission i of a run seed (splitmix64 of base + golden * (i + 1)).
inline std::uint64_t source_seed(std::uint64_t base, std::uint64_t i) {
  std::uint64_t z = base + 0x9E3779B97F4A7C15ULL * (i + 1);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

// Magnitude limits (finite, |x| <= 1e30; positive values also >= 1e-30).
inline bool within_limits(double x) { return std::isfinite(x) && std::abs(x) <= max_magnitude; }
inline bool within_limits(const Vec3 &v) {
  return within_limits(v.x) && within_limits(v.y) && within_limits(v.z);
}
inline bool positive_within_limits(double x) {
  return std::isfinite(x) && x >= min_positive && x <= max_magnitude;
}
inline bool is_unit(const Vec3 &v) { return is_finite(v) && std::abs(norm(v) - 1.0) <= 1e-9; }

namespace detail {

[[noreturn]] inline void reject(const std::string &what) { throw std::invalid_argument(what); }

inline void require(bool ok, const char *what) {
  if (!ok)
    reject(what);
}

} // namespace detail

// Smallest and largest wavelength a spectrum can produce.
inline std::pair<double, double> spectrum_range(const Spectrum &spectrum) {
  if (const auto *f = std::get_if<Fixed>(&spectrum))
    return {f->value, f->value};
  if (const auto *r = std::get_if<UniformRange>(&spectrum))
    return {r->low, r->high};
  if (const auto *i = std::get_if<InverseSquare>(&spectrum))
    return {i->low, i->high};
  const auto &t = std::get<Tabulated>(spectrum);
  return {t.density.lower(), t.density.upper()};
}

inline void validate_spectrum(const Spectrum &spectrum) {
  using detail::require;
  if (const auto *f = std::get_if<Fixed>(&spectrum)) {
    require(positive_within_limits(f->value), "spectrum: wavelength out of range");
  } else if (const auto *r = std::get_if<UniformRange>(&spectrum)) {
    require(positive_within_limits(r->low) && positive_within_limits(r->high) && r->low <= r->high,
            "spectrum: need 0 < low <= high within limits");
  } else if (const auto *i = std::get_if<InverseSquare>(&spectrum)) {
    require(positive_within_limits(i->low) && positive_within_limits(i->high) && i->low < i->high,
            "spectrum: need 0 < low < high within limits");
  } else if (const auto *t = std::get_if<Tabulated>(&spectrum)) {
    require(!t->density.empty() && positive_within_limits(t->density.lower()) &&
                within_limits(t->density.upper()),
            "spectrum: table is empty or out of range");
  }
}

inline void validate_pulse(const Pulse &pulse) {
  using detail::require;
  if (const auto *f = std::get_if<Fixed>(&pulse)) {
    require(within_limits(f->value), "pulse: time out of range");
  } else if (const auto *r = std::get_if<UniformRange>(&pulse)) {
    require(within_limits(r->low) && within_limits(r->high) && r->low <= r->high,
            "pulse: need low <= high within limits");
  } else if (const auto *g = std::get_if<GaussianShape>(&pulse)) {
    require(within_limits(g->mean_value) && within_limits(g->sigma) && g->sigma >= 0.0,
            "pulse: Gaussian mean or sigma out of range");
  } else if (const auto *e = std::get_if<Exponential>(&pulse)) {
    require(within_limits(e->mean_value) && within_limits(e->tau) && e->tau >= 0.0,
            "pulse: exponential mean or decay out of range");
  } else if (const auto *c = std::get_if<ComplexPulse>(&pulse)) {
    require(within_limits(c->flat_width) && c->flat_width >= 0.0 && within_limits(c->gauss_rms) &&
                c->gauss_rms >= 0.0 && within_limits(c->decay) && c->decay >= 0.0 &&
                within_limits(c->delay),
            "pulse: complex pulse widths must be finite and >= 0, delay finite");
  } else if (const auto *t = std::get_if<Tabulated>(&pulse)) {
    require(!t->density.empty() && within_limits(t->density.lower()) &&
                within_limits(t->density.upper()),
            "pulse: table is empty or out of range");
  }
}

inline void validate_angular(const Angular &angular, const char *what) {
  using detail::require;
  if (const auto *a = std::get_if<FixedAngle>(&angular)) {
    require(a->u0 >= 0.0 && a->u0 <= 2.0, what);
  } else if (const auto *c = std::get_if<FilledCone>(&angular)) {
    require(c->u0 > 0.0 && c->u0 <= 2.0, what);
  } else if (const auto *t = std::get_if<TabulatedAngular>(&angular)) {
    require(!t->density.empty() && t->density.lower() >= 0.0 && t->density.upper() <= 2.0 + 1e-12,
            what);
  }
}

namespace detail {

inline void validate_point(const Point &p, const Target &target) {
  require(within_limits(p.position), "point: position out of range");
  require(is_unit(p.axis), "point: axis must be a unit vector");
  validate_angular(p.angular, "point: angular profile must lie in u = 1 - cos(theta) in [0, 2]");
  if (const auto *disc = std::get_if<DiscTarget>(&target)) {
    require(within_limits(disc->centre) && is_unit(disc->normal) &&
                positive_within_limits(disc->radius_m),
            "disc target: centre, normal or radius out of range");
    require(!std::holds_alternative<FixedAngle>(p.angular),
            "a disc target needs a density: fixed-angle emission is not supported");
    require(std::abs(dot(disc->centre - p.position, disc->normal)) > 0.0,
            "disc target: the point lies in the plane of the disc");
  } else if (const auto *sphere = std::get_if<SphereTarget>(&target)) {
    require(within_limits(sphere->centre) && positive_within_limits(sphere->radius_m),
            "sphere target: centre or radius out of range");
    require(std::holds_alternative<Isotropic>(p.angular),
            "a sphere target needs an isotropic point");
    require(norm(sphere->centre - p.position) > sphere->radius_m,
            "sphere target: the point lies inside the sphere");
  }
}

inline void validate_line(const Line &l, const Target &target) {
  require(within_limits(l.start), "line: start out of range");
  require(is_unit(l.direction), "line: direction must be a unit vector");
  require(within_limits(l.length_m) && l.length_m >= 0.0, "line: length out of range");
  require(l.beta >= 0.0 && l.beta <= 1.0, "line: beta must lie in [0, 1]");
  validate_angular(l.angular, "line: angular profile must lie in u = 1 - cos(theta) in [0, 2]");
  if (!l.profile.empty())
    require(l.profile.lower() >= 0.0 && l.profile.upper() <= l.length_m * (1.0 + 1e-12),
            "line: profile must lie within [0, length]");
  require(!std::holds_alternative<DiscTarget>(target), "a line takes no disc target");
  if (const auto *sphere = std::get_if<SphereTarget>(&target)) {
    require(within_limits(sphere->centre) && positive_within_limits(sphere->radius_m),
            "sphere target: centre or radius out of range");
    require(std::holds_alternative<FixedAngle>(l.angular) ||
                std::holds_alternative<Isotropic>(l.angular),
            "a line with a sphere target needs a fixed angle or isotropic emission");
  }
}

inline void validate_plane_wave(const PlaneWave &w, const Target &target) {
  require(std::holds_alternative<std::monostate>(target), "a plane wave takes no target");
  require(is_unit(w.towards_source), "plane wave: towards_source must be a unit vector");
  require(w.towards_source.z != 0.0,
          "plane wave: a horizontal towards_source (z = 0) is not supported yet: horizontal plane "
          "waves need the sphere representation (planned)");
  require(w.distance_m == std::numeric_limits<double>::infinity() ||
              positive_within_limits(w.distance_m),
          "plane wave: distance must be infinite or positive within limits");
  require(positive_within_limits(w.disc_radius_m), "plane wave: disc radius out of range");
  require(within_limits(w.disc_z_m) && within_limits(w.start_z_m),
          "plane wave: disc or start height out of range");
  // The photons move from the start plane to the disc: the plane must be on the source side.
  require((w.start_z_m - w.disc_z_m) * w.towards_source.z >= 0.0,
          "plane wave: the start plane z = start_z must lie on the source side of the disc (above "
          "it for towards_source.z > 0, below it for towards_source.z < 0)");
}

inline void validate_sky_background(const SkyBackground &b, const Target &target) {
  require(std::holds_alternative<std::monostate>(target), "a sky background takes no target");
  require(within_limits(b.centre), "sky background: centre out of range");
  require(positive_within_limits(b.radius_m), "sky background: radius out of range");
  require(is_unit(b.axis), "sky background: axis must be a unit vector");
  require(b.viewcone_deg > 0.0 && b.viewcone_deg <= 180.0,
          "sky background: viewcone must lie in (0, 180] deg");
}

} // namespace detail

// Throws std::invalid_argument for any invalid source. An AutoTarget must be resolved first.
inline void validate(const Source &s) {
  using detail::require;
  require(std::isfinite(s.photons) && s.photons >= 0.0 && s.photons <= max_photon_count,
          "photons must be finite, >= 0 and at most 2^53");
  require(!std::holds_alternative<AutoTarget>(s.target),
          "an automatic target must be resolved from the telescopes before emission");
  validate_spectrum(s.spectrum);
  validate_pulse(s.pulse);
  std::visit(
      [&](const auto &g) {
        using G = std::decay_t<decltype(g)>;
        if constexpr (std::is_same_v<G, Point>)
          detail::validate_point(g, s.target);
        else if constexpr (std::is_same_v<G, Line>)
          detail::validate_line(g, s.target);
        else if constexpr (std::is_same_v<G, PlaneWave>)
          detail::validate_plane_wave(g, s.target);
        else
          detail::validate_sky_background(g, s.target);
      },
      s.geometry);
}

// Fraction of 4 pi that a sphere covers from a point (1 inside the sphere).
inline double sphere_fraction(const Vec3 &point, const SphereTarget &sphere) {
  const double d = norm(sphere.centre - point);
  if (d <= sphere.radius_m)
    return 1.0;
  const double ratio = sphere.radius_m / d;
  return 0.5 * (1.0 - std::sqrt(1.0 - ratio * ratio));
}

// Largest sphere fraction along a line (at the point closest to the centre).
inline double line_sphere_fraction_max(const Line &l, const SphereTarget &sphere) {
  const double s = std::clamp(dot(sphere.centre - l.start, l.direction), 0.0, l.length_m);
  return sphere_fraction(l.start + l.direction * s, sphere);
}

// Azimuths (relative to the basis of the line direction) that reach a sphere: centre +- half;
// fraction half / pi.
struct Arc {
  bool hit{};
  double centre{};
  double half{};
};

namespace detail {

// Arc of a cone of cos_theta about `direction` from `point` towards a sphere.
inline Arc arc_towards(const Vec3 &point, const Vec3 &direction, double cos_theta,
                       const SphereTarget &sphere) {
  const Vec3 to_centre = sphere.centre - point;
  const double distance = norm(to_centre);
  if (distance <= sphere.radius_m)
    return {true, 0.0, pi};
  const double ratio = sphere.radius_m / distance;
  const double cos_alpha = std::sqrt(1.0 - ratio * ratio);
  const Vec3 unit = to_centre * (1.0 / distance);
  const double cos_b = dot(direction, unit);
  const Vec3 perp = unit - direction * cos_b;
  const double sin_b = norm(perp);
  const Basis basis = orthonormal_basis(direction);
  const double phi0 = std::atan2(dot(perp, basis.e2), dot(perp, basis.e1));
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
  const double a = cos_alpha - cos_theta * cos_b;
  const double b = sin_theta * sin_b;
  if (b <= 0.0)
    return a < 0.0 ? Arc{true, 0.0, pi} : Arc{};
  const double x = a / b;
  if (x >= 1.0)
    return {};
  if (x <= -1.0)
    return {true, 0.0, pi};
  return {true, phi0, std::acos(x)};
}

} // namespace detail

// Number of points of the scan for the reference arc fraction of a fixed-angle line.
inline constexpr int line_reference_points = 1025;

// Reference fraction of a line with a sphere target: candidates = photons x reference, weight =
// fraction at the point / reference. Isotropic: the largest fraction of 4 pi (closed form).
// Fixed angle: the largest arc fraction (half / pi) of a uniform scan along the line, so weights
// stay near or below 1; any positive reference gives the same expected weights.
inline double line_reference_fraction(const Line &l, const SphereTarget &sphere) {
  if (std::holds_alternative<Isotropic>(l.angular))
    return line_sphere_fraction_max(l, sphere);
  const double cos_theta = 1.0 - std::get<FixedAngle>(l.angular).u0;
  double best = 0.0;
  for (int j = 0; j < line_reference_points; ++j) {
    const double s = l.length_m * j / (line_reference_points - 1);
    best = std::max(
        best,
        detail::arc_towards(l.start + l.direction * s, l.direction, cos_theta, sphere).half / pi);
  }
  return best;
}

// Mean number of candidate photons (sampled rows before line misses are dropped).
inline double expected_count(const Source &s) {
  validate(s);
  double mean = s.photons;
  if (const auto *p = std::get_if<Point>(&s.geometry)) {
    if (const auto *disc = std::get_if<DiscTarget>(&s.target)) {
      const double h = std::abs(dot(disc->centre - p->position, disc->normal));
      const double per_photon = intensity_maximum(p->angular) / (two_pi * integral(p->angular)) *
                                pi * (disc->radius_m / h) * (disc->radius_m / h);
      mean = s.photons * per_photon;
    } else if (const auto *sphere = std::get_if<SphereTarget>(&s.target)) {
      const double ratio = sphere->radius_m / norm(sphere->centre - p->position);
      mean = s.photons * (1.0 - std::sqrt(1.0 - ratio * ratio)) / 2.0;
    }
  } else if (const auto *l = std::get_if<Line>(&s.geometry)) {
    const auto *sphere = std::get_if<SphereTarget>(&s.target);
    if (sphere)
      mean = s.photons * line_reference_fraction(*l, *sphere);
  }
  if (!std::isfinite(mean) || mean > max_photon_count)
    throw std::invalid_argument("the expected photon count is not finite or exceeds 2^53");
  return mean;
}

// Rounds a mean to a count: floor plus one with the probability of the fraction.
inline std::uint64_t round_count(double mean, std::uint64_t seed, std::uint32_t event,
                                 std::uint32_t telescope) {
  const double whole = std::floor(mean);
  const double u = philox_draw(seed, event, telescope, 0, Stream::count)[0];
  return static_cast<std::uint64_t>(whole) + (u < mean - whole ? 1U : 0U);
}

inline std::uint64_t photon_count(const Source &s, std::uint64_t seed, std::uint32_t event,
                                  std::uint32_t telescope) {
  return round_count(expected_count(s), seed, event, telescope);
}

} // namespace obdeect::sources
