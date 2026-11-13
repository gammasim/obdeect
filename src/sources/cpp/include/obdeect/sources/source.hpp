#pragma once

// Light sources: emitters and stars with angular, spectral and pulse profiles, optional
// targets, photon counts and per-photon emission (docs/spec/source-core.md).

#include "obdeect/sources/density.hpp"
#include "obdeect/sources/streams.hpp"
#include "obdeect/sources/vec3.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/random.hpp"

#include <cmath>
#include <cstddef>
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

// Columns written for the photons of a source (no zem; the stream writer adds it).
inline constexpr std::uint32_t source_columns =
    obdeect::stream::column_bit(obdeect::stream::Column::x) | obdeect::stream::column_bit(obdeect::stream::Column::y) |
    obdeect::stream::column_bit(obdeect::stream::Column::z) | obdeect::stream::column_bit(obdeect::stream::Column::cx) |
    obdeect::stream::column_bit(obdeect::stream::Column::cy) | obdeect::stream::column_bit(obdeect::stream::Column::cz) |
    obdeect::stream::column_bit(obdeect::stream::Column::t) | obdeect::stream::column_bit(obdeect::stream::Column::w) |
    obdeect::stream::column_bit(obdeect::stream::Column::wavelength);

struct Photon {
  Vec3 position{};
  Vec3 direction{};
  double time_ns{};
  double weight{1.0};
  double wavelength_nm{};
};

struct Emitter {
  Vec3 position{};
  Vec3 pointing{0.0, 0.0, 1.0};
  Angular angular{Isotropic{}};
};

// Parallel beam (infinite distance) or point source at a finite distance, seen through a disc.
struct Star {
  Vec3 towards_source{0.0, 0.0, 1.0};
  double distance_m{std::numeric_limits<double>::infinity()};
  double disk_radius_m{};
  double disk_z_m{};
  double start_z_m{};
};

struct DiscTarget {
  Vec3 centre{};
  Vec3 normal{0.0, 0.0, 1.0};
  double radius_m{};
};

struct SphereTarget {
  Vec3 centre{};
  double radius_m{};
};

using Target = std::variant<std::monostate, DiscTarget, SphereTarget>;

struct Source {
  obdeect::stream::Frame frame{obdeect::stream::Frame::sky};
  std::variant<Emitter, Star> geometry{Emitter{}};
  Spectrum spectrum{Fixed{400.0}};
  Pulse pulse{Fixed{0.0}};
  double photons{};
  Target target{};
};

// Draws of one photon for one stream.
inline Uniforms philox_draw(std::uint64_t seed, std::uint32_t event, std::uint32_t telescope,
                            std::uint64_t index, Stream stream) {
  return obdeect::stream::uniform4(
      {seed, stage_sources, event, telescope, index, static_cast<std::uint32_t>(stream)});
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

// Rules for the spectrum of a source.
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

// Throws std::invalid_argument for any invalid source.
inline void validate(const Source &s) {
  using detail::require;
  require(std::isfinite(s.photons) && s.photons >= 0.0 && s.photons <= max_photon_count,
          "photons must be finite, >= 0 and at most 2^53");
  require(s.frame == obdeect::stream::Frame::sky || s.frame == obdeect::stream::Frame::telescope ||
              s.frame == obdeect::stream::Frame::camera,
          "a source frame must be sky, telescope or camera");
  validate_spectrum(s.spectrum);
  validate_pulse(s.pulse);

  if (const auto *star = std::get_if<Star>(&s.geometry)) {
    require(std::holds_alternative<std::monostate>(s.target), "a star takes no target");
    require(is_unit(star->towards_source) && star->towards_source.z > 0.0,
            "star: towards_source must be a unit vector with z > 0");
    require(star->distance_m == std::numeric_limits<double>::infinity() ||
                positive_within_limits(star->distance_m),
            "star: distance must be infinite or positive within limits");
    require(positive_within_limits(star->disk_radius_m), "star: disk radius out of range");
    require(within_limits(star->disk_z_m) && within_limits(star->start_z_m),
            "star: disk or start height out of range");
    return;
  }

  const auto &e = std::get<Emitter>(s.geometry);
  require(within_limits(e.position), "emitter: position out of range");
  require(is_unit(e.pointing), "emitter: pointing must be a unit vector");
  if (const auto *a = std::get_if<FixedAngle>(&e.angular)) {
    require(a->u0 >= 0.0 && a->u0 <= 2.0, "emitter: fixed angle u0 must lie in [0, 2]");
  } else if (const auto *c = std::get_if<FilledCone>(&e.angular)) {
    require(c->u0 > 0.0 && c->u0 <= 2.0, "emitter: filled cone u0 must lie in (0, 2]");
  } else if (const auto *t = std::get_if<TabulatedAngular>(&e.angular)) {
    require(!t->density.empty() && t->density.lower() >= 0.0 && t->density.upper() <= 2.0 + 1e-12,
            "emitter: angular table must be non-empty and lie in u = [0, 2]");
  }
  if (const auto *disc = std::get_if<DiscTarget>(&s.target)) {
    require(within_limits(disc->centre) && is_unit(disc->normal) &&
                positive_within_limits(disc->radius_m),
            "disc target: centre, normal or radius out of range");
    require(!std::holds_alternative<FixedAngle>(e.angular),
            "a disc target needs a density: fixed-angle emission is not supported");
    require(std::abs(dot(disc->centre - e.position, disc->normal)) > 0.0,
            "disc target: the emitter lies in the plane of the disc");
  } else if (const auto *sphere = std::get_if<SphereTarget>(&s.target)) {
    require(within_limits(sphere->centre) && positive_within_limits(sphere->radius_m),
            "sphere target: centre or radius out of range");
    require(std::holds_alternative<Isotropic>(e.angular),
            "a sphere target needs an isotropic emitter");
    require(norm(sphere->centre - e.position) > sphere->radius_m,
            "sphere target: the emitter lies inside the sphere");
  }
}

// Mean number of generated photons (candidates for a disc target).
inline double expected_count(const Source &s) {
  validate(s);
  double mean = s.photons;
  if (const auto *e = std::get_if<Emitter>(&s.geometry)) {
    if (const auto *disc = std::get_if<DiscTarget>(&s.target)) {
      const double h = std::abs(dot(disc->centre - e->position, disc->normal));
      const double per_photon = intensity_maximum(e->angular) / (two_pi * integral(e->angular)) *
                                pi * (disc->radius_m / h) * (disc->radius_m / h);
      mean = s.photons * per_photon;
    } else if (const auto *sphere = std::get_if<SphereTarget>(&s.target)) {
      const double ratio = sphere->radius_m / norm(sphere->centre - e->position);
      mean = s.photons * (1.0 - std::sqrt(1.0 - ratio * ratio)) / 2.0;
    }
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

namespace detail {

inline Vec3 cone_direction(const Vec3 &axis, double cos_theta, double azimuth) {
  const Basis basis = orthonormal_basis(axis);
  const double sin_theta = std::sqrt(std::max(0.0, (1.0 - cos_theta) * (1.0 + cos_theta)));
  return direction_about(axis, basis, cos_theta, sin_theta, azimuth);
}

} // namespace detail

// Photon `index` of a source (no validation: use emit() or validate() first).
inline Photon sample_photon(const Source &s, std::uint64_t seed, std::uint32_t event,
                            std::uint32_t telescope, std::uint64_t index) {
  Photon p;
  const double pulse_time =
      sample(s.pulse, philox_draw(seed, event, telescope, index, Stream::pulse));
  p.wavelength_nm =
      sample(s.spectrum, philox_draw(seed, event, telescope, index, Stream::wavelength));
  p.time_ns = pulse_time;
  p.weight = 1.0;

  if (const auto *star = std::get_if<Star>(&s.geometry)) {
    const Uniforms ut = philox_draw(seed, event, telescope, index, Stream::target);
    const double r = star->disk_radius_m * std::sqrt(ut[0]);
    const double phi = two_pi * ut[1];
    const Vec3 target{r * std::cos(phi), r * std::sin(phi), star->disk_z_m};
    const Vec3 centre{0.0, 0.0, star->disk_z_m};
    p.direction = -star->towards_source;
    if (std::isinf(star->distance_m)) {
      const double along = (star->start_z_m - target.z) / p.direction.z;
      p.position = target + p.direction * along;
      p.time_ns += dot(p.position - centre, p.direction) / speed_of_light_m_per_ns;
    } else {
      const Vec3 source = centre + star->towards_source * star->distance_m;
      p.direction = normalised(target - source);
      if (source.z <= star->start_z_m) {
        p.position = source;
      } else {
        const double along = (star->start_z_m - target.z) / p.direction.z;
        p.position = target + p.direction * along;
      }
      p.time_ns += (norm(p.position - source) - star->distance_m) / speed_of_light_m_per_ns;
    }
    return p;
  }

  const auto &e = std::get<Emitter>(s.geometry);
  p.position = e.position;
  const Uniforms ud = philox_draw(seed, event, telescope, index, Stream::direction);
  if (const auto *disc = std::get_if<DiscTarget>(&s.target)) {
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
    const double cosine = std::clamp(dot(p.direction, e.pointing), -1.0, 1.0);
    const double ratio = h / distance;
    p.weight = std::visit([&](const auto &a) { return a.intensity(1.0 - cosine); }, e.angular) /
               intensity_maximum(e.angular) * ratio * ratio * ratio;
  } else if (const auto *sphere = std::get_if<SphereTarget>(&s.target)) {
    const Vec3 to_centre = sphere->centre - e.position;
    const double d = norm(to_centre);
    const Vec3 axis = to_centre * (1.0 / d);
    const double ratio = sphere->radius_m / d;
    const double cos_max = std::sqrt(1.0 - ratio * ratio);
    p.direction = detail::cone_direction(axis, 1.0 - ud[0] * (1.0 - cos_max), two_pi * ud[1]);
  } else {
    const double u = angle_quantile(e.angular, ud[0]);
    p.direction = detail::cone_direction(e.pointing, 1.0 - u, two_pi * ud[1]);
  }
  return p;
}

// Photons first_index .. first_index + rows - 1 as float32 columns (source_columns).
inline void emit(const Source &s, std::uint64_t seed, std::uint32_t event, std::uint32_t telescope,
                 std::uint64_t first_index, std::size_t rows, obdeect::stream::PhotonChunk &chunk) {
  validate(s);
  if (rows > obdeect::stream::max_chunk_rows)
    throw std::invalid_argument("a chunk holds at most 65536 rows");
  chunk.resize(rows, source_columns);
  chunk.telescope_id = telescope;
  chunk.frame = s.frame;
  chunk.first_index = first_index;
  using obdeect::stream::Column;
  for (std::size_t i = 0; i < rows; ++i) {
    const Photon p = sample_photon(s, seed, event, telescope, first_index + i);
    chunk[Column::x][i] = static_cast<float>(p.position.x);
    chunk[Column::y][i] = static_cast<float>(p.position.y);
    chunk[Column::z][i] = static_cast<float>(p.position.z);
    chunk[Column::cx][i] = static_cast<float>(p.direction.x);
    chunk[Column::cy][i] = static_cast<float>(p.direction.y);
    chunk[Column::cz][i] = static_cast<float>(p.direction.z);
    chunk[Column::t][i] = static_cast<float>(p.time_ns);
    chunk[Column::w][i] = static_cast<float>(p.weight);
    chunk[Column::wavelength][i] = static_cast<float>(p.wavelength_nm);
  }
}

} // namespace obdeect::sources
