#pragma once

#include "dlamb/density.hpp"
#include "dlamb/streams.hpp"
#include "dlamb/vec3.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/random.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace dlamb {

inline constexpr double speed_of_light_m_per_ns = 0.299792458;

// A point emitter (LightEmission LightSource, obdeect Emitter). Directions follow the
// angular profile about `pointing`; positions and directions are in the source frame.
struct Emitter {
  Vec3 position{};
  Vec3 pointing{0.0, 0.0, 1.0};
  Angular angular{Isotropic{}};
};

// Parallel beam from a star (distance infinite) or a point source at a finite distance
// (obdeect PointSource). Photons land uniformly on a disc of radius disk_radius_m at
// height disk_z_m; each ray starts where it crosses the plane start_z_m.
struct Star {
  Vec3 towards_source{0.0, 0.0, 1.0};
  double distance_m{std::numeric_limits<double>::infinity()};
  double disk_radius_m{};
  double disk_z_m{};
  double start_z_m{};
};

// Importance targets: only emitted photons that reach the target are generated.
// DiscTarget keeps every candidate with a weight equal to its probability of reaching
// the disc (obdeect Emitter); it works with any density-based angular profile.
// SphereTarget samples directions in the cone that contains the sphere (LightEmission
// style); it is exact for isotropic emitters only.
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

// One light source. `photons` is the mean number of photons emitted into 4 pi for an
// emitter (before any target) or the mean number of photons on the disc of a star.
struct Source {
  gstrickd::Frame frame{gstrickd::Frame::sky};
  std::variant<Emitter, Star> geometry{Emitter{}};
  Spectrum spectrum{Fixed{400.0}};
  Pulse pulse{Fixed{0.0}};
  double photons{};
  Target target{};
};

// One emitted photon before it is stored as float32.
struct Photon {
  Vec3 position{};
  Vec3 direction{};
  double time_ns{};
  double weight{1.0};
  double wavelength_nm{};
};

inline constexpr std::uint32_t source_columns =
    gstrickd::column_bit(gstrickd::Column::x) | gstrickd::column_bit(gstrickd::Column::y) |
    gstrickd::column_bit(gstrickd::Column::z) | gstrickd::column_bit(gstrickd::Column::cx) |
    gstrickd::column_bit(gstrickd::Column::cy) | gstrickd::column_bit(gstrickd::Column::cz) |
    gstrickd::column_bit(gstrickd::Column::t) | gstrickd::column_bit(gstrickd::Column::w) |
    gstrickd::column_bit(gstrickd::Column::wavelength);

// Largest magnitude of a length, time or wavelength. Squares of lengths and the float32
// photon columns then stay finite, and positive values are at least its inverse, so
// their inverses stay finite too.
inline constexpr double max_magnitude = 1e30;

// Largest photon count that is exact in a double (2^53). Larger means cannot be rounded
// to an integer count without loss.
inline constexpr double max_photon_count = 9007199254740992.0;

[[nodiscard]] inline bool magnitude_ok(double x) {
  return std::isfinite(x) && std::abs(x) <= max_magnitude;
}
[[nodiscard]] inline bool magnitude_ok(Vec3 a) {
  return magnitude_ok(a.x) && magnitude_ok(a.y) && magnitude_ok(a.z);
}
[[nodiscard]] inline bool positive_magnitude_ok(double x) {
  return magnitude_ok(x) && x >= 1.0 / max_magnitude;
}

// Throws std::invalid_argument for an inconsistent source.
inline void validate(const Source &source) {
  const auto fail = [](const char *message) { throw std::invalid_argument(message); };
  if (!(source.photons >= 0.0) || source.photons > max_photon_count)
    fail("photons must be finite, >= 0 and at most 2^53");
  if (source.frame != gstrickd::Frame::sky && source.frame != gstrickd::Frame::telescope &&
      source.frame != gstrickd::Frame::camera)
    fail("frame must be sky, T or C");
  std::visit(
      [&](const auto &s) {
        using S = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<S, Fixed>) {
          if (!positive_magnitude_ok(s.value))
            fail("a fixed wavelength must be in [1e-30, 1e30] nm");
        } else if constexpr (std::is_same_v<S, UniformRange>) {
          if (!positive_magnitude_ok(s.low) || !positive_magnitude_ok(s.high) || !(s.high >= s.low))
            fail("a uniform wavelength range needs 0 < low <= high, within [1e-30, 1e30] nm");
        } else if constexpr (std::is_same_v<S, InverseSquare>) {
          if (!positive_magnitude_ok(s.low) || !positive_magnitude_ok(s.high) || !(s.high > s.low))
            fail(
                "an inverse-square wavelength range needs 0 < low < high, within [1e-30, 1e30] nm");
        } else if constexpr (std::is_same_v<S, Tabulated>) {
          if (s.density.empty() || !positive_magnitude_ok(s.density.lower()) ||
              !magnitude_ok(s.density.upper()))
            fail("a tabulated wavelength density needs nodes in [1e-30, 1e30] nm");
        }
      },
      source.spectrum);
  std::visit(
      [&](const auto &s) {
        using S = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<S, Fixed>) {
          if (!magnitude_ok(s.value))
            fail("a fixed time must be finite and within 1e30 ns");
        } else if constexpr (std::is_same_v<S, UniformRange>) {
          if (!magnitude_ok(s.low) || !magnitude_ok(s.high) || s.high < s.low)
            fail("a uniform pulse needs low <= high, within 1e30 ns");
        } else if constexpr (std::is_same_v<S, GaussianShape>) {
          if (!magnitude_ok(s.mean_value) || !(s.sigma >= 0.0) || !magnitude_ok(s.sigma))
            fail("a Gaussian pulse needs a mean and sigma >= 0, finite and within 1e30 ns");
        } else if constexpr (std::is_same_v<S, Exponential>) {
          if (!magnitude_ok(s.mean_value) || !(s.tau >= 0.0) || !magnitude_ok(s.tau))
            fail("an exponential pulse needs a mean and tau >= 0, finite and within 1e30 ns");
        } else if constexpr (std::is_same_v<S, ComplexPulse>) {
          if (!(s.flat_width >= 0.0) || !(s.gauss_rms >= 0.0) || !(s.decay >= 0.0) ||
              !magnitude_ok(s.flat_width) || !magnitude_ok(s.gauss_rms) || !magnitude_ok(s.decay) ||
              !magnitude_ok(s.delay))
            fail("a complex pulse needs widths >= 0 and a delay, all finite and within 1e30 ns");
        } else if constexpr (std::is_same_v<S, Tabulated>) {
          if (s.density.empty() || !magnitude_ok(s.density.lower()) ||
              !magnitude_ok(s.density.upper()))
            fail("a tabulated pulse needs nodes within 1e30 ns");
        }
      },
      source.pulse);

  if (const auto *star = std::get_if<Star>(&source.geometry)) {
    if (!std::holds_alternative<std::monostate>(source.target))
      fail("a star takes no target (its disc is the target)");
    const double unit = norm(star->towards_source);
    if (!is_finite(star->towards_source) || std::abs(unit - 1.0) > 1e-9 ||
        star->towards_source.z <= 0.0)
      fail("towards_source must be a finite unit vector pointing upwards");
    // Infinite distance is a plane wave; a finite distance must be a valid length.
    const bool plane_wave = std::isinf(star->distance_m) && star->distance_m > 0.0;
    if (!plane_wave && !positive_magnitude_ok(star->distance_m))
      fail("star distance must be positive and finite (or +inf for a plane wave)");
    if (!positive_magnitude_ok(star->disk_radius_m) || !magnitude_ok(star->disk_z_m) ||
        !magnitude_ok(star->start_z_m))
      fail("star disc radius must be positive; heights finite and within 1e30 m");
    return;
  }

  const auto &emitter = std::get<Emitter>(source.geometry);
  if (!magnitude_ok(emitter.position) || !is_finite(emitter.pointing) ||
      std::abs(norm(emitter.pointing) - 1.0) > 1e-9)
    fail("an emitter needs a position within 1e30 m and a unit pointing vector");
  std::visit(
      [&](const auto &a) {
        using A = std::decay_t<decltype(a)>;
        if constexpr (std::is_same_v<A, FixedAngle>) {
          if (!(a.u0 >= 0.0 && a.u0 <= 2.0))
            fail("a fixed angle must lie in [0, 180] deg");
        } else if constexpr (std::is_same_v<A, FilledCone>) {
          if (!(a.u0 > 0.0 && a.u0 <= 2.0))
            fail("a filled cone needs a half angle in (0, 180] deg");
        } else if constexpr (std::is_same_v<A, TabulatedAngular>) {
          if (a.density.empty() || !(a.density.lower() >= 0.0) || a.density.upper() > 2.0 + 1e-12)
            fail("a tabulated angular profile needs nodes with 0 <= 1 - cos(theta) <= 2");
        }
      },
      emitter.angular);
  if (const auto *disc = std::get_if<DiscTarget>(&source.target)) {
    if (!magnitude_ok(disc->centre) || !is_finite(disc->normal) ||
        std::abs(norm(disc->normal) - 1.0) > 1e-9 || !positive_magnitude_ok(disc->radius_m))
      fail("a disc target needs a centre within 1e30 m, a unit normal and a positive radius");
    if (std::holds_alternative<FixedAngle>(emitter.angular))
      fail("a fixed-angle emitter cannot target a disc (it has no density)");
    if (std::abs(dot(disc->centre - emitter.position, disc->normal)) <= 0.0)
      fail("the emitter must not lie in the plane of the disc");
  }
  if (const auto *sphere = std::get_if<SphereTarget>(&source.target)) {
    if (!magnitude_ok(sphere->centre) || !positive_magnitude_ok(sphere->radius_m))
      fail("a sphere target needs a centre within 1e30 m and a positive radius");
    if (!std::holds_alternative<Isotropic>(emitter.angular))
      fail("a sphere target needs an isotropic emitter");
    if (!(norm(sphere->centre - emitter.position) > sphere->radius_m))
      fail("the emitter must lie outside the sphere");
  }
}

// Mean candidate count per emitted photon of a disc target: the expected number of
// photons reaching the disc is photons * candidates_per_photon * mean weight. Equals
// obdeect Emitter::candidates_per_photon.
[[nodiscard]] inline double disc_candidates_per_photon(const Emitter &emitter,
                                                       const DiscTarget &disc) {
  const double h = std::abs(dot(disc.centre - emitter.position, disc.normal));
  const double per_steradian =
      intensity_maximum(emitter.angular) / (two_pi * integral(emitter.angular));
  return per_steradian * pi * disc.radius_m * disc.radius_m / (h * h);
}

// Cosine of the half angle of the cone that contains a sphere target.
[[nodiscard]] inline double sphere_cos_max(const Emitter &emitter, const SphereTarget &sphere) {
  const double ratio = sphere.radius_m / norm(sphere.centre - emitter.position);
  return std::sqrt(1.0 - ratio * ratio);
}

namespace detail {

inline double unchecked_expected_count(const Source &source) {
  if (std::holds_alternative<Star>(source.geometry))
    return source.photons;
  const auto &emitter = std::get<Emitter>(source.geometry);
  if (const auto *disc = std::get_if<DiscTarget>(&source.target))
    return source.photons * disc_candidates_per_photon(emitter, *disc);
  if (const auto *sphere = std::get_if<SphereTarget>(&source.target))
    return source.photons * 0.5 * (1.0 - sphere_cos_max(emitter, *sphere));
  return source.photons;
}

} // namespace detail

// Mean number of generated photons (candidates for a disc target). Throws when it is not
// finite or exceeds 2^53 (a near-degenerate disc target can overflow it).
[[nodiscard]] inline double expected_count(const Source &source) {
  validate(source);
  const double mean = detail::unchecked_expected_count(source);
  if (!(mean >= 0.0) || mean > max_photon_count)
    throw std::invalid_argument("the expected photon count must be finite and at most 2^53");
  return mean;
}

// Number of photons of one source in one (event, telescope): the expected count rounded
// up with probability equal to its fractional part (LightEmission rounds bunches so).
[[nodiscard]] inline std::uint64_t photon_count(const Source &source, std::uint64_t seed,
                                                std::uint32_t event, std::uint32_t telescope) {
  const double mean = expected_count(source);
  const double whole = std::floor(mean);
  const auto u = gstrickd::uniform4(
      {seed, stage_sources, event, telescope, 0, static_cast<std::uint32_t>(Stream::count)});
  return static_cast<std::uint64_t>(whole) + (u[0] < mean - whole ? 1U : 0U);
}

namespace detail {

// Emission angle of an emitter, u = 1 - cos(theta), restricted to u <= umax. The
// emitted fraction below umax maps the uniform deviate onto the profile.
[[nodiscard]] inline double angle_u(const Angular &angular, double deviate, double umax) {
  const double fraction = std::visit([umax](const auto &a) { return a.cdf(umax); }, angular);
  return angle_quantile(angular, deviate * fraction);
}

// Unit vector at polar angle cos_theta (u = 1 - cos_theta) and azimuth phi around axis.
[[nodiscard]] inline Vec3 around_axis(Vec3 axis, double cos_theta, double phi) {
  const Basis b = orthonormal_basis(axis);
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
  return axis * cos_theta + (b.e1 * std::cos(phi) + b.e2 * std::sin(phi)) * sin_theta;
}

inline void star_photon(const Star &star, const Uniforms &u, Photon &photon) {
  const double r = star.disk_radius_m * std::sqrt(u[0]);
  const double phi = two_pi * u[1];
  const Vec3 centre{0.0, 0.0, star.disk_z_m};
  const Vec3 target{r * std::cos(phi), r * std::sin(phi), star.disk_z_m};
  if (std::isinf(star.distance_m)) {
    // Plane wave: the wavefront through the disc centre is at t = 0.
    const Vec3 direction = -star.towards_source;
    const Vec3 start = target + direction * ((star.start_z_m - target.z) / direction.z);
    photon.position = start;
    photon.direction = direction;
    photon.time_ns += dot(start - centre, direction) / speed_of_light_m_per_ns;
    return;
  }
  // Point source: emits from its own position if that is below the start plane.
  const Vec3 position = centre + star.towards_source * star.distance_m;
  const Vec3 direction = normalised(target - position);
  const Vec3 start = position.z <= star.start_z_m
                         ? position
                         : target + direction * ((star.start_z_m - target.z) / direction.z);
  photon.position = start;
  photon.direction = direction;
  photon.time_ns += (norm(start - position) - star.distance_m) / speed_of_light_m_per_ns;
}

inline void emitter_photon(const Emitter &emitter, const Target &target,
                           const Uniforms &direction_u, const Uniforms &target_u, Photon &photon) {
  photon.position = emitter.position;
  if (const auto *disc = std::get_if<DiscTarget>(&target)) {
    // Uniform point on the disc; the photon weight is I(theta) / I_max (h / r)^3, the
    // probability that the emitted photon reaches the disc (obdeect Emitter).
    const Basis b = orthonormal_basis(disc->normal);
    const double r = disc->radius_m * std::sqrt(target_u[0]);
    const double phi = two_pi * target_u[1];
    const Vec3 point = disc->centre + b.e1 * (r * std::cos(phi)) + b.e2 * (r * std::sin(phi));
    const Vec3 offset = point - emitter.position;
    const double distance = norm(offset);
    const Vec3 direction = offset * (1.0 / distance);
    const double cosine = std::clamp(dot(direction, emitter.pointing), -1.0, 1.0);
    const double h = std::abs(dot(disc->centre - emitter.position, disc->normal));
    const double ratio = h / distance;
    const double intensity =
        std::visit([u = 1.0 - cosine](const auto &a) { return a.intensity(u); }, emitter.angular);
    photon.weight = intensity / intensity_maximum(emitter.angular) * ratio * ratio * ratio;
    photon.direction = direction;
    return;
  }
  if (const auto *sphere = std::get_if<SphereTarget>(&target)) {
    // Isotropic emission restricted to the cone that contains the sphere.
    const Vec3 offset = sphere->centre - emitter.position;
    const double cos_max = sphere_cos_max(emitter, *sphere);
    const double cos_theta = 1.0 - direction_u[0] * (1.0 - cos_max);
    photon.direction = around_axis(normalised(offset), cos_theta, two_pi * direction_u[1]);
    photon.weight = 1.0;
    return;
  }
  const double u_angle = angle_u(emitter.angular, direction_u[0], 2.0);
  photon.direction = around_axis(emitter.pointing, 1.0 - u_angle, two_pi * direction_u[1]);
  photon.weight = 1.0;
}

} // namespace detail

// Photon number `photon_index` of one source in one (event, telescope). Depends only on
// the keys, so chunks may be generated in any order.
[[nodiscard]] inline Photon sample_photon(const Source &source, std::uint64_t seed,
                                          std::uint32_t event, std::uint32_t telescope,
                                          std::uint64_t photon_index) {
  const auto key = [&](Stream stream) {
    return gstrickd::DrawKey{seed,      stage_sources, event,
                             telescope, photon_index,  static_cast<std::uint32_t>(stream)};
  };
  Photon photon;
  const Uniforms pulse_u = gstrickd::uniform4(key(Stream::pulse));
  photon.time_ns = std::visit([&](const auto &p) { return p.sample(pulse_u); }, source.pulse);
  const Uniforms wavelength_u = gstrickd::uniform4(key(Stream::wavelength));
  photon.wavelength_nm =
      std::visit([&](const auto &s) { return s.sample(wavelength_u); }, source.spectrum);
  if (const auto *star = std::get_if<Star>(&source.geometry)) {
    detail::star_photon(*star, gstrickd::uniform4(key(Stream::target)), photon);
  } else {
    detail::emitter_photon(std::get<Emitter>(source.geometry), source.target,
                           gstrickd::uniform4(key(Stream::direction)),
                           gstrickd::uniform4(key(Stream::target)), photon);
  }
  return photon;
}

// Photons first_index .. first_index + rows - 1 of one source, in the source frame, with
// columns x y z cx cy cz t w wavelength. Directions are unit vectors; cz < 0 is
// downward (propagation direction, as obdeect). Weights are 1 except for disc targets.
inline void emit(const Source &source, std::uint64_t seed, std::uint32_t event,
                 std::uint32_t telescope, std::uint64_t first_index, std::size_t rows,
                 gstrickd::PhotonChunk &out) {
  validate(source);
  if (rows > gstrickd::max_chunk_rows)
    throw std::invalid_argument("a chunk holds at most 2^16 photons");
  out.telescope_id = telescope;
  out.frame = source.frame;
  out.first_index = first_index;
  out.resize(rows, source_columns);
  auto &x = out[gstrickd::Column::x], &y = out[gstrickd::Column::y], &z = out[gstrickd::Column::z],
       &cx = out[gstrickd::Column::cx], &cy = out[gstrickd::Column::cy],
       &cz = out[gstrickd::Column::cz], &t = out[gstrickd::Column::t],
       &w = out[gstrickd::Column::w], &wl = out[gstrickd::Column::wavelength];
  for (std::size_t i = 0; i < rows; ++i) {
    const Photon p = sample_photon(source, seed, event, telescope, first_index + i);
    x[i] = static_cast<float>(p.position.x);
    y[i] = static_cast<float>(p.position.y);
    z[i] = static_cast<float>(p.position.z);
    cx[i] = static_cast<float>(p.direction.x);
    cy[i] = static_cast<float>(p.direction.y);
    cz[i] = static_cast<float>(p.direction.z);
    t[i] = static_cast<float>(p.time_ns);
    w[i] = static_cast<float>(p.weight);
    wl[i] = static_cast<float>(p.wavelength_nm);
  }
}

} // namespace dlamb
