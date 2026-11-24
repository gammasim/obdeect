#pragma once

#include "obdeect/ray_tracing/model.hpp"
#include "obdeect/ray_tracing/random.hpp"
#include "obdeect/ray_tracing/window.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <limits>
#include <mutex>
#include <optional>
#include <string_view>
#include <thread>
#include <vector>

namespace obdeect::ray_tracing {

enum class Status : std::uint8_t {
  detected,              // Reached the focal surface (the final leg, or an asphere from below).
  missed_primary,        // Never met the primary mirror.
  missed_secondary,      // Reflected on the primary of a dual mirror, then met nothing.
  missed_focal_surface,  // Met nothing after the last reflection (or from a flasher).
  blocked,               // Absorbed by an obscurer or a back face (terminal_surface names it).
  invalid_input,         // Non-finite input, wavelength <= 0 or non-unit direction.
  numerical_failure,     // An intersection did not converge.
  non_sequential,        // Retired: the reflecting face of the secondary is traced on (no longer
                         // produced).
  direct_focal,          // A plane focal surface (no extent) met from below before the last mirror.
  bounce_limit,          // kMaxReflections reflections reached.
  blocked_camera_window, // Lost at the rim or the outer ring of the camera window.
  camera_window_total_reflection, // Total internal reflection at a camera window surface.
};
inline constexpr std::size_t kStatusCount = 12;

[[nodiscard]] constexpr std::string_view to_string(Status status) {
  constexpr std::string_view names[kStatusCount] = {
      "detected",     "missed_primary", "missed_secondary",      "missed_focal_surface",
      "blocked",      "invalid_input",  "numerical_failure",     "non_sequential",
      "direct_focal", "bounce_limit",   "blocked_camera_window", "camera_window_total_reflection"};
  return names[static_cast<std::size_t>(status)];
}

// Reflections traced per photon (the first, the second, ...); a photon that would reflect
// again ends as bounce_limit at that mirror.
inline constexpr std::uint32_t kMaxReflections = 16;

inline constexpr std::uint32_t kNoSurface = 0xFFFFFFFFU;

// Result of one photon. Mirror hits not reached keep NaN coordinates and
// kNoSurface; `weight` is the input weight times the applied responses.
struct TraceResult {
  Status status{Status::invalid_input};
  std::uint32_t terminal_surface{kNoSurface};
  Vec3 position, direction;
  double path_length_m{}, time_ns{}, weight{};
  std::uint32_t primary_surface{kNoSurface}, secondary_surface{kNoSurface};
  Vec3 primary_point{nan(), nan(), nan()}, secondary_point{nan(), nan(), nan()};
  double primary_incidence_deg{nan()}, secondary_incidence_deg{nan()}, focal_incidence_deg{nan()};
  // Reflections of the photon (on either mirror).
  std::uint32_t reflections{};

  static constexpr double nan() { return std::numeric_limits<double>::quiet_NaN(); }
};

namespace detail {

constexpr double kDegree = 57.29577951308232087680;

inline double incidence_deg(const Vec3 &direction, const Vec3 &normal) {
  return std::acos(std::min(1.0, std::abs(dot(direction, normal)))) * kDegree;
}

struct MirrorHit {
  SurfaceHit hit;
  std::uint32_t facet_index{};
};

// Nearest hit on a mirror (lowest id on exact ties). `failed` reports a
// non-converged aspheric intersection.
// With `any_side` a flat or curved facet is hit from its back too (photons from below).
inline std::optional<MirrorHit> intersect_mirror(const Mirror &mirror, const Ray &ray, bool &failed,
                                                 bool any_side = false) {
  failed = false;
  if (mirror.asphere) {
    const auto &a = *mirror.asphere;
    if (!a.may_meet(ray, kMinDistanceM))
      return std::nullopt;
    // The lowest segment id at the accepted point: accept() is the last call
    // that returned true when intersect() returns, and it saw the hit point.
    std::uint32_t segment_id = a.id;
    const auto accept = [&](const Vec3 &p) {
      if (mirror.segments.empty())
        return true;
      const auto id = mirror.segment_grid.find(mirror.segments, mirror.segment_frames, p);
      if (id)
        segment_id = *id;
      return id.has_value();
    };
    auto hit = a.intersect(ray, kMinDistanceM, accept, failed);
    if (!hit)
      return std::nullopt;
    hit->id = segment_id;
    return MirrorHit{*hit, 0};
  }
  std::optional<MirrorHit> nearest;
  mirror.grid.candidates(ray, mirror.facets.size(), [&](std::size_t i) {
    const auto &f = mirror.facets[i];
    const auto hit =
        any_side ? f.intersect_any_side(ray, kMinDistanceM) : f.intersect(ray, kMinDistanceM);
    if (hit)
      if (!nearest || hit->t < nearest->hit.t ||
          (hit->t == nearest->hit.t && hit->id < nearest->hit.id))
        nearest = MirrorHit{*hit, static_cast<std::uint32_t>(i)};
  });
  return nearest;
}

// Rotation of v about the unit axis k by `angle` (Rodrigues).
inline Vec3 rotate(const Vec3 &v, const Vec3 &k, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return v * c + cross(k, v) * s + k * (dot(k, v) * (1.0 - c));
}

// Reflected direction, with the mirror's random reflection error.
inline Vec3 reflect_at(const Mirror &mirror, const MirrorHit &mh, const Vec3 &direction,
                       std::uint64_t photon_id, const RandomKey &seed, Stream stream,
                       std::uint32_t sub_photon = 0) {
  const Vec3 &n = mh.hit.normal;
  if (!mirror.scatter)
    return reflect(direction, n);
  const auto &s = *mirror.scatter;
  const auto u = uniform4(seed, photon_id, stream, sub_photon);
  const double sigma = u[0] < s.fraction2 ? s.sigma2_rad : s.sigma1_rad;
  if (sigma == 0.0)
    return reflect(direction, n);
  const double radius = sigma * std::sqrt(-2.0 * std::log(u[1]));
  const double phi = 6.283185307179586476925 * u[2];
  const double a1 = radius * std::cos(phi), a2 = radius * std::sin(phi);
  if (s.surface_slopes) {
    // Graph-surface normal (-dz/dx, -dz/dy, 1): perturb each slope by half
    // the angle, so the reflected ray deviates by about the full angle where
    // the surface is flat. A slope is not an angle: on a tilted surface the
    // effective width is smaller (docs/ray-tracing/optics.md, mirror scatter conventions).
    const Vec3 tilted{n.x / n.z + 0.5 * a1, n.y / n.z + 0.5 * a2, 1.0};
    return reflect(direction, normalized(tilted));
  }
  // Rotate the reflected ray about two axes in the mirror's tangent plane.
  // The facet tangent is perpendicular to the normal at the vertex only; on a
  // curved facet remove its component along the hit normal, so that e1, e2 and n
  // are orthonormal (the rotations need unit axes).
  Vec3 e1 = cross(std::abs(n.z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0}, n);
  if (!mirror.facets.empty()) {
    const Vec3 &u = mirror.facets[mh.facet_index].u;
    const Vec3 projected = u - n * dot(u, n);
    // Near a steep cap edge the normal can align with u; keep the generic axis then.
    if (dot(projected, projected) > 1e-12)
      e1 = projected;
  }
  e1 = normalized(e1);
  const Vec3 e2 = cross(n, e1);
  return normalized(rotate(rotate(reflect(direction, n), e2, a1), e1, a2));
}

// Does the mirror's reflectivity depend on the incidence angle?
inline bool angle_matters(const Mirror &mirror) {
  return mirror.reflectivity && mirror.reflectivity->incidence_deg.size() > 1;
}

inline double mirror_response(const Mirror &mirror, const MirrorHit &mh, double wavelength_nm,
                              double incidence) {
  double r = mirror.reflectivity ? mirror.reflectivity->at(wavelength_nm, incidence) : 1.0;
  if (mirror.degradation) {
    const Vec3 &p = mh.hit.point;
    if (mirror.facets.empty()) {
      r *= mirror.degradation->at(p.x, p.y);
    } else {
      const auto &f = mirror.facets[mh.facet_index];
      r *= mirror.degradation->at(dot(p - f.centre, f.map_u), dot(p - f.centre, f.map_v));
    }
  }
  return r;
}

} // namespace detail

// Per-photon switches of trace_photon (see the comment there). `sub_photon` numbers the photons
// that share one photon_id; `from_below` traces a photon that moves up from below the optics.
struct TraceOptions {
  bool apply_responses{true};  // Reflectivities, degradation and telescope transmission.
  bool direct{false};          // Start inside the telescope on the leg to the focal surface.
  bool incidence_angles{true}; // Compute the incidence angles (NaN otherwise).
  bool classify_losses{true};  // Report a photon that misses the primary as missed_primary.
  std::uint32_t sub_photon{0}; // Photon number within a bunch (its random numbers).
  bool from_below{false};      // Traces up from below the optics.
};

// Trace one photon through the optical model: one loop for every photon.
// Each step takes the nearest hit along the ray among the primary, the secondary
// (both sides), the focal surface and the obscurers of the current leg. Ties go
// to the surfaces, in the order primary, secondary, focal surface, obscurer.
//
// Legs: `incoming` until the first reflection, then `after_primary` or
// `after_secondary` after the last reflection on that mirror. Outcomes:
// - obscurer, or the back of a mirror: blocked (terminal_surface names it);
// - the reflecting face of a mirror: its response and reflection, then the next step;
// - the focal surface on the final leg (after the last mirror): detected, from either
//   side (the sequential path, unchanged);
// - the focal surface on another leg, from its receiving side (the side the light of
//   the sequence arrives from: below for a single mirror, above for a dual one): an
//   asphere (with extent) is detected, a plane (no extent) gives direct_focal (a counted
//   loss); from the other side it is passed (the camera body obscurers give its shadow);
// - nothing: missed_primary before any reflection, missed_secondary after a
//   reflection on the primary of a dual mirror, missed_focal_surface otherwise;
// - kMaxReflections reflections reached: bounce_limit.
//
// The secondary is tested on the incoming leg only where the model blocks from
// the sky (`secondary_blocks_incoming`) or the photon comes from below. A photon
// from below meets the obscurers of every leg before its first reflection (it
// crosses them all on its way up).
//
// The random numbers of the first reflection on each mirror are the ones of the
// sequential path (Stream::primary_scatter, Stream::secondary_scatter); every
// later reflection draws from reflection_stream(k), k its reflection number.
//
// `direct` traces a photon that starts inside the telescope on the leg to the
// focal surface (e.g. from a flasher): no mirror and no telescope transmission act,
// only the focal surface and the obscurers of that leg. With `apply_responses`
// false the weight stays at its input value (geometry-only transport). With
// `incidence_angles` false the incidence angles stay NaN (saving an arccosine each)
// except where a response table depends on them. With `classify_losses` false a
// photon that misses the primary on the incoming leg is reported as missed_primary
// even where it would be blocked first; detected photons are the same either way.
// `sub_photon` numbers the photons that share one `photon_id` (the random
// counter's block); 0 gives the random numbers of a single photon.
//
// `from_below` traces a photon that moves up from a start below the optics (the
// caller decides, e.g. by the sign of cz). Its first hit decides as above.
//
// The flags are the fields of TraceOptions (defaults: a photon from the sky, with responses).
[[nodiscard]] inline TraceResult trace_photon(const OpticalModel &model, Ray ray,
                                              double wavelength_nm, double time_ns, double weight,
                                              std::uint64_t photon_id, const RandomKey &seed,
                                              const TraceOptions &options = {}) {
  using namespace detail;
  const bool apply_responses = options.apply_responses, direct = options.direct;
  const bool incidence_angles = options.incidence_angles;
  const bool classify_losses = options.classify_losses, from_below = options.from_below;
  const std::uint32_t sub_photon = options.sub_photon;
  TraceResult out;
  out.position = ray.origin;
  out.direction = ray.direction;
  const double length = norm(ray.direction);
  if (!finite(ray.origin) || !std::isfinite(length) || std::abs(length - 1.0) > 1e-6 ||
      !(wavelength_nm > 0.0) || !std::isfinite(wavelength_nm) || !std::isfinite(time_ns) ||
      !(weight >= 0.0) || !std::isfinite(weight))
    return out;
  ray.direction = ray.direction * (1.0 / length);
  out.weight = weight;
  double path = 0.0;
  std::uint32_t reflections = 0;
  // Extra optical path in the camera window glass: time of (n - n_group) per metre.
  double time_path_extra = 0.0;
  const auto finish = [&](Status status, const Vec3 &position, std::uint32_t surface) {
    out.status = status;
    out.position = position;
    out.direction = ray.direction;
    out.terminal_surface = surface;
    out.path_length_m = path;
    out.reflections = reflections;
    out.time_ns =
        time_ns + (path * model.propagation_group_index + time_path_extra) / kSpeedOfLightMPerNs;
    if (status != Status::detected)
      out.weight = 0.0;
    return out;
  };
  bool failed = false;
  const bool dual = model.secondary.has_value();
  const Leg final_leg = dual ? Leg::after_secondary : Leg::after_primary;
  // The focal surface receives on the side that the light of the sequence arrives from: the
  // light leaves the last mirror upwards (single mirror) or downwards (dual mirror, whose focal
  // surface lies below the secondary).
  const double receiving = dual ? -1.0 : 1.0;
  if (!direct && apply_responses && model.telescope_transmission)
    out.weight *= model.telescope_transmission->at(ray.direction);

  Leg leg = direct ? final_leg : Leg::incoming;
  // Camera window (after the final mirror, before the focal surface): the photon refracts
  // through the outer surface, crosses the glass and refracts out at the inner surface.
  // Returns the terminal outcome when the photon is lost in the window. Otherwise it returns
  // nullopt with `ray` set to the emerged direction (or, outside the footprint, unchanged).
  // `t_max` is the distance to the focal surface on the current ray: a window hit beyond it
  // (e.g. the far side of the sphere for a photon between window and focal surface) is ignored.
  const auto through_window = [&](double t_max) -> std::optional<TraceResult> {
    const Window &window = *model.camera_window;
    const auto outer = intersect_window_surface(window, ray, false, kMinDistanceM);
    if (!outer || !(outer->t < t_max))
      return std::nullopt;
    const double radius = std::hypot(outer->point.x, outer->point.y);
    // Outside the modeled window footprint, the path remains in air.
    if (radius > window.outer_radius)
      return std::nullopt;
    path += outer->t;
    if (radius > window.aperture_radius)
      return finish(Status::blocked_camera_window, outer->point, window.id);
    const auto transmit = [](const Vec3 &direction, const Vec3 &surface_normal, double n_from,
                             double n_to, Vec3 &out_direction, double &transmission) {
      Vec3 normal = surface_normal;
      if (dot(direction, normal) > 0.0)
        normal = normal * -1.0;
      const double cos_i = -dot(direction, normal);
      const double eta = n_from / n_to;
      const double sin2_t = eta * eta * std::max(0.0, 1.0 - cos_i * cos_i);
      if (sin2_t >= 1.0)
        return false;
      const double cos_t = std::sqrt(1.0 - sin2_t);
      const double rs = (n_from * cos_i - n_to * cos_t) / (n_from * cos_i + n_to * cos_t);
      const double rp = (n_to * cos_i - n_from * cos_t) / (n_to * cos_i + n_from * cos_t);
      transmission = 1.0 - 0.5 * (rs * rs + rp * rp);
      out_direction = normalized(direction * eta + normal * (eta * cos_i - cos_t));
      return true;
    };
    const double n = window.refractive_index.at(wavelength_nm, 0.0);
    Vec3 refracted;
    double t0 = 1.0;
    if (!transmit(ray.direction, outer->normal, 1.0, n, refracted, t0))
      return finish(Status::camera_window_total_reflection, outer->point, window.id);
    if (apply_responses)
      out.weight *= t0;
    ray = {outer->point, refracted};
    const auto inner = intersect_window_surface(window, ray, true, kMinDistanceM);
    if (!inner)
      return finish(Status::missed_focal_surface, ray.origin, kNoSurface);
    if (const auto b = model.obscurers.nearest(ray, leg, inner->t)) {
      path += b->t;
      time_path_extra += b->t * (n - model.propagation_group_index);
      return finish(Status::blocked, ray.at(b->t), b->id);
    }
    path += inner->t;
    time_path_extra += inner->t * (n - model.propagation_group_index);
    // Leaving through the non-spherical transition zone counts as a rim loss.
    if (std::hypot(inner->point.x, inner->point.y) > window.aperture_radius)
      return finish(Status::blocked_camera_window, inner->point, window.id);
    if (apply_responses) {
      const double tau = window.internal_transmittance.at(wavelength_nm, 0.0);
      out.weight *= std::pow(tau, inner->t / window.thickness);
    }
    Vec3 emerged;
    double t1 = 1.0;
    if (!transmit(ray.direction, inner->normal, n, 1.0, emerged, t1))
      return finish(Status::camera_window_total_reflection, inner->point, window.id);
    if (apply_responses)
      out.weight *= t1;
    ray = {inner->point, emerged};
    return std::nullopt;
  };
  // The mirror that reflected the ray last (its surface holds the ray's origin).
  bool last_primary = false, last_secondary = false;
  // The focal surface: the first intersection along the ray (sets `failed`).
  const auto focal_hit = [&]() -> std::optional<SurfaceHit> {
    const auto &fs = model.focal_surface;
    if (fs.asphere) {
      if (!fs.asphere->may_meet(ray, kMinDistanceM))
        return std::nullopt;
      return fs.asphere->intersect(ray, kMinDistanceM, [](const Vec3 &) { return true; }, failed);
    }
    if (ray.direction.z != 0.0) {
      const double t = (fs.plane_z - ray.origin.z) / ray.direction.z;
      if (t > kMinDistanceM)
        return SurfaceHit{t, ray.at(t), Vec3{0.0, 0.0, 1.0}, fs.id};
    }
    return std::nullopt;
  };

  for (;;) {
    // Mirror hits on both sides. Without a mirror step (direct) none is tested.
    std::optional<MirrorHit> m1, m2;
    // A mirror that just reflected the ray is not tested again when it can be crossed only once
    // (see Asphere::single_crossing): the ray starts on it and meets it nowhere else.
    if (!direct) {
      const bool skip_primary = last_primary && model.primary.asphere &&
                                model.primary.asphere->single_crossing(ray.direction);
      if (!skip_primary) {
        m1 = intersect_mirror(model.primary, ray, failed, true);
        if (failed)
          return finish(Status::numerical_failure, ray.origin, kNoSurface);
      }
      if (!classify_losses && !from_below && leg == Leg::incoming && !m1)
        return finish(Status::missed_primary, ray.origin, kNoSurface);
      const bool skip_secondary = last_secondary && model.secondary && model.secondary->asphere &&
                                  model.secondary->asphere->single_crossing(ray.direction);
      if (model.secondary && !skip_secondary &&
          (leg != Leg::incoming || model.secondary_blocks_incoming || from_below)) {
        m2 = intersect_mirror(*model.secondary, ray, failed, true);
        if (failed)
          return finish(Status::numerical_failure, ray.origin, kNoSurface);
      }
    }
    // The focal surface is met on the final leg (from either side), and on the other legs from
    // its receiving side only (from behind it is transparent: the camera body obscurers give its
    // shadow, and a focal surface between the mirrors is not a surface the light meets).
    std::optional<SurfaceHit> focal;
    if (leg == final_leg || receiving * ray.direction.z > 0.0) {
      focal = focal_hit();
      if (failed)
        return finish(Status::numerical_failure, ray.origin, kNoSurface);
      if (focal && leg != final_leg && !(receiving * dot(ray.direction, focal->normal) > 0.0))
        focal.reset();
    }

    // The nearest surface, then the obscurers strictly nearer than it.
    enum class Next : std::uint8_t { none, primary, secondary, focal, obscurer };
    Next next = Next::none;
    double t_best = std::numeric_limits<double>::infinity();
    if (m1) {
      t_best = m1->hit.t;
      next = Next::primary;
    }
    if (m2 && m2->hit.t < t_best) {
      t_best = m2->hit.t;
      next = Next::secondary;
    }
    if (focal && focal->t < t_best) {
      t_best = focal->t;
      next = Next::focal;
    }
    std::uint32_t obscurer_id = kNoSurface;
    const auto obscurers_on = [&](Leg l) {
      if (const auto b = model.obscurers.nearest(ray, l, t_best)) {
        t_best = b->t;
        obscurer_id = b->id;
        next = Next::obscurer;
      }
    };
    if (from_below && reflections == 0) {
      obscurers_on(Leg::incoming);
      obscurers_on(Leg::after_primary);
      obscurers_on(Leg::after_secondary);
    } else {
      obscurers_on(leg);
    }

    switch (next) {
    case Next::none:
      if (direct)
        return finish(Status::missed_focal_surface, ray.origin, kNoSurface);
      if (reflections == 0)
        return finish(Status::missed_primary, ray.origin, kNoSurface);
      if (leg == Leg::after_primary && dual)
        return finish(Status::missed_secondary, ray.origin, kNoSurface);
      return finish(Status::missed_focal_surface, ray.origin, kNoSurface);
    case Next::obscurer:
      path += t_best;
      return finish(Status::blocked, ray.at(t_best), obscurer_id);
    case Next::primary: {
      const MirrorHit &m = *m1;
      path += m.hit.t;
      // The first hit on each mirror is recorded (the sequential path has one per mirror).
      const bool first = out.primary_surface == kNoSurface;
      if (first) {
        out.primary_point = m.hit.point;
        out.primary_surface = m.hit.id;
      }
      // The reflecting face of the primary faces up: dot(direction, normal) < 0.
      if (dot(ray.direction, m.hit.normal) >= 0.0)
        return finish(Status::blocked, m.hit.point, m.hit.id);
      if (reflections >= kMaxReflections)
        return finish(Status::bounce_limit, m.hit.point, m.hit.id);
      const bool angles = incidence_angles || angle_matters(model.primary);
      const double incidence = angles ? incidence_deg(ray.direction, m.hit.normal) : 0.0;
      if (first && angles)
        out.primary_incidence_deg = incidence;
      if (apply_responses)
        out.weight *= mirror_response(model.primary, m, wavelength_nm, incidence);
      const Stream stream = first ? Stream::primary_scatter : reflection_stream(reflections + 1);
      ray = {m.hit.point,
             reflect_at(model.primary, m, ray.direction, photon_id, seed, stream, sub_photon)};
      leg = Leg::after_primary;
      last_primary = true;
      last_secondary = false;
      ++reflections;
      break;
    }
    case Next::secondary: {
      const MirrorHit &m = *m2;
      path += m.hit.t;
      const bool first = out.secondary_surface == kNoSurface;
      if (first) {
        out.secondary_point = m.hit.point;
        out.secondary_surface = m.hit.id;
      }
      // The secondary is a graph surface (normal +z) whose reflecting face faces down:
      // dot(direction, normal) > 0 on that face.
      if (dot(ray.direction, m.hit.normal) <= 0.0)
        return finish(Status::blocked, m.hit.point, m.hit.id);
      if (reflections >= kMaxReflections)
        return finish(Status::bounce_limit, m.hit.point, m.hit.id);
      const bool angles = incidence_angles || angle_matters(*model.secondary);
      const double incidence = angles ? incidence_deg(ray.direction, m.hit.normal) : 0.0;
      if (first && angles)
        out.secondary_incidence_deg = incidence;
      if (apply_responses)
        out.weight *= mirror_response(*model.secondary, m, wavelength_nm, incidence);
      const Stream stream = first ? Stream::secondary_scatter : reflection_stream(reflections + 1);
      ray = {m.hit.point,
             reflect_at(*model.secondary, m, ray.direction, photon_id, seed, stream, sub_photon)};
      leg = Leg::after_secondary;
      last_secondary = true;
      last_primary = false;
      ++reflections;
      break;
    }
    case Next::focal: {
      const auto &fs = model.focal_surface;
      if (leg != final_leg && !fs.asphere)
        return finish(Status::direct_focal, focal->point, fs.id);
      if (leg == final_leg && model.camera_window) {
        // The window lies between the last mirror and the focal surface. The focal surface
        // is then met again on the photon's refracted ray (the nearest obscurer before the
        // window was already taken above).
        if (auto lost = through_window(focal->t))
          return *lost;
        focal = focal_hit();
        if (failed)
          return finish(Status::numerical_failure, ray.origin, kNoSurface);
        if (!focal)
          return finish(Status::missed_focal_surface, ray.origin, kNoSurface);
        if (const auto b = model.obscurers.nearest(ray, leg, focal->t)) {
          path += b->t;
          return finish(Status::blocked, ray.at(b->t), b->id);
        }
      }
      path += focal->t;
      if (incidence_angles)
        out.focal_incidence_deg = incidence_deg(ray.direction, focal->normal);
      return finish(Status::detected, focal->point, fs.id);
    }
    }
  }
}

// Structure-of-arrays photon block. Vectors are (n, 3) row-major.
struct PhotonBlock {
  const double *position{}, *direction{}, *wavelength_nm{}, *time_ns{}, *weight{};
  const std::uint64_t *photon_id{};
  std::size_t size{};
};

// Output columns. A null pointer means "not requested": nothing is written.
// Vectors are (n, 3) row-major. `index` receives the input row of each output
// row (useful with detected_only).
struct ResultBlock {
  std::uint8_t *status{};
  std::uint32_t *terminal_surface{}, *primary_surface{}, *secondary_surface{};
  double *position{}, *direction{}, *path_length_m{}, *time_ns{}, *weight{};
  double *primary_point{}, *secondary_point{};
  double *primary_incidence_deg{}, *secondary_incidence_deg{}, *focal_incidence_deg{};
  std::uint64_t *index{};
  std::uint8_t *reflections{};
};

// Chunk size of the parallel loops: photons cost very different amounts, so threads
// take fixed chunks from a shared counter instead of one static share each.
inline constexpr std::size_t kChunkSize = 4096;

// A fixed set of worker threads that runs one job at a time; the calling
// thread works too. Owned by the caller (no global pool). A task must never
// call `run` on the same pool from inside a running job (it would deadlock).
class ThreadPool {
public:
  explicit ThreadPool(unsigned threads) {
    try {
      for (unsigned i = 1; i < std::max(1U, threads); ++i)
        workers_.emplace_back([this] { loop(); });
    } catch (...) { // Joinable threads must not be destroyed.
      {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
      }
      wake_.notify_all();
      for (auto &w : workers_)
        w.join();
      throw;
    }
  }
  ~ThreadPool() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
    }
    wake_.notify_all();
    for (auto &w : workers_)
      w.join();
  }
  ThreadPool(const ThreadPool &) = delete;
  ThreadPool &operator=(const ThreadPool &) = delete;

  [[nodiscard]] unsigned size() const { return static_cast<unsigned>(workers_.size()) + 1; }

  // Calls task(c) for every c in [0, count) and returns when all calls are done.
  // If `task` throws, the remaining calls are skipped and the first exception
  // is rethrown here once all workers have left the job.
  template <class Task> void run(std::size_t count, const Task &task) {
    if (count == 0)
      return;
    std::lock_guard<std::mutex> one_job(run_mutex_);
    Job job;
    job.count = count;
    job.call = [&task](std::size_t c) { task(c); };
    {
      std::lock_guard<std::mutex> lock(mutex_);
      job_ = &job;
      ++generation_;
    }
    wake_.notify_all();
    work(job);
    std::unique_lock<std::mutex> lock(mutex_);
    job_ = nullptr; // Workers that have not joined yet will not touch the job.
    done_.wait(lock, [&] { return job.active == 0; });
    if (job.error)
      std::rethrow_exception(job.error);
  }

private:
  struct Job {
    std::size_t count{};
    std::function<void(std::size_t)> call;
    std::atomic<std::size_t> next{0};
    std::size_t active{0}; // Workers inside work(job); guarded by mutex_.
    std::exception_ptr error;
    std::mutex error_mutex;
  };
  static void work(Job &job) {
    for (std::size_t c = job.next.fetch_add(1); c < job.count; c = job.next.fetch_add(1))
      try {
        job.call(c);
      } catch (...) {
        {
          std::lock_guard<std::mutex> lock(job.error_mutex);
          if (!job.error)
            job.error = std::current_exception();
        }
        job.next.store(job.count); // Skip the remaining chunks.
      }
  }
  void loop() {
    std::uint64_t seen = 0;
    for (;;) {
      Job *job = nullptr;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        wake_.wait(lock, [&] { return stop_ || (job_ != nullptr && generation_ != seen); });
        if (stop_)
          return;
        seen = generation_;
        job = job_;
        ++job->active;
      }
      work(*job);
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (--job->active == 0)
          done_.notify_all();
      }
    }
  }
  std::vector<std::thread> workers_;
  std::mutex mutex_, run_mutex_;
  std::condition_variable wake_, done_;
  Job *job_{nullptr};
  std::uint64_t generation_{0};
  bool stop_{false};
};

// Run body(begin, end) over [0, n) in chunks of kChunkSize items (the last one shorter) on the
// pool: the chunks are the pool's tasks, and the body must write only the items of its range.
template <class Body> void parallel_for(std::size_t n, ThreadPool &pool, const Body &body) {
  const std::size_t chunks = (n + kChunkSize - 1) / kChunkSize;
  pool.run(chunks, [&](std::size_t c) { body(c * kChunkSize, std::min(n, (c + 1) * kChunkSize)); });
}

// Threads of a pool for n items: at least 1, at most `threads`, and no more than there are chunks.
inline unsigned pool_size(std::size_t n, unsigned threads) {
  const std::size_t chunks = (n + kChunkSize - 1) / kChunkSize;
  return static_cast<unsigned>(std::max<std::size_t>(1, std::min<std::size_t>(threads, chunks)));
}

// The same on `threads` threads: a pool for this call (see pool_size).
template <class Body> void parallel_for(std::size_t n, unsigned threads, const Body &body) {
  ThreadPool pool(pool_size(n, threads));
  parallel_for(n, pool, body);
}

namespace detail {

// An output column for compact_chunks: its data and the bytes of one row.
struct Column {
  void *data;
  std::size_t bytes;
};

// Add a column of `per_row` elements per row; a null pointer is not requested.
template <class T> void add_column(std::vector<Column> &columns, T *pointer, std::size_t per_row) {
  if (pointer)
    columns.push_back({pointer, per_row * sizeof(T)});
}

// Chunk c wrote found[c] rows from row chunk_start[c] on (starts increasing);
// move them down in order so that all rows are contiguous, and return their
// number. Moving in order never overwrites rows not yet moved. The columns are
// independent: `run(count, task)` calls task(k) for every column k, on any
// threads.
template <class Run>
inline std::size_t compact_chunks(const std::vector<Column> &columns,
                                  const std::vector<std::size_t> &chunk_start,
                                  const std::vector<std::size_t> &found, const Run &run) {
  std::vector<std::size_t> destination(found.size());
  std::size_t total = 0;
  for (std::size_t c = 0; c < found.size(); ++c) {
    destination[c] = total;
    total += found[c];
  }
  run(columns.size(), [&](std::size_t k) {
    auto *bytes = static_cast<char *>(columns[k].data);
    for (std::size_t c = 0; c < found.size(); ++c)
      if (destination[c] != chunk_start[c])
        std::memmove(bytes + destination[c] * columns[k].bytes,
                     bytes + chunk_start[c] * columns[k].bytes, found[c] * columns[k].bytes);
  });
  return total;
}

// One block to trace: the per-chunk work and the final compaction, shared by
// every way of running the chunks (threads, pool, several blocks in one job).
class BlockTask {
public:
  BlockTask(const OpticalModel &model, const PhotonBlock &in, const ResultBlock &out,
            std::uint64_t seed, bool apply_responses, bool direct, bool detected_only,
            bool classify_losses)
      : model_(model), in_(in), out_(out), seed_(seed), apply_responses_(apply_responses),
        direct_(direct), detected_only_(detected_only), classify_losses_(classify_losses),
        wants_angles_(out.primary_incidence_deg || out.secondary_incidence_deg ||
                      out.focal_incidence_deg),
        found_(detected_only ? (in.size + kChunkSize - 1) / kChunkSize : 0, 0) {
    for (std::size_t k = 0; k < in.size; ++k)
      check_photon_ids(in.photon_id[k], 1);
  }

  [[nodiscard]] std::size_t size() const { return in_.size; }

  // Trace the photons [begin, end); the range must not cross a kChunkSize
  // boundary unless it starts at a chunk boundary (as parallel_for gives).
  void run(std::size_t begin, std::size_t end) {
    if (!detected_only_) {
      for (std::size_t i = begin; i < end; ++i)
        write_row(i, i, trace_one(i));
      return;
    }
    // Each chunk writes its detected rows, in input order, from the chunk's
    // first row on and records their number.
    for (std::size_t first = begin; first < end; first += kChunkSize) {
      std::size_t row = first;
      for (std::size_t i = first; i < std::min(end, first + kChunkSize); ++i) {
        const auto r = trace_one(i);
        if (r.status == Status::detected)
          write_row(row++, i, r);
      }
      found_[first / kChunkSize] = row - first;
    }
  }

  // Number of rows written; with detected_only the chunks are moved down, the
  // columns in parallel through `run` (see compact_chunks).
  template <class Run> std::size_t finish(const Run &run) const {
    if (!detected_only_)
      return in_.size;
    std::vector<Column> columns;
    add_column(columns, out_.status, 1);
    add_column(columns, out_.terminal_surface, 1);
    add_column(columns, out_.primary_surface, 1);
    add_column(columns, out_.secondary_surface, 1);
    add_column(columns, out_.position, 3);
    add_column(columns, out_.direction, 3);
    add_column(columns, out_.path_length_m, 1);
    add_column(columns, out_.time_ns, 1);
    add_column(columns, out_.weight, 1);
    add_column(columns, out_.primary_point, 3);
    add_column(columns, out_.secondary_point, 3);
    add_column(columns, out_.primary_incidence_deg, 1);
    add_column(columns, out_.secondary_incidence_deg, 1);
    add_column(columns, out_.focal_incidence_deg, 1);
    add_column(columns, out_.index, 1);
    add_column(columns, out_.reflections, 1);
    std::vector<std::size_t> chunk_start(found_.size());
    for (std::size_t c = 0; c < found_.size(); ++c)
      chunk_start[c] = c * kChunkSize;
    return compact_chunks(columns, chunk_start, found_, run);
  }

private:
  [[nodiscard]] TraceResult trace_one(std::size_t i) const {
    const Ray ray{{in_.position[3 * i], in_.position[3 * i + 1], in_.position[3 * i + 2]},
                  {in_.direction[3 * i], in_.direction[3 * i + 1], in_.direction[3 * i + 2]}};
    return trace_photon(model_, ray, in_.wavelength_nm[i], in_.time_ns[i], in_.weight[i],
                        in_.photon_id[i], seed_,
                        {.apply_responses = apply_responses_,
                         .direct = direct_,
                         .incidence_angles = wants_angles_,
                         .classify_losses = classify_losses_,
                         .from_below = ray.direction.z > 0.0});
  }

  void write_row(std::size_t row, std::size_t input, const TraceResult &r) const {
    const auto put = [row](double *target, const Vec3 &v) {
      if (!target)
        return;
      target[3 * row] = v.x;
      target[3 * row + 1] = v.y;
      target[3 * row + 2] = v.z;
    };
    const auto set = [row](auto *target, auto value) {
      if (target)
        target[row] = value;
    };
    set(out_.status, static_cast<std::uint8_t>(r.status));
    set(out_.terminal_surface, r.terminal_surface);
    set(out_.primary_surface, r.primary_surface);
    set(out_.secondary_surface, r.secondary_surface);
    put(out_.position, r.position);
    put(out_.direction, r.direction);
    put(out_.primary_point, r.primary_point);
    put(out_.secondary_point, r.secondary_point);
    set(out_.path_length_m, r.path_length_m);
    set(out_.time_ns, r.time_ns);
    set(out_.weight, r.weight);
    set(out_.primary_incidence_deg, r.primary_incidence_deg);
    set(out_.secondary_incidence_deg, r.secondary_incidence_deg);
    set(out_.focal_incidence_deg, r.focal_incidence_deg);
    set(out_.index, static_cast<std::uint64_t>(input));
    set(out_.reflections, static_cast<std::uint8_t>(r.reflections));
  }

  const OpticalModel &model_;
  PhotonBlock in_;
  ResultBlock out_;
  std::uint64_t seed_;
  bool apply_responses_, direct_, detected_only_, classify_losses_, wants_angles_;
  std::vector<std::size_t> found_; // Detected rows per chunk (detected_only).
};

} // namespace detail

// Trace a block on a persistent pool. Every photon writes only its own row,
// and its random numbers depend only on (key, photon_id), so results do not
// depend on the thread count or the order of photons. Returns the number of
// rows written: in.size, or with `detected_only` the number of detected
// photons, whose rows come first, in input order (the rest of the output is
// left unwritten).
inline std::size_t trace_block(const OpticalModel &model, const PhotonBlock &in,
                               const ResultBlock &out, std::uint64_t seed, bool apply_responses,
                               ThreadPool &pool, bool direct = false, bool detected_only = false,
                               bool classify_losses = true) {
  detail::BlockTask task(model, in, out, seed, apply_responses, direct, detected_only,
                         classify_losses);
  parallel_for(in.size, pool, [&](std::size_t begin, std::size_t end) { task.run(begin, end); });
  return task.finish([&](std::size_t count, const auto &each) { pool.run(count, each); });
}

// The same on `threads` threads: a pool for this call (see pool_size). The results are equal
// to those of the pool version.
inline std::size_t trace_block(const OpticalModel &model, const PhotonBlock &in,
                               const ResultBlock &out, std::uint64_t seed, bool apply_responses,
                               unsigned threads, bool direct = false, bool detected_only = false,
                               bool classify_losses = true) {
  ThreadPool pool(pool_size(in.size, threads));
  return trace_block(model, in, out, seed, apply_responses, pool, direct, detected_only,
                     classify_losses);
}

// One block of trace_blocks.
struct BlockRequest {
  const OpticalModel *model{};
  PhotonBlock in;
  ResultBlock out;
  bool apply_responses{true}, direct{false}, detected_only{false}, classify_losses{true};
};

// Trace several blocks, each with its own model, in one pool job: all chunks of
// all blocks form one flat list. Returns the rows written per block; the
// results equal one trace_block call per block.
inline std::vector<std::size_t> trace_blocks(const std::vector<BlockRequest> &requests,
                                             std::uint64_t seed, ThreadPool &pool) {
  std::vector<detail::BlockTask> tasks;
  tasks.reserve(requests.size());
  std::vector<std::size_t> first_chunk(requests.size() + 1, 0); // Chunks before each block.
  for (std::size_t b = 0; b < requests.size(); ++b) {
    const BlockRequest &r = requests[b];
    tasks.emplace_back(*r.model, r.in, r.out, seed, r.apply_responses, r.direct, r.detected_only,
                       r.classify_losses);
    first_chunk[b + 1] = first_chunk[b] + (r.in.size + kChunkSize - 1) / kChunkSize;
  }
  pool.run(first_chunk.back(), [&](std::size_t c) {
    const std::size_t b = static_cast<std::size_t>(
        std::upper_bound(first_chunk.begin(), first_chunk.end(), c) - first_chunk.begin() - 1);
    const std::size_t begin = (c - first_chunk[b]) * kChunkSize;
    tasks[b].run(begin, std::min(tasks[b].size(), begin + kChunkSize));
  });
  std::vector<std::size_t> rows;
  rows.reserve(tasks.size());
  for (const auto &task : tasks)
    rows.push_back(
        task.finish([&](std::size_t count, const auto &each) { pool.run(count, each); }));
  return rows;
}

} // namespace obdeect::ray_tracing
