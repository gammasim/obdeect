// Reference copy of the sequential trace of downward photons as it was before the
// non-sequential loop (ray tracing phase 1, downward branch). Used only by
// test_trace to prove that the loop gives bit-identical results for sequential photons.
#pragma once

#include "obdeect/ray_tracing/trace.hpp"

#include <limits>
#include <optional>

namespace obdeect::ray_tracing::legacy {

inline TraceResult trace_sequential(const OpticalModel &model, Ray ray, double wavelength_nm,
                                    double time_ns, double weight, std::uint64_t photon_id,
                                    const RandomKey &seed, bool apply_responses = true,
                                    bool direct = false, bool incidence_angles = true,
                                    bool classify_losses = true, std::uint32_t sub_photon = 0) {
  using namespace detail;
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
  const auto finish = [&](Status status, const Vec3 &position, std::uint32_t surface) {
    out.status = status;
    out.position = position;
    out.direction = ray.direction;
    out.terminal_surface = surface;
    out.path_length_m = path;
    out.time_ns = time_ns + path * model.propagation_group_index / kSpeedOfLightMPerNs;
    if (status != Status::detected)
      out.weight = 0.0;
    return out;
  };
  const auto blocked_on = [&](Leg leg, double t_limit) -> std::optional<TraceResult> {
    if (const auto b = model.obscurers.nearest(ray, leg, t_limit)) {
      path += b->t;
      return finish(Status::blocked, ray.at(b->t), b->id);
    }
    return std::nullopt;
  };
  bool failed = false;
  const auto continue_from_primary = [&](const MirrorHit &m1) -> std::optional<TraceResult> {
    path += m1.hit.t;
    out.primary_point = m1.hit.point;
    out.primary_surface = m1.hit.id;
    if (incidence_angles || angle_matters(model.primary))
      out.primary_incidence_deg = incidence_deg(ray.direction, m1.hit.normal);
    if (apply_responses)
      out.weight *=
          mirror_response(model.primary, m1, wavelength_nm,
                          std::isnan(out.primary_incidence_deg) ? 0.0 : out.primary_incidence_deg);
    ray = {m1.hit.point, reflect_at(model.primary, m1, ray.direction, photon_id, seed,
                                    Stream::primary_scatter, sub_photon)};
    if (model.secondary) {
      const auto m2 = intersect_mirror(*model.secondary, ray, failed);
      if (failed)
        return finish(Status::numerical_failure, ray.origin, kNoSurface);
      if (auto b = blocked_on(Leg::after_primary,
                              m2 ? m2->hit.t : std::numeric_limits<double>::infinity()))
        return b;
      if (!m2)
        return finish(Status::missed_secondary, ray.origin, kNoSurface);
      path += m2->hit.t;
      out.secondary_point = m2->hit.point;
      out.secondary_surface = m2->hit.id;
      if (incidence_angles || angle_matters(*model.secondary))
        out.secondary_incidence_deg = incidence_deg(ray.direction, m2->hit.normal);
      if (apply_responses)
        out.weight *= mirror_response(
            *model.secondary, *m2, wavelength_nm,
            std::isnan(out.secondary_incidence_deg) ? 0.0 : out.secondary_incidence_deg);
      ray = {m2->hit.point, reflect_at(*model.secondary, *m2, ray.direction, photon_id, seed,
                                       Stream::secondary_scatter, sub_photon)};
    }
    return std::nullopt;
  };
  const auto focal_hit = [&]() -> std::optional<SurfaceHit> {
    const auto &fs = model.focal_surface;
    std::optional<SurfaceHit> focal;
    if (fs.asphere) {
      focal = fs.asphere->intersect(ray, kMinDistanceM, [](const Vec3 &) { return true; }, failed);
    } else if (ray.direction.z != 0.0) {
      const double t = (fs.plane_z - ray.origin.z) / ray.direction.z;
      if (t > kMinDistanceM)
        focal = SurfaceHit{t, ray.at(t), Vec3{0.0, 0.0, 1.0}, fs.id};
    }
    return focal;
  };

  if (!direct) {
    if (apply_responses && model.telescope_transmission)
      out.weight *= model.telescope_transmission->at(ray.direction);
    const auto m1 = intersect_mirror(model.primary, ray, failed);
    if (failed)
      return finish(Status::numerical_failure, ray.origin, kNoSurface);
    double t_m1 = m1 ? m1->hit.t : std::numeric_limits<double>::infinity();
    if (!classify_losses && !m1)
      return finish(Status::missed_primary, ray.origin, kNoSurface);
    if (model.secondary_blocks_incoming) {
      const auto back = intersect_mirror(*model.secondary, ray, failed);
      if (failed)
        return finish(Status::numerical_failure, ray.origin, kNoSurface);
      if (back && back->hit.t < t_m1) {
        if (auto b = blocked_on(Leg::incoming, back->hit.t))
          return *b;
        path += back->hit.t;
        return finish(Status::blocked, back->hit.point, back->hit.id);
      }
    }
    if (auto b = blocked_on(Leg::incoming, t_m1))
      return *b;
    if (!m1)
      return finish(Status::missed_primary, ray.origin, kNoSurface);
    if (auto r = continue_from_primary(*m1))
      return *r;
  }
  const Leg leg = model.secondary ? Leg::after_secondary : Leg::after_primary;
  const auto &fs = model.focal_surface;
  const auto focal = focal_hit();
  if (failed)
    return finish(Status::numerical_failure, ray.origin, kNoSurface);
  if (auto b = blocked_on(leg, focal ? focal->t : std::numeric_limits<double>::infinity()))
    return *b;
  if (!focal)
    return finish(Status::missed_focal_surface, ray.origin, kNoSurface);
  path += focal->t;
  if (incidence_angles)
    out.focal_incidence_deg = incidence_deg(ray.direction, focal->normal);
  return finish(Status::detected, focal->point, fs.id);
}

} // namespace obdeect::ray_tracing::legacy
