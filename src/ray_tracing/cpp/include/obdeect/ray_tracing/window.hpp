#pragma once

#include "obdeect/ray_tracing/model.hpp"

#include <cmath>
#include <optional>

namespace obdeect::ray_tracing {

// Intersect either side of a window surface. Unlike a mirror facet this is a
// transmitting surface, so both signs of the normal are valid.
[[nodiscard]] inline std::optional<SurfaceHit>
intersect_window_surface(const Window &window, const Ray &ray, bool inner, double t_min) {
  if (window.shape == Window::Shape::flat) {
    const double z = window.outer_apex_z + (inner ? window.thickness : 0.0);
    if (ray.direction.z == 0.0)
      return std::nullopt;
    const double t = (z - ray.origin.z) / ray.direction.z;
    if (!(t > t_min))
      return std::nullopt;
    return SurfaceHit{t, ray.at(t), {0.0, 0.0, 1.0}, window.id};
  }
  const double radius = window.radius_of_curvature - (inner ? window.thickness : 0.0);
  const Vec3 centre{0.0, 0.0, window.outer_apex_z + window.radius_of_curvature};
  const Vec3 offset = ray.origin - centre;
  const double b = dot(offset, ray.direction);
  const double discriminant = b * b - (dot(offset, offset) - radius * radius);
  if (discriminant < 0.0)
    return std::nullopt;
  const double root = std::sqrt(discriminant);
  for (const double t : {-b - root, -b + root}) {
    if (t > t_min) {
      const Vec3 point = ray.at(t);
      return SurfaceHit{t, point, normalized(point - centre), window.id};
    }
  }
  return std::nullopt;
}

} // namespace obdeect::ray_tracing
