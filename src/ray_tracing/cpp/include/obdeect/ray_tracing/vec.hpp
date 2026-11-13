#pragma once

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/vec3.hpp"

#include <cmath>

namespace obdeect::ray_tracing {

// Speed of light in vacuum, metres per nanosecond.
inline constexpr double kSpeedOfLightMPerNs = obdeect::stream::speed_of_light_m_per_ns;
// Smallest accepted propagation distance; rejects self-intersections at the
// point a ray left.
constexpr double kMinDistanceM = 1e-9;

// The vector type is shared with the other components.
using obdeect::stream::cross;
using obdeect::stream::dot;
using obdeect::stream::norm;
using obdeect::stream::Vec3;
inline Vec3 normalized(const Vec3 &a) { return obdeect::stream::normalised(a); }
inline bool finite(const Vec3 &a) { return obdeect::stream::is_finite(a); }

// Ray with a unit direction.
struct Ray {
  Vec3 origin;
  Vec3 direction;
  [[nodiscard]] Vec3 at(double t) const { return origin + direction * t; }
};

// Specular reflection of a unit direction about a unit normal.
inline Vec3 reflect(const Vec3 &d, const Vec3 &n) { return d - n * (2.0 * dot(d, n)); }

} // namespace obdeect::ray_tracing
