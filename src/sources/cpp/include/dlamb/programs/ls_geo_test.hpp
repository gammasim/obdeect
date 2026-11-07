#pragma once

// Port of LightEmission ls-geo-test.cc: geometry checks (points, lines, planes, rotations)
// with the same test points as the original. It emits no photons, so it writes a report
// instead of a photon stream. Lengths are in metres here (the original uses cm; the
// numbers are the same, only the unit label differs).

#include "dlamb/density.hpp"
#include "dlamb/vec3.hpp"

#include <cmath>
#include <optional>
#include <ostream>

namespace dlamb::ls_geo_test {

// Line through `point` with unit direction `direction`.
struct Line {
  Vec3 point{};
  Vec3 direction{0.0, 0.0, 1.0};
};

// Plane through `point` with unit normal `normal`.
struct Plane {
  Vec3 point{};
  Vec3 normal{0.0, 0.0, 1.0};
};

[[nodiscard]] inline Line line(Vec3 point, Vec3 direction) {
  return {point, normalised(direction)};
}

[[nodiscard]] inline Plane plane(Vec3 point, Vec3 normal) { return {point, normalised(normal)}; }

// Distance from a point to a line.
[[nodiscard]] inline double distance(Vec3 p, const Line &l) {
  return norm(cross(p - l.point, l.direction));
}

// Distance between two lines (skew: common perpendicular; parallel: point to line).
[[nodiscard]] inline double distance(const Line &a, const Line &b) {
  const Vec3 n = cross(a.direction, b.direction);
  if (norm(n) < 1e-12)
    return distance(b.point, a);
  return std::abs(dot(b.point - a.point, n)) / norm(n);
}

// Distance from a point to a plane.
[[nodiscard]] inline double distance(Vec3 p, const Plane &pl) {
  return std::abs(dot(p - pl.point, pl.normal));
}

// Intersection of a line and a plane; empty if the line is parallel to the plane.
[[nodiscard]] inline std::optional<Vec3> intersection(const Line &l, const Plane &pl) {
  const double denominator = dot(l.direction, pl.normal);
  if (std::abs(denominator) < 1e-12)
    return std::nullopt;
  const double t = dot(pl.point - l.point, pl.normal) / denominator;
  return l.point + l.direction * t;
}

// Right-handed rotations by `angle` [rad] about the coordinate axes.
[[nodiscard]] inline Vec3 rotate_x(Vec3 v, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return {v.x, c * v.y - s * v.z, s * v.y + c * v.z};
}
[[nodiscard]] inline Vec3 rotate_y(Vec3 v, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return {c * v.x + s * v.z, v.y, -s * v.x + c * v.z};
}
[[nodiscard]] inline Vec3 rotate_z(Vec3 v, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return {c * v.x - s * v.y, s * v.x + c * v.y, v.z};
}

// Writes the report of ls-geo-test.cc (same points, same quantities).
inline void write_report(std::ostream &os) {
  const Vec3 o{0.0, 0.0, 0.0}, px{10.0, 0.0, 0.0}, py{0.0, 10.0, 0.0}, pz{0.0, 0.0, 10.0};
  const Vec3 p345{3.0, 4.0, 5.0};
  const Line lx = line(o, {1, 0, 0}), ly = line(o, {0, 1, 0}), lz = line(o, {0, 0, 1});
  const Line l1 = line(px, {0, 0, 1}), l2 = line(py, {1, 0, 0}), l3 = line(pz, {0, 1, 0});
  const Plane a = plane(p345, {1.0, 2.0, 0.0});
  os << "O  = (0, 0, 0)\n";
  os << "Distance O - Px = " << norm(px - o) << " m\n";
  os << "Distance O - Py = " << norm(py - o) << " m\n";
  os << "Distance O - Pz = " << norm(pz - o) << " m\n";
  os << "Distance y axis - L1 = " << distance(l1, ly) << " m\n";
  os << "Distance z axis - L2 = " << distance(l2, lz) << " m\n";
  os << "Distance x axis - L3 = " << distance(l3, lx) << " m\n";
  os << "Distance x axis - P345 = " << distance(p345, lx) << " m\n";
  os << "Distance y axis - P345 = " << distance(p345, ly) << " m\n";
  os << "Distance z axis - P345 = " << distance(p345, lz) << " m\n";
  os << "Length of P345 vector = " << norm(p345) << " m\n";
  const double angle_x = std::acos(dot(normalised(p345), lx.direction)) * 180.0 / pi;
  os << "Angle of P345 with x axis = " << angle_x << " deg\n";
  os << "P345 rotated 90 deg around x axis = (" << rotate_x(p345, pi / 2).x << ", "
     << rotate_x(p345, pi / 2).y << ", " << rotate_x(p345, pi / 2).z << ")\n";
  os << "Distance x/y plane - P345 = " << distance(p345, plane(o, {0, 0, 1})) << " m\n";
  os << "Distance O - A = " << distance(o, a) << " m\n";
  os << "Distance P345 - A = " << distance(p345, a) << " m\n";
  if (const auto p = intersection(lx, a))
    os << "A intersects x axis at (" << p->x << ", " << p->y << ", " << p->z << ")\n";
  else
    os << "A has no intersection with x axis\n";
  if (const auto p = intersection(ly, a))
    os << "A intersects y axis at (" << p->x << ", " << p->y << ", " << p->z << ")\n";
  else
    os << "A has no intersection with y axis\n";
  if (const auto p = intersection(lz, a))
    os << "A intersects z axis at (" << p->x << ", " << p->y << ", " << p->z << ")\n";
  else
    os << "A has no intersection with z axis\n";
}

} // namespace dlamb::ls_geo_test
