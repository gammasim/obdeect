#pragma once

#include <cmath>

namespace dlamb {

// Three-vector in metres (positions) or unit-free components (directions).
struct Vec3 {
  double x{};
  double y{};
  double z{};
};

[[nodiscard]] constexpr Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
[[nodiscard]] constexpr Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
[[nodiscard]] constexpr Vec3 operator-(Vec3 a) { return {-a.x, -a.y, -a.z}; }
[[nodiscard]] constexpr Vec3 operator*(Vec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
[[nodiscard]] constexpr Vec3 operator*(double s, Vec3 a) { return a * s; }

[[nodiscard]] constexpr double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

[[nodiscard]] constexpr Vec3 cross(Vec3 a, Vec3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

[[nodiscard]] inline double norm(Vec3 a) { return std::sqrt(dot(a, a)); }

[[nodiscard]] inline Vec3 normalised(Vec3 a) { return a * (1.0 / norm(a)); }

// Two unit vectors orthogonal to the unit vector n, right-handed with n.
struct Basis {
  Vec3 e1;
  Vec3 e2;
};

[[nodiscard]] inline Basis orthonormal_basis(Vec3 n) {
  const Vec3 e1 =
      normalised(std::abs(n.x) < 0.9 ? cross(n, {1.0, 0.0, 0.0}) : cross(n, {0.0, 1.0, 0.0}));
  return {e1, cross(n, e1)};
}

[[nodiscard]] inline bool is_finite(Vec3 a) {
  return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}

} // namespace dlamb
