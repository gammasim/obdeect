#pragma once

// Three-component vectors and the orthonormal basis used for azimuths, shared by sources and
// ray_tracing. Templated on the scalar type so that kernels built on them can later run on other
// scalar types.

#include <cmath>

namespace obdeect::stream {

template <class T> struct BasicVec3 {
  T x{}, y{}, z{};
};
using Vec3 = BasicVec3<double>;

template <class T> constexpr BasicVec3<T> operator+(BasicVec3<T> a, BasicVec3<T> b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <class T> constexpr BasicVec3<T> operator-(BasicVec3<T> a, BasicVec3<T> b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
template <class T> constexpr BasicVec3<T> operator-(BasicVec3<T> a) { return {-a.x, -a.y, -a.z}; }
template <class T> constexpr BasicVec3<T> operator*(BasicVec3<T> a, T s) {
  return {a.x * s, a.y * s, a.z * s};
}
template <class T> constexpr BasicVec3<T> operator*(T s, BasicVec3<T> a) { return a * s; }

template <class T> constexpr T dot(BasicVec3<T> a, BasicVec3<T> b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <class T> constexpr BasicVec3<T> cross(BasicVec3<T> a, BasicVec3<T> b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <class T> T norm(BasicVec3<T> a) { return std::sqrt(dot(a, a)); }
template <class T> BasicVec3<T> normalised(BasicVec3<T> a) { return a * (T{1} / norm(a)); }
// Non-template overload so that braced vectors work: normalised({0, 0, 1}).
inline Vec3 normalised(Vec3 a) { return normalised<double>(a); }
template <class T> bool is_finite(BasicVec3<T> a) {
  return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}

// Two unit vectors that complete a unit vector n to a right-handed basis (e1, e2, n).
// The choice fixes the zero point of azimuths measured around n.
template <class T> struct BasicBasis {
  BasicVec3<T> e1, e2;
};
using Basis = BasicBasis<double>;

template <class T> BasicBasis<T> orthonormal_basis(BasicVec3<T> n) {
  const BasicVec3<T> helper =
      std::abs(n.x) < T{0.9} ? BasicVec3<T>{T{1}, T{0}, T{0}} : BasicVec3<T>{T{0}, T{1}, T{0}};
  const BasicVec3<T> e1 = normalised(cross(n, helper));
  return {e1, cross(n, e1)};
}

// Direction at polar angle theta (cosine `cos_theta`, sine `sin_theta`) and azimuth phi around
// the unit vector n, with the azimuth zero point of orthonormal_basis(n).
template <class T>
BasicVec3<T> direction_about(BasicVec3<T> n, const BasicBasis<T> &basis, T cos_theta, T sin_theta,
                             T phi) {
  return n * cos_theta + (basis.e1 * std::cos(phi) + basis.e2 * std::sin(phi)) * sin_theta;
}

} // namespace obdeect::stream
