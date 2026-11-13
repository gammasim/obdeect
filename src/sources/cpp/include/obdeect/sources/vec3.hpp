#pragma once

// The vector type lives in obdeect/stream/vec3.hpp (shared with ray_tracing); this header only
// brings it into obdeect::sources together with the constants the sources use.

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/vec3.hpp"

namespace obdeect::sources {

inline constexpr double pi = obdeect::stream::pi;
inline constexpr double two_pi = obdeect::stream::two_pi;
inline constexpr double speed_of_light_m_per_ns = obdeect::stream::speed_of_light_m_per_ns;

using obdeect::stream::BasicBasis;
using obdeect::stream::BasicVec3;
using obdeect::stream::Basis;
using obdeect::stream::cross;
using obdeect::stream::direction_about;
using obdeect::stream::dot;
using obdeect::stream::is_finite;
using obdeect::stream::norm;
using obdeect::stream::normalised;
using obdeect::stream::orthonormal_basis;
using obdeect::stream::Vec3;

} // namespace obdeect::sources
