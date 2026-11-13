#pragma once

// Version of the C++ core: the obdeect stream package version (one version for all components).

#include "obdeect/stream/version.hpp"

namespace obdeect::propagation {

inline constexpr int version_major = obdeect::stream::version_major;
inline constexpr int version_minor = obdeect::stream::version_minor;
inline constexpr int version_patch = obdeect::stream::version_patch;

} // namespace obdeect::propagation
