#pragma once

// --version output shared by the command-line tools.

#include "gstrickd/format.hpp"

#include <cstdio>

#ifndef GSTRICKD_PACKAGE_VERSION
#error "GSTRICKD_PACKAGE_VERSION is defined by CMake for each tool"
#endif

// Prints "<tool> <package version> (photon-stream format <major>.<minor>)".
inline void print_version(const char *tool) {
  std::printf("%s %s (photon-stream format %d.%d)\n", tool, GSTRICKD_PACKAGE_VERSION,
              gstrickd::format_major, gstrickd::format_minor);
}
