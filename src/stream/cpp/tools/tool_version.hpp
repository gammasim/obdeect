#pragma once

// --version output shared by the command-line tools.

#include "obdeect/stream/format.hpp"

#include <cstdio>

#ifndef OBDEECT_STREAM_PACKAGE_VERSION
#error "OBDEECT_STREAM_PACKAGE_VERSION is defined by CMake for each tool"
#endif

// Prints "<tool> <package version> (photon-stream format <major>.<minor>)".
inline void print_version(const char *tool) {
  std::printf("%s %s (photon-stream format %d.%d)\n", tool, OBDEECT_STREAM_PACKAGE_VERSION,
              obdeect::stream::format_major, obdeect::stream::format_minor);
}
