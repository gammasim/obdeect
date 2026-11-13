#pragma once

// run_tool for the stream tools: the version line carries the CMake package version.

#include "obdeect/stream/cli.hpp"

#ifndef OBDEECT_STREAM_PACKAGE_VERSION
#error "OBDEECT_STREAM_PACKAGE_VERSION is defined by CMake for each tool"
#endif

template <class Body>
int run_stream_tool(const char *tool, int argc, char **argv, const std::function<void()> &usage,
                    Body body) {
  return obdeect::stream::cli::run_tool(tool, OBDEECT_STREAM_PACKAGE_VERSION, argc, argv, usage,
                                        body);
}
