#pragma once

// Stream-file conversion of the Python API: open the input, write the output to
// `<output>.partial` and rename it to `<output>` only when the conversion succeeds, so a failed
// run leaves no partial output behind. Used by the propagate and to_sky bindings.

#include "obdeect/propagation/propagator.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>

namespace obdeect_propagation_bindings {

// Runs convert(in, out) and returns its result. The output exists only after a successful run.
template <class Convert>
auto convert_stream_file(const std::string &input_path, const std::string &output_path,
                         Convert &&convert) {
  if (obdeect::propagation::same_file(input_path, output_path)) {
    throw std::invalid_argument("input and output are the same file: " + input_path);
  }
  std::ifstream in(input_path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("cannot open input stream: " + input_path);
  }
  const std::string partial = output_path + ".partial";
  const auto remove_partial = [&partial] {
    std::error_code ignored;
    std::filesystem::remove(partial, ignored);
  };
  try {
    std::ofstream out(partial, std::ios::binary);
    if (!out) {
      throw std::runtime_error("cannot open output stream: " + output_path);
    }
    auto result = convert(in, out);
    out.close();
    if (!out) {
      throw std::runtime_error("cannot write output stream: " + output_path);
    }
    std::filesystem::rename(partial, output_path);
    return result;
  } catch (...) {
    remove_partial();
    throw;
  }
}

} // namespace obdeect_propagation_bindings
