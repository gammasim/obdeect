#pragma once

// Plumbing shared by the tool executables: argument cursor, output file, seed text, and the
// top-level error handling (see stream/cli.hpp: usage error 2, other failure 1).

#include "obdeect/sources/cli.hpp"

#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::sources::cli {

// A seed: negative values convert to uint64 modulo 2^64.
inline std::uint64_t parse_seed(const std::string &text, const std::string &option) {
  const long long value = parse_integer(text, option);
  if (value < 0)
    return static_cast<std::uint64_t>(value);
  return std::strtoull(text.c_str(), nullptr, 10);
}

using stream::cli::Args;

// Output to a file (opened, truncated, at construction) or to stdout.
class Output {
public:
  explicit Output(const std::string &path) {
    if (path.empty() || path == "-")
      return;
    file_ = std::make_unique<std::ofstream>(path, std::ios::binary | std::ios::trunc);
    if (!*file_)
      throw std::runtime_error("cannot open output file " + path);
  }
  std::ostream &stream() { return file_ ? static_cast<std::ostream &>(*file_) : std::cout; }
  void close() {
    if (file_)
      file_->close();
    else
      std::cout.flush();
  }

private:
  std::unique_ptr<std::ofstream> file_;
};

// Runs a tool body (stream::cli::run_tool). Option validation in the sources core throws
// std::invalid_argument; here it is a usage error (exit status 2).
template <class Body>
int run_tool(const std::string &prefix, int argc, char **argv, Body body,
             const std::function<void()> &usage = {}) {
  return stream::cli::run_tool(prefix, stream::cli::core_version(), argc, argv, usage, [&] {
    try {
      return body();
    } catch (const UsageError &) {
      throw;
    } catch (const std::invalid_argument &error) {
      throw UsageError(error.what());
    }
  });
}

} // namespace obdeect::sources::cli
