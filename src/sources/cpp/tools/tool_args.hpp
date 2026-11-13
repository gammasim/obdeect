#pragma once

// Plumbing shared by the tool executables: argument cursor, output file, seed text, and the
// top-level error handling (one line to stderr, exit 1).

#include "obdeect/sources/cli.hpp"

#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
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

// Walks argv; value() takes the argument after an option or fails with `<arg>: missing value`.
class ArgCursor {
public:
  ArgCursor(int argc, char **argv) {
    for (int i = 1; i < argc; ++i)
      args_.emplace_back(argv[i]);
  }
  bool next(std::string &arg) {
    if (i_ >= args_.size())
      return false;
    arg = args_[i_++];
    return true;
  }
  std::string value(const std::string &arg) {
    if (i_ >= args_.size())
      throw std::invalid_argument(arg + ": missing value");
    return args_[i_++];
  }

private:
  std::vector<std::string> args_;
  std::size_t i_{0};
};

// Output to a file (opened, truncated, at construction) or to stdout.
class Output {
public:
  explicit Output(const std::string &path) {
    if (path.empty() || path == "-")
      return;
    file_ = std::make_unique<std::ofstream>(path, std::ios::binary | std::ios::trunc);
    if (!*file_)
      throw std::invalid_argument("cannot open output file " + path);
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

// Runs a tool body; an exception prints `<prefix>: <message>` and gives exit status 1.
template <class Body> int run_tool(const std::string &prefix, Body body) {
  try {
    return body();
  } catch (const std::exception &error) {
    std::cerr << prefix << ": " << error.what() << '\n';
    return 1;
  }
}

} // namespace obdeect::sources::cli
