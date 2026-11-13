#pragma once

// Command-line scaffolding shared by every native tool: argument cursor, option parsing built on
// parse.hpp, `--version`, and one exit-status rule.
//
//   0  success (also --help and --version)
//   1  runtime failure (unreadable input, bad data, write error)
//   2  usage error (UsageError)
//
// Messages go to stderr as `<tool>: <message>`; a usage error is followed by the tool's usage text.

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/parse.hpp"
#include "obdeect/stream/version.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::stream::cli {

// A bad command line. Library errors (std::invalid_argument and others) are runtime failures.
struct UsageError : std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};

// Walks argv; value() takes the argument after an option or fails with `<arg>: missing value`.
class Args {
public:
  Args(int argc, char **argv) {
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
      throw UsageError(arg + ": missing value");
    return args_[i_++];
  }

private:
  std::vector<std::string> args_;
  std::size_t i_{0};
};

// A finite real number; the whole text must be a number.
inline double number(const std::string &text, const std::string &option) {
  double value = 0.0;
  if (!parse::finite_number(text, value))
    throw UsageError(option + ": not a finite number: " + text);
  return value;
}

// A decimal integer filling the whole text and lying in lo..hi.
inline long long integer(const std::string &text, const std::string &option,
                         long long lo = std::numeric_limits<long long>::min(),
                         long long hi = std::numeric_limits<long long>::max()) {
  long long value = 0;
  if (!parse::integer(text, value))
    throw UsageError(option + ": not an integer: " + text);
  if (value < lo || value > hi)
    throw UsageError(option + ": integer out of range: " + text);
  return value;
}

// The package version as "major.minor.patch" (sources and propagation have no version of their
// own).
inline std::string core_version() {
  return std::to_string(version_major) + "." + std::to_string(version_minor) + "." +
         std::to_string(version_patch);
}

// "<tool> <package version> (photon-stream format <major>.<minor>)".
inline std::string version_line(const std::string &tool, const std::string &version) {
  return tool + " " + version + " (photon-stream format " + std::to_string(format_major) + "." +
         std::to_string(format_minor) + ")\n";
}

// Runs a tool body. `--version` anywhere on the command line prints the version line and returns
// 0 before the body runs. A usage error prints `<tool>: <message>`, then `usage` (may be empty),
// and returns 2; any other exception prints `<tool>: <message>` and returns 1.
template <class Body>
int run_tool(const std::string &tool, const std::string &version, int argc, char **argv,
             const std::function<void()> &usage, Body body) {
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--version") {
      std::fputs(version_line(tool, version).c_str(), stdout);
      return 0;
    }
  }
  try {
    return body();
  } catch (const UsageError &error) {
    std::cerr << tool << ": " << error.what() << '\n';
    if (usage)
      usage();
    return 2;
  } catch (const std::exception &error) {
    std::cerr << tool << ": " << error.what() << '\n';
    return 1;
  }
}

} // namespace obdeect::stream::cli
