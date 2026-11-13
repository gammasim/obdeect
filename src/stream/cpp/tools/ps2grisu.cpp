// obdeect stream to-grisu: ground-frame photon stream (file or stdin) to GrISU/GrOptics photon text on stdout or
// -o FILE.

#include "obdeect/stream/grisu.hpp"
#include "obdeect/stream/stdin_buffer.hpp"
#include "tool_version.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr std::size_t buffer_size = std::size_t{1} << 20;

void usage() {
  std::fprintf(stderr, "usage: obdeect stream to-grisu [--telescope-config FILE] [--qe FACTOR] "
                       "[--expand-bunches] [--event-details] [-o FILE|-] [input|-]\n");
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  std::string input = "-";
  std::string output = "-";
  std::string config_path;
  obdeect::stream::GrisuOptions options;
  options.producer = std::string("obdeect stream to-grisu ") + OBDEECT_STREAM_PACKAGE_VERSION;
  options.warnings = &std::cerr;
  bool have_input = false;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--version") {
      print_version("obdeect stream to-grisu");
      return 0;
    }
    if (arg == "-h" || arg == "--help") {
      usage();
      return 0;
    }
    if (arg == "--expand-bunches") {
      options.expand_bunches = true;
      continue;
    }
    if (arg == "--event-details") {
      options.event_details = true;
      continue;
    }
    if (arg == "--telescope-config" || arg == "--qe" || arg == "-o") {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "obdeect stream to-grisu: %s needs an argument\n", arg.c_str());
        usage();
        return 2;
      }
      const std::string value = argv[++i];
      if (arg == "-o") {
        output = value;
      } else if (arg == "--telescope-config") {
        config_path = value;
      } else {
        char *end = nullptr;
        options.qe = std::strtod(value.c_str(), &end);
        if (value.empty() || end != value.c_str() + value.size()) {
          std::fprintf(stderr, "obdeect stream to-grisu: invalid number for --qe: %s\n", value.c_str());
          return 2;
        }
      }
      continue;
    }
    if (arg.size() > 1 && arg[0] == '-') {
      std::fprintf(stderr, "obdeect stream to-grisu: unknown option %s\n", arg.c_str());
      usage();
      return 2;
    }
    if (have_input) {
      std::fprintf(stderr, "obdeect stream to-grisu: only one input is accepted\n");
      usage();
      return 2;
    }
    input = arg;
    have_input = true;
  }

  try {
    if (!config_path.empty()) {
      std::ifstream config(config_path);
      if (!config) {
        throw std::runtime_error("cannot open telescope config " + config_path);
      }
      options.telescopes = obdeect::stream::read_tlloc(config);
    }

    static std::vector<char> in_buffer(buffer_size);
    static std::vector<char> out_buffer(buffer_size);
    obdeect::stream::StdinBuffer stdin_buffer;
    std::istream stdin_stream(&stdin_buffer);
    std::ifstream in_file;
    std::istream *in = &stdin_stream;
    if (input != "-") {
      in_file.rdbuf()->pubsetbuf(in_buffer.data(), static_cast<std::streamsize>(in_buffer.size()));
      in_file.open(input, std::ios::binary);
      if (!in_file) {
        throw std::runtime_error("cannot open " + input);
      }
      in = &in_file;
    }

    std::ofstream out_file;
    std::ostream *out = &std::cout;
    if (output == "-") {
      std::cout.rdbuf()->pubsetbuf(out_buffer.data(),
                                   static_cast<std::streamsize>(out_buffer.size()));
    } else {
      out_file.rdbuf()->pubsetbuf(out_buffer.data(),
                                  static_cast<std::streamsize>(out_buffer.size()));
      out_file.open(output, std::ios::binary);
      if (!out_file) {
        throw std::runtime_error("cannot open " + output + " for writing");
      }
      out = &out_file;
    }

    obdeect::stream::stream_to_grisu(*in, *out, std::move(options));
    out->flush();
    if (!*out) {
      throw std::runtime_error("write error on output");
    }
    if (out_file.is_open()) {
      out_file.close();
    }
  } catch (const std::exception &error) {
    std::fprintf(stderr, "obdeect stream to-grisu: %s\n", error.what());
    return 1;
  }
  return 0;
}
