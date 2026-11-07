// dlamb-ls-test: figure-eight light source with 4e7 photons per position, as a gstrickd photon
// stream (port of LightEmission ls-test). Writes to stdout, or to -o FILE.

#include "dlamb/programs/ls_test.hpp"
#include "dlamb/cli.hpp"

#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void usage(std::FILE *to) {
  std::fprintf(to, "usage: dlamb-ls-test [options]\n"
                   "Writes a gstrickd photon stream (sky frame) to stdout, or to -o FILE.\n"
                   "  --photons N   photons into 4 pi per position [4e7]\n"
                   "  --seed S      random seed [1]\n"
                   "  -o FILE       output file (default stdout)\n"
                   "  --help        this text\n");
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  dlamb::ls_test::Options options;
  std::string output;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--help") {
        usage(stdout);
        return 0;
      }
      if (i + 1 >= argc)
        throw std::invalid_argument("option " + arg + " needs a value");
      const std::string value = argv[++i];
      if (arg == "--photons")
        options.photons_per_position = dlamb::cli::parse_number(arg, value);
      else if (arg == "--seed")
        options.seed = static_cast<std::uint64_t>(dlamb::cli::parse_integer(arg, value));
      else if (arg == "-o" || arg == "--output")
        output = value;
      else
        throw std::invalid_argument("unknown option " + arg);
    }
    dlamb::StreamSummary summary;
    if (output.empty()) {
      summary = dlamb::ls_test::write(std::cout, options);
      std::cout.flush();
    } else {
      std::ofstream file(output, std::ios::binary);
      if (!file)
        throw std::runtime_error("cannot open " + output);
      summary = dlamb::ls_test::write(file, options);
    }
    std::fprintf(stderr, "dlamb-ls-test: expected weight %.6g, wrote %llu photons, weight %.6g\n",
                 summary.expected, static_cast<unsigned long long>(summary.photons),
                 summary.weight);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "dlamb-ls-test: %s\n", error.what());
    usage(stderr);
    return 1;
  }
  return 0;
}
