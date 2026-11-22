// obdeect propagation to-sky: converts the frame G photons of a photon stream (CORSIKA, after
// `obdeect stream from-corsika`) to frame sky, the inverse of the atmosphere stage's sky-to-G step
// (docs/propagation/ground-to-sky.md). Reads a stream from a file or stdin, writes the converted
// stream to a file or stdout, and prints a summary to stderr.

#include "obdeect/propagation/ground_to_sky.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/propagator.hpp"
#include "obdeect/stream/cli.hpp"
#include "obdeect/stream/stdin_buffer.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

using namespace obdeect::propagation;

namespace {

void usage() {
  std::fprintf(stderr, "usage: obdeect propagation to-sky [--profile FILE]\n"
                       "                       [-i FILE|-] [-o FILE|-]\n");
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  return obdeect::stream::cli::run_tool(
      "obdeect propagation to-sky", obdeect::stream::cli::core_version(), argc, argv, usage, [&] {
        std::optional<std::string> profile;
        std::string input = "-";
        std::string output = "-";

        obdeect::stream::cli::Args args(argc, argv);
        std::string arg;
        while (args.next(arg)) {
          if (arg == "--profile") {
            profile = args.value("--profile");
          } else if (arg == "-i") {
            input = args.value("-i");
          } else if (arg == "-o") {
            output = args.value("-o");
          } else if (arg == "-h" || arg == "--help") {
            usage();
            return 0;
          } else {
            throw obdeect::stream::cli::UsageError("unknown argument: " + arg);
          }
        }

        std::optional<AtmosphereProfile> atmosphere;
        if (profile) {
          atmosphere.emplace(AtmosphereProfile::from_file(*profile));
        }
        if (input != "-" && output != "-" && same_file(input, output)) {
          throw obdeect::stream::cli::UsageError("input and output are the same file: " + input);
        }
        std::ifstream in_file;
        if (input != "-") {
          in_file.open(input, std::ios::binary);
          if (!in_file) {
            throw std::runtime_error("cannot open input stream: " + input);
          }
        }
        std::ofstream out_file;
        if (output != "-") {
          out_file.open(output, std::ios::binary);
          if (!out_file) {
            throw std::runtime_error("cannot open output stream: " + output);
          }
        }
        // std::cin is read one character at a time by libc++ once sync_with_stdio is off.
        obdeect::stream::StdinBuffer stdin_buffer;
        std::istream stdin_stream(&stdin_buffer);
        std::istream &in = input == "-" ? stdin_stream : static_cast<std::istream &>(in_file);
        std::ostream &out = output == "-" ? static_cast<std::ostream &>(std::cout)
                                          : static_cast<std::ostream &>(out_file);
        const ToSkySummary summary =
            ground_to_sky_stream(in, out, atmosphere, "obdeect propagation to-sky");
        out.flush();
        std::fprintf(stderr,
                     "obdeect propagation to-sky: photons in %llu, out %llu, dropped %llu\n",
                     static_cast<unsigned long long>(summary.photons_in),
                     static_cast<unsigned long long>(summary.photons_out),
                     static_cast<unsigned long long>(summary.dropped));
        return 0;
      });
}
