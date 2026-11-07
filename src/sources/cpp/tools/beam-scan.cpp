// dlamb-beam-scan: parallel light beam scanned along two lines, as sim_telarray
// LightEmission/xyls.cc (test program). Writes a gstrickd photon stream (sky frame) to stdout or -o
// FILE.

#include "dlamb/programs/beam-scan.hpp"

#include "sky_args.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

namespace {

const char *usage_text =
    "Usage: dlamb-beam-scan [options]\n"
    "Parallel downward light beam scanned along x and y; photon stream on stdout.\n"
    "Options (original xyls.cc option, unit, default; defaults from r = 123 cm):\n"
    "  -r rcam     camera fiducial radius [cm] (123)\n"
    "  -x xpos     start of the scan [cm] (0.3 x 123)\n"
    "  -y ypos     start of the scan [cm] (0.2 x 123)\n"
    "  -h height   source height above the observation level [cm] (0.05 x 123)\n"
    "  -o fname    output file (stdout)\n"
    "  --seed n    random seed (1, new option)\n"
    "  --help      this text\n";

int run(int argc, char **argv) {
  dlamb::programs::XylsOptions o;
  std::string output;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto next = [&]() -> std::string {
      if (i + 1 >= argc)
        throw std::invalid_argument(arg + ": missing value");
      return argv[++i];
    };
    if (arg == "--help") {
      std::fputs(usage_text, stdout);
      return 0;
    } else if (arg == "-r") {
      o.radius_cm = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-x") {
      o.x_cm = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-y") {
      o.y_cm = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-h") {
      o.height_cm = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--seed") {
      o.seed = static_cast<std::uint64_t>(dlamb::cli::to_int(next(), arg));
    } else if (arg == "-o") {
      output = next();
    } else {
      throw std::invalid_argument("unknown option " + arg + " (see --help)");
    }
  }
  if (output.empty()) {
    std::ios::sync_with_stdio(false);
    dlamb::programs::write_beam_scan(std::cout, o);
    std::cout.flush();
  } else {
    std::ofstream file(output, std::ios::binary);
    if (!file)
      throw std::invalid_argument("cannot write " + output);
    dlamb::programs::write_beam_scan(file, o);
  }
  return 0;
}

} // namespace

int main(int argc, char **argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "dlamb-beam-scan: %s\n", error.what());
    return 1;
  }
}
