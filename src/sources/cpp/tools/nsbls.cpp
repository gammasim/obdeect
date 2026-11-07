// dlamb-nsbls: night-sky background photons on one telescope sphere, as sim_telarray
// LightEmission nsbls. Writes a gstrickd photon stream (sky frame) to stdout or -o FILE.

#include "dlamb/programs/nsbls.hpp"

#include "sky_args.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

namespace {

const char *usage_text =
    "Usage: dlamb-nsbls [options]\n"
    "Night-sky background photons on the telescope sphere; photon stream on stdout.\n"
    "Options (original nsbls.cc option, unit, default):\n"
    "  -r, --radius, --telescope-radius rtel_m  fiducial radius [m] (15)\n"
    "  -v, --viewcone deg    half opening angle of the emission cone [deg] (10)\n"
    "  -a, --altitude, --obslev, --observation-level m  observation level [m] (2150)\n"
    "  -A, --atmosphere n    profile number, recorded only (-1 = none)\n"
    "  --photons n           photons per event (1e6)\n"
    "  -b, --bs, --bunch-size f  photons per bunch (2; values below 1 are used as 1)\n"
    "  --events n            events (1)\n"
    "  --run n               run number (1)\n"
    "  -s, --spectrum s      wavelength [nm] (200..2000) or table file, default\n"
    "                        Benn_LaPalma_sky_converted.lis read from the working directory\n"
    "                        (copy it from sim_telarray/cfg/common or give -s)\n"
    "  -w, --pulse-width w   width [ns] of the flat time distribution (0). As the original,\n"
    "                        the width enters divided by 2.35482\n"
    "  -d, --pulse-delay d   delay of emission times [ns] (0)\n"
    "  --seed n              random seed (1, new option)\n"
    "  -o, --output fname    output file (stdout)\n"
    "  --help                this text\n";

int run(int argc, char **argv) {
  dlamb::programs::NsbOptions o;
  std::string output;
  std::string spectrum_text = "Benn_LaPalma_sky_converted.lis";
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
    } else if (arg == "-r" || arg == "--radius" || arg == "--telescope-radius") {
      o.telescope_radius_m = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-v" || arg == "--viewcone") {
      o.viewcone_deg = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-a" || arg == "--altitude" || arg == "--obslev" ||
               arg == "--observation-level") {
      o.altitude_m = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-A" || arg == "--atmosphere") {
      o.atmosphere = dlamb::cli::to_int(next(), arg);
    } else if (arg == "--photons") {
      o.photons = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-b" || arg == "--bs" || arg == "--bunch-size") {
      o.bunch_size = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--events") {
      o.events = dlamb::cli::to_count(next(), arg);
    } else if (arg == "--run") {
      o.run_number = dlamb::cli::to_int(next(), arg);
    } else if (arg == "-s" || arg == "--spectrum") {
      spectrum_text = next();
    } else if (arg == "-w" || arg == "--pulse-width") {
      o.pulse_width_ns = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-d" || arg == "--pulse-delay") {
      o.pulse_delay_ns = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--seed") {
      o.seed = static_cast<std::uint64_t>(dlamb::cli::to_int(next(), arg));
    } else if (arg == "-o" || arg == "--output") {
      output = next();
    } else {
      throw std::invalid_argument("unknown option " + arg + " (see --help)");
    }
  }
  if (o.atmosphere < -1)
    throw std::invalid_argument("--atmosphere: profile number (or -1 for none)");
  o.spectrum = dlamb::cli::parse_spectrum(spectrum_text, "--spectrum");

  if (output.empty()) {
    std::ios::sync_with_stdio(false);
    dlamb::programs::write_nsbls(std::cout, o);
    std::cout.flush();
  } else {
    std::ofstream file(output, std::ios::binary);
    if (!file)
      throw std::invalid_argument("cannot write " + output);
    dlamb::programs::write_nsbls(file, o);
  }
  return 0;
}

} // namespace

int main(int argc, char **argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "dlamb-nsbls: %s\n", error.what());
    return 1;
  }
}
