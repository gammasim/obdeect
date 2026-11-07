// dlamb-xyzls: point light source above the telescope, as sim_telarray LightEmission
// xyzls. Writes a gstrickd photon stream (sky frame) to stdout or to -o FILE.

#include "dlamb/programs/xyzls.hpp"

#include "sky_args.hpp"
#include "telpos_file.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

namespace {

const char *usage_text =
    "Usage: dlamb-xyzls [options]\n"
    "Point light source at (x, y, z) above the observation level; photon stream on stdout.\n"
    "Options (original xyzls.cc option, unit, default):\n"
    "  -R run          run number (1)\n"
    "  -N nevents      events per case (1)\n"
    "  -r rtelcam      fiducial radius of the default telescope [cm] (1500)\n"
    "  -x xpos         source x [cm] (0)\n"
    "  -y ypos         source y [cm] (0)\n"
    "  -z zpos         source height above the observation level [cm] (1000e2)\n"
    "  -n nphot[,...]  photons into 4 pi per event, one case per value (1e7)\n"
    "  -b bsize        photons per bunch (1; values below 1 are used as 1)\n"
    "  -s spectrum     nm (200..2000), Cherenkov:a:b, or table file (400)\n"
    "  -p pulse        bare number = Gauss sigma [ns], or Simple:w, Gauss:s, Gauss-fwhm:f,\n"
    "                  Exp:tau, Complex:wr,wg,we, or table file, with optional ,delay (Gauss:0.0)\n"
    "  -a angdist      isotropic, Gauss:s, Gauss-fwhm:f, Rayleigh, Cone:a, Parallel,\n"
    "                  FilledCone:a [deg], or table file (isotropic)\n"
    "  -h alt_m        observation level a.s.l. [m] (2150)\n"
    "  -A atmosphere   profile number, recorded in the run keys (1)\n"
    "  -d cx,cy,cz     pointing (0,0,-1), normalised (0,0,-1)\n"
    "  --telpos-file f      telescopes 'x y z r' per line [cm] (one telescope at origin)\n"
    "  --corsika-telpos-file f  TELESCOPE lines of CORSIKA inputs [cm]\n"
    "  --seed n        random seed (1, new option)\n"
    "  -o fname        output file (stdout)\n"
    "  --help          this text\n";

int run(int argc, char **argv) {
  dlamb::programs::XyzlsOptions o;
  std::string output;
  std::string pulse_text = "Gauss:0.0";
  std::string spectrum_text = "400";
  std::string angular_text = "isotropic";
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
    } else if (arg == "-R") {
      o.run_number = dlamb::cli::to_int(next(), arg);
    } else if (arg == "-N") {
      o.events = dlamb::cli::to_count(next(), arg);
    } else if (arg == "-r") {
      o.telcam_radius_cm = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-x") {
      o.x_cm = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-y") {
      o.y_cm = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-z") {
      o.z_cm = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-n") {
      o.photons = dlamb::cli::to_list(next(), arg);
    } else if (arg == "-b") {
      o.bunch_size = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-s") {
      spectrum_text = next();
    } else if (arg == "-p") {
      pulse_text = next();
    } else if (arg == "-a") {
      angular_text = next();
    } else if (arg == "-h") {
      o.altitude_m = dlamb::cli::to_double(next(), arg);
    } else if (arg == "-A") {
      o.atmosphere = dlamb::cli::to_int(next(), arg);
    } else if (arg == "-d") {
      o.direction = dlamb::cli::to_vec3(next(), arg);
    } else if (arg == "--telpos-file") {
      o.telescopes = dlamb::cli::read_telpos(next(), arg, false);
    } else if (arg == "--corsika-telpos-file") {
      o.telescopes = dlamb::cli::read_telpos(next(), arg, true);
    } else if (arg == "--seed") {
      o.seed = static_cast<std::uint64_t>(dlamb::cli::to_int(next(), arg));
    } else if (arg == "-o") {
      output = next();
    } else {
      throw std::invalid_argument("unknown option " + arg + " (see --help)");
    }
  }
  // Numbers as pulse names are Gauss sigma, as in xyzls.cc (main: isdigit -> "Gauss:").
  o.pulse = dlamb::cli::parse_pulse(pulse_text, "-p", true);
  o.spectrum = dlamb::cli::parse_spectrum(spectrum_text, "-s");
  o.angular = dlamb::cli::parse_angular(angular_text, "-a");

  if (output.empty()) {
    std::ios::sync_with_stdio(false);
    dlamb::programs::write_xyzls(std::cout, o);
    std::cout.flush();
  } else {
    std::ofstream file(output, std::ios::binary);
    if (!file)
      throw std::invalid_argument("cannot write " + output);
    dlamb::programs::write_xyzls(file, o);
  }
  return 0;
}

} // namespace

int main(int argc, char **argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "dlamb-xyzls: %s\n", error.what());
    return 1;
  }
}
