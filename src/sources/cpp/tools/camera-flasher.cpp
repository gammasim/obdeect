// dlamb-camera-flasher: flat-fielding flasher photon stream (port of LightEmission/ff-1m.cc).
// Writes a gstrickd stream to stdout or to -o FILE.

#include "dlamb/programs/camera-flasher.hpp"

#include <exception>
#include <iostream>
#include <string>

using namespace dlamb::programs;

namespace {

void usage(std::ostream &out, const char *name) {
  out << "usage: " << name << " [options]\n"
      << "Options (default in braces; names as in camera-flasher):\n"
      << "  --events N[,...]          events per photon level {1}\n"
      << "  --photons N[,...]         photons emitted into 4 pi per flasher {2.5e6}\n"
      << "  --bunchsize F             photons per bunch (at least 1) {1}\n"
      << "  --xy X,Y                  flasher offset from the axis [cm] {0,0}\n"
      << "  --distance D              flasher distance from the camera front [cm] {1600}\n"
      << "  --camera-radius R         fiducial camera radius [cm] {120}\n"
      << "  --spectrum WLEN|Cherenkov:A:B  wavelength [nm] {400}\n"
      << "  --lightpulse SPEC         Simple:W[,D], Gauss:S, Gauss-fwhm:F, Exp:T, Complex:A,B,C[,D]"
         " {Simple:0}\n"
      << "  --angular-distribution SPEC  isotropic, Parallel, Cone:DEG, FilledCone:DEG, "
         "Gauss:DEG, Gauss-fwhm:DEG, Rayleigh {isotropic}\n"
      << "  --run N                   run number {1}\n"
      << "  --altitude M              observation level [m] {2150}\n"
      << "  --atmosphere N            atmospheric profile number {1}\n"
      << "  --seed N                  random seed (new) {1}\n"
      << "  -o FILE                   output file {stdout}\n"
      << "  -h, --help                this text\n";
}

} // namespace

int main(int argc, char **argv) {
  camera_flasher::Options options;
  std::string output;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "-h" || arg == "--help") {
        usage(std::cout, argv[0]);
        return 0;
      } else if (arg == "--events") {
        options.events = parse_integer_list(option_value(argc, argv, i), arg);
      } else if (arg == "--photons") {
        options.photons = parse_real_list(option_value(argc, argv, i), arg);
      } else if (arg == "--bunchsize" || arg == "--bs" || arg == "-b") {
        options.bunch_size = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--xy" || arg == "--flasher-xy") {
        const auto xy = parse_real_list(option_value(argc, argv, i), arg);
        if (xy.size() != 2) {
          bad("--xy needs x,y");
        }
        options.x_cm = xy[0];
        options.y_cm = xy[1];
      } else if (arg == "--distance" || arg == "-z") {
        options.distance_cm = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--camera-radius" || arg == "-r") {
        options.camera_radius_cm = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--spectrum") {
        options.spectrum = option_value(argc, argv, i);
      } else if (arg == "--lightpulse" || arg == "--pulse") {
        options.lightpulse = option_value(argc, argv, i);
      } else if (arg == "--angular-distribution" || arg == "--angular") {
        options.angular = option_value(argc, argv, i);
      } else if (arg == "--run") {
        options.run_number = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--altitude" || arg == "-a") {
        options.altitude_m = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--atmosphere" || arg == "-A") {
        options.atmosphere = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--seed") {
        options.seed = static_cast<std::uint64_t>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "-o") {
        output = option_value(argc, argv, i);
      } else {
        bad("unknown option " + arg);
      }
    }
    write_output(output, [&](std::ostream &out) { camera_flasher::write_stream(out, options); });
  } catch (const std::exception &error) {
    std::cerr << argv[0] << ": " << error.what() << "\n";
    return 1;
  }
  return 0;
}
