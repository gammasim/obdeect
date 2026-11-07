// dlamb-ff-gct: flasher photon stream for the camera corners of a dual-mirror telescope
// (port of LightEmission ff-gct). Writes a gstrickd stream to stdout or to -o FILE.

#include "dlamb/programs/ff-gct.hpp"

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>

using namespace dlamb::programs;

namespace {

void usage(std::ostream &out, const char *name) {
  out << "usage: " << name << " [options]\n"
      << "Options (default in braces; names as in ff-gct):\n"
      << "  --events N                events {1}\n"
      << "  --photons N[,...]         photons emitted into 4 pi per flasher {2.5e6}\n"
      << "  --bunchsize F             photons per bunch (at least 1) {1}\n"
      << "  --flasher-xy X            flasher offset on the diagonal [cm] {14.62}\n"
      << "  --flasher-depth D         flasher depth behind the camera front [cm] {2.46}\n"
      << "  --flasher-inclination A   outward inclination [deg] {0}\n"
      << "  --mirror-camera-distance D  focal front to secondary [cm] {51.06949}\n"
      << "  --spectrum WLEN|Cherenkov:A:B  wavelength [nm] {400}\n"
      << "  --lightpulse SPEC         pulse, as in ff-1m {Simple:0}\n"
      << "  --angular-distribution SPEC  angular profile, as in ff-1m {isotropic}\n"
      << "  --delays T1,T2,T3,T4[,T5] delays of the flashers, and of the common flash [ns] "
         "{0,0,0,0,0}\n"
      << "  --fire all|alternate|demo|1..4  flashers to fire {all}\n"
      << "  --run N                   run number {1}\n"
      << "  --altitude M              observation level [m] {2100}\n"
      << "  --atmosphere N            atmospheric profile number {1}\n"
      << "  --seed N                  random seed (new) {1}\n"
      << "  -o FILE                   output file {stdout}\n"
      << "  -h, --help                this text\n";
}

} // namespace

int main(int argc, char **argv) {
  ff_gct::Options options;
  std::string output;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "-h" || arg == "--help") {
        usage(std::cout, argv[0]);
        return 0;
      } else if (arg == "--events") {
        options.events = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--photons") {
        options.photons = parse_real_list(option_value(argc, argv, i), arg);
      } else if (arg == "--bunchsize") {
        options.bunch_size = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--xy" || arg == "--flasher-xy") {
        options.flasher_xy_cm = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "-z" || arg == "--z" || arg == "--flasher-depth") {
        options.flasher_depth_cm = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--inclination" || arg == "--flasher-inclination") {
        options.inclination_deg = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--mirror-distance" || arg == "--mirror-camera-distance") {
        options.mirror_camera_distance_cm = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--spectrum") {
        options.spectrum = option_value(argc, argv, i);
      } else if (arg == "--lightpulse" || arg == "--pulse") {
        options.lightpulse = option_value(argc, argv, i);
      } else if (arg == "--angular-distribution" || arg == "--angular") {
        options.angular = option_value(argc, argv, i);
      } else if (arg == "--delays") {
        const auto delays = parse_real_list(option_value(argc, argv, i), arg);
        if (delays.size() < 4 || delays.size() > 5) {
          bad("--delays needs 4 or 5 values");
        }
        options.delays = {0.0, 0.0, 0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < delays.size(); ++k) {
          options.delays[k] = delays[k];
        }
      } else if (arg == "--fire") {
        options.fire = option_value(argc, argv, i);
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
    write_output(output, [&](std::ostream &out) { ff_gct::write_stream(out, options); });
  } catch (const std::exception &error) {
    std::cerr << argv[0] << ": " << error.what() << "\n";
    return 1;
  }
  return 0;
}
