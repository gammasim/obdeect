// dlamb-rotating-led: rotating focal-plane LED photon stream (port of LightEmission/fpls.cc).
// Writes a gstrickd stream to stdout or to -o FILE.

#include "dlamb/programs/rotating-led.hpp"

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>

using namespace dlamb::programs;

namespace {

void usage(std::ostream &out, const char *name) {
  out << "usage: " << name << " [options]\n"
      << "The original fpls.cc has no options; these are its constants and defaults and defaults:\n"
      << "  --altitude M              observation level [m] {2100}\n"
      << "  --camera-radius R         fiducial camera radius [m]; LED at 0.7 R {1.4}\n"
      << "  --photons N               photons per time step {5000}\n"
      << "  --bunchsize F             photons per bunch (at least 1) {1}\n"
      << "  --rotation-period T       rotation period [ns] {60}\n"
      << "  --time-step DT            time step [ns] {0.25}\n"
      << "  --spectrum WLEN|Cherenkov:A:B  wavelength [nm] {400} (original: ls_spec.lis)\n"
      << "  --lightpulse SPEC         pulse, as in camera-flasher {Simple:0}\n"
      << "  --angular-distribution SPEC  beam profile {isotropic} (original: ls_beamed.lis)\n"
      << "  --run N                   run number {1}\n"
      << "  --atmosphere N            atmospheric profile number {1}\n"
      << "  --seed N                  random seed (new) {1}\n"
      << "  -o FILE                   output file {stdout}\n"
      << "  -h, --help                this text\n";
}

} // namespace

int main(int argc, char **argv) {
  rotating_led::Options options;
  std::string output;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "-h" || arg == "--help") {
        usage(std::cout, argv[0]);
        return 0;
      } else if (arg == "--altitude" || arg == "-a") {
        options.altitude_m = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--camera-radius" || arg == "-r") {
        options.camera_radius_m = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--photons") {
        options.photons = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--bunchsize") {
        options.bunch_size = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--rotation-period") {
        options.rotation_period_ns = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--time-step") {
        options.time_step_ns = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--spectrum") {
        options.spectrum = option_value(argc, argv, i);
      } else if (arg == "--lightpulse" || arg == "--pulse") {
        options.lightpulse = option_value(argc, argv, i);
      } else if (arg == "--angular-distribution" || arg == "--angular") {
        options.angular = option_value(argc, argv, i);
      } else if (arg == "--run") {
        options.run_number = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
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
    write_output(output, [&](std::ostream &out) { rotating_led::write_stream(out, options); });
  } catch (const std::exception &error) {
    std::cerr << argv[0] << ": " << error.what() << "\n";
    return 1;
  }
  return 0;
}
