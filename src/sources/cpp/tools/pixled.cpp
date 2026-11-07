// dlamb-pixled: LED photon stream for selected camera pixels (port of LightEmission
// pixled). Writes a gstrickd stream to stdout or to -o FILE; --list-only prints the pixel
// combinations instead.

#include "dlamb/programs/pixled.hpp"

#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

using namespace dlamb::programs;

namespace {

void usage(std::ostream &out, const char *name) {
  out << "usage: " << name << " --camera-file FILE (--use LIST | --use-trg-rand N | "
      << "--use-trg-all) [options]\n"
      << "Options (default in braces; names as in pixled):\n"
      << "  --camera-file FILE        sim_telarray camera definition (required)\n"
      << "  --use LIST                pixels, '+' required, [a,b] pre-sum slaves\n"
      << "  --use-trg-rand N          N random trigger groups of the camera file\n"
      << "  --use-trg-all             all trigger groups of the camera file\n"
      << "  --require N               pixels to fire per event; 0 = all {0}\n"
      << "  --random N                N random distinct pixels per event (uses --require)\n"
      << "  --slaves                  also fire the pre-sum slaves\n"
      << "  --events N                events per pixel combination {1}\n"
      << "  --photons N               photons per pixel {500}\n"
      << "  --bunchsize F             photons per bunch (at least 1) {1}\n"
      << "  --camera-radius R         fiducial camera radius [cm] {123}\n"
      << "  --spectrum WLEN|Cherenkov:A:B  wavelength [nm] {400}\n"
      << "  --lightpulse SPEC         pulse, as in ff-1m {Simple:0}\n"
      << "  --run N                   run number {1}\n"
      << "  --altitude M              observation level [m] {2150}\n"
      << "  --atmosphere N            atmospheric profile number {1}\n"
      << "  --seed N                  random seed (new) {1}\n"
      << "  --list-only, -L           print the pixel combinations, no photons\n"
      << "  --quiet, -q               accepted; this program prints nothing else\n"
      << "  -o FILE                   output file {stdout}\n"
      << "  -h, --help                this text\n";
}

} // namespace

int main(int argc, char **argv) {
  pixled::Options options;
  std::string camera_file;
  std::string output;
  bool list_only = false;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "-h" || arg == "--help") {
        usage(std::cout, argv[0]);
        return 0;
      } else if (arg == "--camera-file" || arg == "-c") {
        camera_file = option_value(argc, argv, i);
      } else if (arg == "--use" || arg == "-u") {
        options.use = option_value(argc, argv, i);
      } else if (arg == "--use-trg-rand") {
        options.trigger_random = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--use-trg-all") {
        options.trigger_all = true;
      } else if (arg == "--require" || arg == "--required") {
        options.require = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--random") {
        options.random = true;
        options.require = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--slaves") {
        options.slaves = true;
      } else if (arg == "--events") {
        options.events = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--photons") {
        options.photons = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--bunchsize" || arg == "-b") {
        options.bunch_size = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--camera-radius" || arg == "-r") {
        options.camera_radius_cm = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--spectrum") {
        options.spectrum = option_value(argc, argv, i);
      } else if (arg == "--lightpulse" || arg == "--pulse") {
        options.lightpulse = option_value(argc, argv, i);
      } else if (arg == "--run" || arg == "-R") {
        options.run_number = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--altitude" || arg == "-a") {
        options.altitude_m = parse_real(option_value(argc, argv, i), arg);
      } else if (arg == "--atmosphere" || arg == "-A") {
        options.atmosphere = static_cast<int>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--seed") {
        options.seed = static_cast<std::uint64_t>(parse_integer(option_value(argc, argv, i), arg));
      } else if (arg == "--list-only" || arg == "-L") {
        list_only = true;
      } else if (arg == "--quiet" || arg == "-q") {
        continue;
      } else if (arg == "-o") {
        output = option_value(argc, argv, i);
      } else {
        bad("unknown option " + arg);
      }
    }
    if (camera_file.empty()) {
      bad("no camera definition file defined (--camera-file)");
    }
    std::ifstream camera(camera_file);
    if (!camera) {
      bad("cannot open camera configuration file " + camera_file);
    }
    options.camera = pixled::read_camera(camera, options.slaves);

    if (list_only) {
      for (const auto &event : pixled::build_sources(options)) {
        pixled::print_labels(std::cout, event);
      }
    } else {
      write_output(output, [&](std::ostream &out) { pixled::write_stream(out, options); });
    }
  } catch (const std::exception &error) {
    std::cerr << argv[0] << ": " << error.what() << "\n";
    return 1;
  }
  return 0;
}
