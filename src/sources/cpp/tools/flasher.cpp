// dlamb-flasher: flasher photon stream. --type dish (default) is a port of LightEmission/ff-1m.cc,
// --type secondary a port of LightEmission/ff-gct.cc. Writes a gstrickd stream to stdout or -o
// FILE.

#include "dlamb/programs/flasher.hpp"

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace dlamb::programs;

namespace {

void usage(std::ostream &out, const char *name) {
  out << "usage: " << name << " [--type dish|secondary] [options]\n"
      << "Options (default in braces; names as in ff-1m.cc for dish, ff-gct.cc for secondary):\n"
      << "  --type dish|secondary     dish: flasher above the dish centre lighting the camera;\n"
      << "                            secondary: four camera-corner flashers onto the secondary "
         "{dish}\n"
      << "All types:\n"
      << "  --photons N[,...]         photons emitted into 4 pi per flasher {2.5e6}\n"
      << "  --bunchsize F             photons per bunch (at least 1) {1}\n"
      << "  --spectrum WLEN|Cherenkov:A:B  wavelength [nm] {400}\n"
      << "  --lightpulse SPEC         Simple:W[,D], Gauss:S, Gauss-fwhm:F, Exp:T, Complex:A,B,C[,D]"
         " {Simple:0}\n"
      << "  --angular-distribution SPEC  isotropic, Parallel, Cone:DEG, FilledCone:DEG, "
         "Gauss:DEG, Gauss-fwhm:DEG, Rayleigh {isotropic}\n"
      << "  --run N                   run number {1}\n"
      << "  --altitude M              observation level [m] {2150 dish, 2100 secondary}\n"
      << "  --atmosphere N            atmospheric profile number {1}\n"
      << "  --seed N                  random seed (new) {1}\n"
      << "  -o FILE                   output file {stdout}\n"
      << "  -h, --help                this text\n"
      << "--type dish only:\n"
      << "  --events N[,...]          events per photon level {1}; a short list repeats its last "
         "value\n"
      << "  --xy X,Y (--flasher-xy X,Y)  flasher offset from the axis [cm] {0,0}\n"
      << "  --distance D (-z D)       flasher distance above the telescope centre [cm] {1600}\n"
      << "  --camera-radius R (-r R)  fiducial camera radius [cm] {120}\n"
      << "--type secondary only:\n"
      << "  --events N                events {1}\n"
      << "  --flasher-xy X (--xy X)   flasher offset on the diagonal [cm] {14.62}\n"
      << "  --flasher-depth D (-z D, --z D)  flasher depth behind the camera front [cm] {2.46}\n"
      << "  --flasher-inclination A (--inclination A)  outward inclination [deg] {0}\n"
      << "  --mirror-camera-distance D (--mirror-distance D)  focal front to secondary [cm] "
         "{51.06949}\n"
      << "  --delays T1,T2,T3,T4[,T5] delays of the flashers and of the common flash [ns] "
         "{0,0,0,0,0}\n"
      << "  --fire all|alternate|demo|1..4  flashers to fire {all}\n";
}

} // namespace

int main(int argc, char **argv) {
  std::string output;
  try {
    // Collect name and value pairs first: the type decides how --xy, -z and the others apply.
    std::vector<std::pair<std::string, std::string>> items;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "-h" || arg == "--help") {
        usage(std::cout, argv[0]);
        return 0;
      }
      items.emplace_back(arg, option_value(argc, argv, i));
    }

    flasher::Type type = flasher::Type::dish;
    for (const auto &[arg, value] : items) {
      if (arg == "--type") {
        type = flasher::parse_type(value);
      }
    }

    flasher::Options options;
    options.type = type;
    const auto dish_only = [&](const std::string &arg) {
      flasher::require_type(type, flasher::Type::dish, arg);
    };
    const auto secondary_only = [&](const std::string &arg) {
      flasher::require_type(type, flasher::Type::secondary, arg);
    };

    for (const auto &[arg, value] : items) {
      if (arg == "--type") {
        continue;
      } else if (arg == "--events") {
        options.events = parse_integer_list(value, arg);
      } else if (arg == "--photons") {
        options.photons = parse_real_list(value, arg);
      } else if (arg == "--bunchsize" || arg == "--bs" || arg == "-b") {
        options.bunch_size = parse_real(value, arg);
      } else if (arg == "--xy" || arg == "--flasher-xy") {
        if (type == flasher::Type::dish) {
          const auto xy = parse_real_list(value, arg);
          if (xy.size() != 2) {
            bad(arg + " needs x,y with --type dish");
          }
          options.x_cm = xy[0];
          options.y_cm = xy[1];
        } else {
          options.flasher_xy_cm = parse_real(value, arg);
        }
      } else if (arg == "--distance") {
        dish_only(arg);
        options.distance_cm = parse_real(value, arg);
      } else if (arg == "-z") {
        if (type == flasher::Type::dish) {
          options.distance_cm = parse_real(value, arg);
        } else {
          options.flasher_depth_cm = parse_real(value, arg);
        }
      } else if (arg == "--camera-radius" || arg == "-r") {
        dish_only(arg);
        options.camera_radius_cm = parse_real(value, arg);
      } else if (arg == "--flasher-depth" || arg == "--z") {
        secondary_only(arg);
        options.flasher_depth_cm = parse_real(value, arg);
      } else if (arg == "--inclination" || arg == "--flasher-inclination") {
        secondary_only(arg);
        options.inclination_deg = parse_real(value, arg);
      } else if (arg == "--mirror-distance" || arg == "--mirror-camera-distance") {
        secondary_only(arg);
        options.mirror_camera_distance_cm = parse_real(value, arg);
      } else if (arg == "--delays") {
        secondary_only(arg);
        const auto delays = parse_real_list(value, arg);
        if (delays.size() < 4 || delays.size() > 5) {
          bad("--delays needs 4 or 5 values");
        }
        options.delays = {0.0, 0.0, 0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < delays.size(); ++k) {
          options.delays[k] = delays[k];
        }
      } else if (arg == "--fire") {
        secondary_only(arg);
        options.fire = value;
      } else if (arg == "--spectrum") {
        options.spectrum = value;
      } else if (arg == "--lightpulse" || arg == "--pulse") {
        options.lightpulse = value;
      } else if (arg == "--angular-distribution" || arg == "--angular") {
        options.angular = value;
      } else if (arg == "--run") {
        options.run_number = static_cast<int>(parse_integer(value, arg));
      } else if (arg == "--altitude" || arg == "-a") {
        options.altitude_m = parse_real(value, arg);
      } else if (arg == "--atmosphere" || arg == "-A") {
        options.atmosphere = static_cast<int>(parse_integer(value, arg));
      } else if (arg == "--seed") {
        options.seed = static_cast<std::uint64_t>(parse_integer(value, arg));
      } else if (arg == "-o") {
        output = value;
      } else {
        bad("unknown option " + arg);
      }
    }
    write_output(output, [&](std::ostream &out) { flasher::write_stream(out, options); });
  } catch (const std::exception &error) {
    std::cerr << argv[0] << ": " << error.what() << "\n";
    return 1;
  }
  return 0;
}
