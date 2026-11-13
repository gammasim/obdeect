// dlamb-laser-beam: primary photons of a laser beam inside a telescope's field of view.

#include "obdeect/sources/programs/laser-beam.hpp"

#include "sky_args.hpp"
#include "tool_args.hpp"

#include <iostream>

namespace {

const char *usage = R"(usage: dlamb-laser-beam [options] [-o FILE]
  --events n  --step m  --bunchsize f  --primary-per-m f
  --spectrum s  --lightpulse p  --angular text
  --run n  --altitude m  --atmosphere n
  --laser-position x,y[,z]  --laser-theta deg  --laser-phi deg (astronomical azimuth)
  --telescope-theta deg  --telescope-phi deg  --telescope-radius m  --fov deg
  --seed n  -o FILE  --help
)";

// Astronomical azimuth to the mathematical sense.
double mathematical_azimuth(double astronomical) {
  return astronomical <= 0.0 ? -astronomical : 360.0 - astronomical;
}

} // namespace

int main(int argc, char **argv) {
  using namespace obdeect::sources;
  return cli::run_tool("dlamb-laser-beam", [&] {
    cli::ArgCursor cursor(argc, argv);
    programs::LsBeamOptions o;
    std::string spectrum = "355", pulse = "Simple:0", output, arg;
    while (cursor.next(arg)) {
      if (arg == "--help") {
        std::cout << usage;
        return 0;
      }
      const auto value = [&] { return cursor.value(arg); };
      if (arg == "--events")
        o.events = cli::to_count(value(), arg);
      else if (arg == "--step" || arg == "--steplength")
        o.step_m = cli::to_double(value(), arg);
      else if (arg == "--bunchsize")
        o.bunch_size = cli::to_double(value(), arg);
      else if (arg == "--primary-per-m")
        o.primary_photons_per_m = cli::to_double(value(), arg);
      else if (arg == "--brightness")
        throw std::invalid_argument("--brightness is replaced by --primary-per-m");
      else if (arg == "--spectrum")
        spectrum = value();
      else if (arg == "--lightpulse" || arg == "--pulse")
        pulse = value();
      else if (arg == "--angular" || arg == "--angular-distribution")
        o.angular = value();
      else if (arg == "--run")
        o.run_number = cli::parse_integer(value(), arg);
      else if (arg == "--altitude")
        o.altitude_m = cli::to_double(value(), arg);
      else if (arg == "--atmosphere")
        o.atmosphere = cli::parse_integer(value(), arg);
      else if (arg == "--laser-position") {
        const std::string text = value();
        const auto items = cli::split_items(text, arg);
        if (items.size() != 2 && items.size() != 3)
          throw std::invalid_argument("--laser-position needs x,y or x,y,z");
        o.laser_position_m = {cli::parse_number(items[0], arg), cli::parse_number(items[1], arg),
                              items.size() == 3 ? cli::parse_number(items[2], arg) : 0.0};
      } else if (arg == "--laser-theta")
        o.laser_theta_deg = cli::to_double(value(), arg);
      else if (arg == "--laser-phi")
        o.laser_phi_deg = mathematical_azimuth(cli::to_double(value(), arg));
      else if (arg == "--telescope-theta")
        o.telescope_theta_deg = cli::to_double(value(), arg);
      else if (arg == "--telescope-phi")
        o.telescope_phi_deg = mathematical_azimuth(cli::to_double(value(), arg));
      else if (arg == "--telescope-radius")
        o.telescope_radius_m = cli::to_double(value(), arg);
      else if (arg == "--fov")
        o.fov_deg = cli::to_double(value(), arg);
      else if (arg == "--seed")
        o.seed = cli::parse_seed(value(), arg);
      else if (arg == "-o")
        output = value();
      else
        throw std::invalid_argument("unknown option " + arg + " (see --help)");
    }
    o.spectrum = cli::parse_spectrum(spectrum, "--spectrum");
    o.pulse = cli::parse_pulse(pulse, "--lightpulse", false);
    cli::Output out(output);
    programs::write_laser_beam(out.stream(), o);
    out.close();
    return 0;
  });
}
