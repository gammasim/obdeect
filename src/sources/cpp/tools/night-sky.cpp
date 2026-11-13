// obdeect sources night-sky: night-sky background photons on a telescope sphere (sky frame).

#include "obdeect/sources/programs/night-sky.hpp"

#include "sky_args.hpp"
#include "tool_args.hpp"

#include <iostream>

namespace {

const char *usage = R"(usage: obdeect sources night-sky [options] [-o FILE]
  -r, --radius, --telescope-radius m   telescope sphere radius (15)
  -v, --viewcone deg                   half opening angle (10)
  -a, --altitude, --obslev m           observation level (2150)
  -A, --atmosphere n                   profile number or -1
  --photons n  -b, --bs, --bunch-size f  --events n  --run n
  -s, --spectrum s                     wavelength, Cherenkov:a:b or table (default
                                       Benn_LaPalma_sky_converted.lis in the working directory)
  -w, --pulse-width ns  -d, --pulse-delay ns  --seed n  -o, --output FILE  --help
)";

} // namespace

int main(int argc, char **argv) {
  using namespace obdeect::sources;
  return cli::run_tool("obdeect sources night-sky", argc, argv, [&] {
    cli::Args cursor(argc, argv);
    programs::NsbOptions o;
    std::string spectrum = "Benn_LaPalma_sky_converted.lis", output, arg;
    while (cursor.next(arg)) {
      if (arg == "--help") {
        std::cout << usage;
        return 0;
      }
      const auto value = [&] { return cursor.value(arg); };
      if (arg == "-r" || arg == "--radius" || arg == "--telescope-radius")
        o.telescope_radius_m = cli::parse_number(value(), arg);
      else if (arg == "-v" || arg == "--viewcone")
        o.viewcone_deg = cli::parse_number(value(), arg);
      else if (arg == "-a" || arg == "--altitude" || arg == "--obslev" ||
               arg == "--observation-level")
        o.altitude_m = cli::parse_number(value(), arg);
      else if (arg == "-A" || arg == "--atmosphere") {
        o.atmosphere = cli::parse_integer(value(), arg);
        if (o.atmosphere < -1)
          throw std::invalid_argument("--atmosphere: profile number (or -1 for none)");
      } else if (arg == "--photons")
        o.photons = cli::parse_number(value(), arg);
      else if (arg == "-b" || arg == "--bs" || arg == "--bunch-size")
        o.bunch_size = cli::parse_number(value(), arg);
      else if (arg == "--events")
        o.events = cli::to_count(value(), arg);
      else if (arg == "--run")
        o.run_number = cli::parse_integer(value(), arg);
      else if (arg == "-s" || arg == "--spectrum")
        spectrum = value();
      else if (arg == "-w" || arg == "--pulse-width")
        o.pulse_width_ns = cli::parse_number(value(), arg);
      else if (arg == "-d" || arg == "--pulse-delay")
        o.pulse_delay_ns = cli::parse_number(value(), arg);
      else if (arg == "--seed")
        o.seed = cli::parse_seed(value(), arg);
      else if (arg == "-o" || arg == "--output")
        output = value();
      else
        throw std::invalid_argument("unknown option " + arg + " (see --help)");
    }
    o.spectrum = cli::parse_spectrum(spectrum, "--spectrum");
    cli::Output out(output);
    programs::write_night_sky(out.stream(), o);
    out.close();
    return 0;
  });
}
