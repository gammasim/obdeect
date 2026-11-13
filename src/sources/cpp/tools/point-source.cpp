// obdeect sources point-source: point light source above the observation level (sky frame).

#include "obdeect/sources/programs/point-source.hpp"

#include "sky_args.hpp"
#include "telpos_file.hpp"
#include "tool_args.hpp"

#include <fstream>
#include <iostream>

namespace {

const char *usage = R"(usage: obdeect sources point-source [options] [-o FILE]
  -R run     run number (1)              -N n      events per case (1)
  -r r       telescope radius [cm]       -x X -y Y source position [cm]
  -z z       height [cm]                 -n N[,..] photons into 4 pi, one case each
  -b b       photons per bunch           -s spec   wavelength, Cherenkov:a:b or table file
  -p pulse   Gauss:s, Simple:w, ... or bare number (Gauss sigma)
  -a ang     isotropic, Cone:a, Gauss:s, Rayleigh, ... or table file
  -h alt     observation level [m]       -A n      atmosphere id
  -d x,y,z   pointing                    --telpos-file F, --corsika-telpos-file F
  --seed n   -o FILE   --help
)";

} // namespace

int main(int argc, char **argv) {
  using namespace obdeect::sources;
  return cli::run_tool("obdeect sources point-source", argc, argv, [&] {
    cli::Args cursor(argc, argv);
    programs::XyzlsOptions o;
    std::string spectrum = "400", pulse = "Gauss:0.0", angular = "isotropic", output, telpos,
                corsika_telpos;
    std::string telescope_option;
    std::string arg;
    while (cursor.next(arg)) {
      if (arg == "--help") {
        std::cout << usage;
        return 0;
      }
      const auto value = [&] { return cursor.value(arg); };
      if (arg == "-R")
        o.run_number = cli::parse_integer(value(), arg);
      else if (arg == "-N")
        o.events = cli::to_count(value(), arg);
      else if (arg == "-r")
        o.telcam_radius_cm = cli::parse_number(value(), arg);
      else if (arg == "-x")
        o.x_cm = cli::parse_number(value(), arg);
      else if (arg == "-y")
        o.y_cm = cli::parse_number(value(), arg);
      else if (arg == "-z")
        o.z_cm = cli::parse_number(value(), arg);
      else if (arg == "-n")
        o.photons = cli::parse_list(value(), arg);
      else if (arg == "-b")
        o.bunch_size = cli::parse_number(value(), arg);
      else if (arg == "-s")
        spectrum = value();
      else if (arg == "-p")
        pulse = value();
      else if (arg == "-a")
        angular = value();
      else if (arg == "-h")
        o.altitude_m = cli::parse_number(value(), arg);
      else if (arg == "-A")
        o.atmosphere = cli::parse_integer(value(), arg);
      else if (arg == "-d")
        o.direction = cli::parse_vector(value(), arg);
      else if (arg == "--telpos-file") {
        telpos = value();
        telescope_option = arg;
      } else if (arg == "--corsika-telpos-file") {
        telpos = value();
        telescope_option = arg;
      } else if (arg == "--seed")
        o.seed = cli::parse_seed(value(), arg);
      else if (arg == "-o")
        output = value();
      else
        throw std::invalid_argument("unknown option " + arg + " (see --help)");
    }
    o.spectrum = cli::parse_spectrum(spectrum, "-s");
    o.pulse = cli::parse_pulse(pulse, "-p", true);
    o.angular = cli::parse_angular(angular, "-a");
    if (!telescope_option.empty()) {
      std::ifstream file(telpos);
      if (!file)
        throw std::invalid_argument(telescope_option + ": cannot open telescope file " + telpos);
      o.telescopes = cli::read_telpos(file, telescope_option, telescope_option != "--telpos-file");
    }
    cli::Output out(output);
    programs::write_point_source(out.stream(), o);
    out.close();
    return 0;
  });
}
