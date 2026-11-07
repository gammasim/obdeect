// dlamb-fake-muon: Cherenkov light of a single muon-like track as a gstrickd photon stream
// (port of LightEmission fake-muon). Writes to stdout, or to -o FILE.

#include "dlamb/programs/fake_muon.hpp"
#include "dlamb/cli.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void usage(std::FILE *to) {
  std::fprintf(to,
               "usage: dlamb-fake-muon [options]\n"
               "Writes a gstrickd photon stream (sky frame) to stdout, or to -o FILE.\n"
               "All lengths in metres, angles in degrees, wavelengths in nm.\n"
               "  -r, --radius, --telescope-radius R  telescope sphere radius [15]\n"
               "  -x X, -y Y                          impact position [0, 0]\n"
               "  -h, --h-start H                     starting altitude above obs. level [3000]\n"
               "  --theta DEG, --zenith-angle DEG     track zenith angle [0]\n"
               "  --phi DEG, --orientation, --azimuth track azimuth angle [0]\n"
               "  -n, --n1, --refraction N            index of refraction (N or N-1) [1.000292]\n"
               "  -A, --atmosphere const|-1           constant index (default)\n"
               "  --index-profile FILE                altitude [m], n-1 table (two columns)\n"
               "  -a, --altitude, --obslev H          observation level [2150]\n"
               "  -b, --bs, --bunch-size B            accepted, no effect on the photons [5]\n"
               "  -l, --wlmin NM, -L, --wlmax NM      emission band [240, 1000]\n"
               "  -q, --charge Q                      charge number [-1]\n"
               "  --beta B                            speed in units of c, (0, 1] [1]\n"
               "  -m, --mass GEV                      particle mass [0.105658]\n"
               "  -p, --momentum GEV                  momentum (sets beta from the mass)\n"
               "  --primary ID|NAME                   CORSIKA id or gamma, electron, muon, ... "
               "(sets mass, charge)\n"
               "  -N, --events N                      number of events [1]\n"
               "  --run N                             run number [1]\n"
               "  --seed S                            random seed [1]\n"
               "  -o, --output FILE                   output file (default stdout)\n"
               "  --help                              this text\n"
               "Not ported: multiple scattering (--enable-scatter), --rscat, -I.\n");
}

struct Parsed {
  dlamb::fake_muon::Options options;
  std::string output;
  bool help{false};
};

Parsed parse(int argc, char **argv) {
  Parsed parsed;
  auto &o = parsed.options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help") {
      parsed.help = true;
      return parsed;
    }
    if (arg == "--enable-scatter" || arg == "--scatter") {
      o.multiple_scattering = true;
      continue;
    }
    if (arg == "--disable-scatter" || arg == "--no-scatter") {
      o.multiple_scattering = false;
      continue;
    }
    if (arg.size() > 1 && arg[0] == '-' && arg[1] == 'I')
      throw std::invalid_argument("-I (atmosphere table search path) is not supported");
    if (i + 1 >= argc)
      throw std::invalid_argument("option " + arg + " needs a value");
    const std::string value = argv[++i];
    using dlamb::cli::parse_integer;
    using dlamb::cli::parse_number;
    if (arg == "-r" || arg == "--radius" || arg == "--telescope-radius")
      o.telescope_radius_m = parse_number(arg, value);
    else if (arg == "-x")
      o.impact_x_m = parse_number(arg, value);
    else if (arg == "-y")
      o.impact_y_m = parse_number(arg, value);
    else if (arg == "--nscat")
      o.nscat = static_cast<int>(parse_integer(arg, value));
    else if (arg == "--rscat")
      o.rscat_m = parse_number(arg, value);
    else if (arg == "-h" || arg == "--h-start")
      o.h_start_m = parse_number(arg, value);
    else if (arg == "-N" || arg == "--events")
      o.events = dlamb::fake_muon::event_count(parse_integer(arg, value));
    else if (arg == "-n" || arg == "--n1" || arg == "--refraction") {
      const double v = parse_number(arg, value);
      if (v > 1.0 && v < 1.0003)
        o.index_minus_one = v - 1.0;
      else if (v > 0.0 && v < 0.0003)
        o.index_minus_one = v;
      else
        throw std::invalid_argument("invalid refraction index: " + value);
    } else if (arg == "-A" || arg == "--atmosphere") {
      if (value != "const" && value != "-1")
        throw std::invalid_argument("atmosphere tables are not part of dlamb; use "
                                    "--index-profile FILE or --atmosphere const");
    } else if (arg == "--index-profile") {
      std::ifstream in(value);
      if (!in)
        throw std::invalid_argument("cannot open index profile " + value);
      o.index_profile = dlamb::fake_muon::read_index_profile(in, value);
    } else if (arg == "-a" || arg == "--altitude" || arg == "--obslev" ||
               arg == "--observation-level")
      o.observation_level_m = parse_number(arg, value);
    else if (arg == "-b" || arg == "--bs" || arg == "--bunch-size")
      o.bunch_size = parse_number(arg, value);
    else if (arg == "--run")
      o.run_number = static_cast<int>(parse_integer(arg, value));
    else if (arg == "--seed")
      o.seed = static_cast<std::uint64_t>(parse_integer(arg, value));
    else if (arg == "--primary") {
      if (!value.empty() && value[0] >= '0' && value[0] <= '9') {
        o.primary_id = static_cast<int>(parse_integer(arg, value));
      } else {
        o.primary_id = dlamb::fake_muon::primary_id_from_name(value);
      }
      dlamb::fake_muon::mass_charge(o.primary_id, o.mass_gev, o.charge);
    } else if (arg == "-l" || arg == "--wlmin")
      o.wavelength_min_nm = parse_number(arg, value);
    else if (arg == "-L" || arg == "--wlmax")
      o.wavelength_max_nm = parse_number(arg, value);
    else if (arg == "-q" || arg == "--charge")
      o.charge = parse_number(arg, value);
    else if (arg == "-m" || arg == "--mass")
      o.mass_gev = parse_number(arg, value);
    else if (arg == "-p" || arg == "--momentum") {
      const double p = parse_number(arg, value);
      const double ratio = p / o.mass_gev;
      o.beta = std::sqrt(ratio * ratio / (1.0 + ratio * ratio));
    } else if (arg == "--beta")
      o.beta = parse_number(arg, value);
    else if (arg == "--theta" || arg == "--zenith-angle")
      o.theta_deg = parse_number(arg, value);
    else if (arg == "--phi" || arg == "--orientation" || arg == "--azimuth")
      o.phi_deg = parse_number(arg, value);
    else if (arg == "-o" || arg == "--output")
      parsed.output = value;
    else
      throw std::invalid_argument("unknown option " + arg);
  }
  return parsed;
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  try {
    const Parsed parsed = parse(argc, argv);
    if (parsed.help) {
      usage(stdout);
      return 0;
    }
    dlamb::StreamSummary summary;
    if (parsed.output.empty()) {
      summary = dlamb::fake_muon::write(std::cout, parsed.options);
      std::cout.flush();
    } else {
      std::ofstream file(parsed.output, std::ios::binary);
      if (!file)
        throw std::runtime_error("cannot open " + parsed.output);
      summary = dlamb::fake_muon::write(file, parsed.options);
    }
    std::fprintf(stderr, "dlamb-fake-muon: expected weight %.6g, wrote %llu photons, weight %.6g\n",
                 summary.expected, static_cast<unsigned long long>(summary.photons),
                 summary.weight);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "dlamb-fake-muon: %s\n", error.what());
    usage(stderr);
    return 1;
  }
  return 0;
}
