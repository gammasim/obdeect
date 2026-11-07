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
#include <vector>

namespace {

void usage(std::FILE *to) {
  std::fprintf(
      to,
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
      "  -p, --momentum GEV                  momentum (sets beta from the mass; not with --beta)\n"
      "  --primary ID|NAME                   CORSIKA id or gamma, electron, muon, ... "
      "(sets mass, charge)\n"
      "  -N, --events N                      number of events [1]\n"
      "  --run N                             run number [1]\n"
      "  --seed S                            random seed [1]\n"
      "  -o, --output FILE                   output file (default stdout)\n"
      "  --help                              this text\n"
      "Not ported: multiple scattering (--enable-scatter), --rscat, -I.\n");
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  try {
    const auto parsed =
        dlamb::fake_muon::parse_command_line(std::vector<std::string>(argv + 1, argv + argc));
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
