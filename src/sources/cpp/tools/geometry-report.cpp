// dlamb-geometry-report: geometry checks of LightEmission/ls-geo-test.cc. The original emits no
// photons; this tool writes its report to stdout and no photon stream.

#include "dlamb/programs/geometry-report.hpp"

#include <cstdio>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help") {
      std::printf("usage: dlamb-geometry-report\n"
                  "Prints the geometry report of LightEmission/ls-geo-test.cc (no photons).\n");
      return 0;
    }
    std::fprintf(stderr, "dlamb-geometry-report: unknown option %s\n", arg.c_str());
    return 1;
  }
  dlamb::geometry_report::write_report(std::cout);
  return 0;
}
