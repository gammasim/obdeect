// dlamb-ls-geo-test: geometry checks of LightEmission ls-geo-test. The original emits no
// photons; this tool writes its report to stdout and no photon stream.

#include "dlamb/programs/ls_geo_test.hpp"

#include <cstdio>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help") {
      std::printf("usage: dlamb-ls-geo-test\n"
                  "Prints the geometry report of LightEmission ls-geo-test (no photons).\n");
      return 0;
    }
    std::fprintf(stderr, "dlamb-ls-geo-test: unknown option %s\n", arg.c_str());
    return 1;
  }
  dlamb::ls_geo_test::write_report(std::cout);
  return 0;
}
