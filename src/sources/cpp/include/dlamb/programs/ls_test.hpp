#pragma once

// Port of LightEmission ls-test.cc: the octo figure-eight with 4e7 photons per position
// (ls-test.cc:53, passed as the photon number of Emit, ls-test.cc:87). The original
// does not build as written (main() uses argc and argv without parameters).

#include "dlamb/programs/octo.hpp"

#include <cstdint>
#include <ostream>
#include <string>

namespace dlamb::ls_test {

// Options of ls-test.cc. The file name (octo.dat.gz, ls-test.cc:49) and the bunch size of
// 2 photons (ls-test.cc:87) do not change the output; see the octo header.
struct Options {
  double photons_per_position{4e7}; // nb_per_pos (ls-test.cc:53)
  std::uint64_t seed{1};            // new: explicit random seed
};

// Writes the complete stream (run, one event, sky frame, sphere at the origin of radius
// 12.5 m). Same geometry as octo.
inline StreamSummary write(std::ostream &out, const Options &o,
                           const std::string &producer = "dlamb-ls-test") {
  octo::Options octo_options;
  octo_options.photons_per_position = o.photons_per_position;
  octo_options.seed = o.seed;
  return octo::write(out, octo_options, producer);
}

} // namespace dlamb::ls_test
