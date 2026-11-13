#pragma once

// Reader of CORSIKA/sim_telarray atmospheric tables (atmprof format), shared by the stream
// adapters and propagation. Header-only, standard library only.

#include "obdeect/stream/eventio.hpp"
#include "obdeect/stream/parse.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace obdeect::stream {

// Profile number of the custom table (CORSIKA ATMOFILE).
inline constexpr std::uint32_t custom_profile_id = 99;

// Parses an atmospheric table as CORSIKA and sim_telarray read it:
// '#' starts a comment, blank lines are skipped, and the first four columns are altitude
// [km a.s.l.], density [g/cm^3], thickness [g/cm^2] and n-1. Further columns are ignored.
// Altitudes must ascend; 2 to 120 levels. The profile number is
// custom_profile_id; the observation level is set by the caller. Throws
// std::invalid_argument with the line number.
inline AtmosphereTable read_atmprof_text(std::istream &in, const std::string &file_name) {
  AtmosphereTable table;
  table.profile_id = custom_profile_id;
  table.file_name = file_name;
  std::string line;
  std::size_t number = 0;
  const auto fail = [&number](const std::string &what) {
    return std::invalid_argument("atmospheric table line " + std::to_string(number) + ": " + what);
  };
  while (std::getline(in, line)) {
    ++number;
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') {
      continue;
    }
    std::istringstream fields(line);
    std::array<std::string, 4> tokens;
    for (std::string &token : tokens) {
      if (!(fields >> token)) {
        throw fail("expected 4 columns: altitude_km density thickness n_minus_1");
      }
    }
    std::array<double, 4> v{};
    for (std::size_t k = 0; k < v.size(); ++k) {
      if (!parse::finite_number(tokens[k], v[k])) {
        throw fail("not a finite number: " + tokens[k]);
      }
    }
    if (v[1] < 0.0 || v[2] < 0.0) {
      throw fail("density and thickness must not be negative");
    }
    if (!table.alt_km.empty() && !(v[0] > table.alt_km.back())) {
      throw fail("altitudes must ascend");
    }
    table.alt_km.push_back(v[0]);
    table.rho.push_back(v[1]);
    table.thick.push_back(v[2]);
    table.refidx_m1.push_back(v[3]);
  }
  if (table.alt_km.size() < 2) {
    throw std::invalid_argument("atmospheric table needs at least 2 levels");
  }
  if (table.alt_km.size() > 120) {
    throw std::invalid_argument("atmospheric table has more than 120 levels");
  }
  return table;
}

} // namespace obdeect::stream
