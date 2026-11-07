#pragma once

// Telescope lists of dlamb-xyzls: --telpos-file (lines 'x y z r', cm) and
// --corsika-telpos-file (CORSIKA 'TELESCOPE x y z r' lines; other lines are CORSIKA input
// and are skipped). '#' starts a comment. A file without telescope records, or a telescope
// line that is not four numbers (r > 0), is an error: the caller never falls back to the
// default telescope.

#include "dlamb/programs/xyzls.hpp"

#include "sky_args.hpp"

#include <cstddef>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dlamb::cli {

// Telescope records of `in`; `option` names the option in error messages.
inline std::vector<programs::TelposCm> read_telpos(std::istream &in, const std::string &option,
                                                   bool corsika) {
  std::vector<programs::TelposCm> list;
  std::string line;
  std::size_t number = 0;
  while (std::getline(in, line)) {
    ++number;
    const auto hash = line.find('#');
    if (hash != std::string::npos)
      line.resize(hash);
    std::istringstream fields(line);
    std::vector<std::string> words;
    for (std::string word; fields >> word;)
      words.push_back(word);
    if (words.empty())
      continue;
    if (corsika) {
      if (words[0] != "TELESCOPE")
        continue;
      words.erase(words.begin());
    }
    const std::string where = option + " line " + std::to_string(number);
    if (words.size() != 4)
      throw std::invalid_argument(where + ": needs four values x y z r, got " +
                                  std::to_string(words.size()));
    programs::TelposCm t;
    t.x = to_double(words[0], where);
    t.y = to_double(words[1], where);
    t.z = to_double(words[2], where);
    t.r = to_double(words[3], where);
    if (!(t.r > 0.0))
      throw std::invalid_argument(where + ": telescope radius must be positive");
    list.push_back(t);
  }
  if (list.empty())
    throw std::invalid_argument(option + ": no telescope records" +
                                (corsika ? " (TELESCOPE lines)" : "") + " in the file");
  return list;
}

inline std::vector<programs::TelposCm> read_telpos(const std::string &path,
                                                   const std::string &option, bool corsika) {
  std::ifstream in(path);
  if (!in)
    throw std::invalid_argument(option + ": cannot open telescope file " + path);
  return read_telpos(in, option, corsika);
}

} // namespace dlamb::cli
