#pragma once

// Telescope list files of dlamb-point-source (docs/spec/cli-conventions.md).

#include "obdeect/sources/cli.hpp"
#include "obdeect/sources/programs/point-source.hpp"

#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::sources::cli {

// Reads lines "x y z r" [cm]; with `corsika` only lines "TELESCOPE x y z r".
inline std::vector<programs::TelposCm> read_telpos(std::istream &in, const std::string &option,
                                                   bool corsika) {
  std::vector<programs::TelposCm> list;
  std::string line;
  int number = 0;
  while (std::getline(in, line)) {
    ++number;
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos)
      line.erase(hash);
    std::istringstream stream(line);
    std::vector<std::string> words;
    std::string word;
    while (stream >> word)
      words.push_back(word);
    if (words.empty())
      continue;
    if (corsika) {
      if (words[0] != "TELESCOPE")
        continue;
      words.erase(words.begin());
    }
    const std::string where = option + " line " + std::to_string(number) + ": ";
    if (words.size() != 4)
      throw std::invalid_argument(where + "needs four values x y z r, got " +
                                  std::to_string(words.size()));
    double v[4];
    for (int i = 0; i < 4; ++i) {
      try {
        v[i] = parse_number(words[i], option);
      } catch (const std::invalid_argument &) {
        throw std::invalid_argument(where + "not a number: " + words[i]);
      }
    }
    if (!(v[3] > 0.0))
      throw std::invalid_argument(where + "the radius must be > 0");
    list.push_back({v[0], v[1], v[2], v[3]});
  }
  if (list.empty())
    throw std::invalid_argument(option + ": no telescope records " +
                                (corsika ? "(TELESCOPE lines) " : "") + "in the file");
  return list;
}

} // namespace obdeect::sources::cli
