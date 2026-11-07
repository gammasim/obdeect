#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace dlamb::cli {

// Parses a finite number for option `option`; the whole text must be a number.
[[nodiscard]] inline double parse_number(const std::string &option, const std::string &text) {
  std::size_t used = 0;
  double value = 0.0;
  try {
    value = std::stod(text, &used);
  } catch (const std::exception &) {
    used = 0;
  }
  if (used == 0 || used != text.size() || !std::isfinite(value))
    throw std::invalid_argument("invalid number for " + option + ": '" + text + "'");
  return value;
}

// Parses an integer for option `option`; the whole text must be an integer.
[[nodiscard]] inline long long parse_integer(const std::string &option, const std::string &text) {
  std::size_t used = 0;
  long long value = 0;
  try {
    value = std::stoll(text, &used);
  } catch (const std::exception &) {
    used = 0;
  }
  if (used == 0 || used != text.size())
    throw std::invalid_argument("invalid integer for " + option + ": '" + text + "'");
  return value;
}

} // namespace dlamb::cli
