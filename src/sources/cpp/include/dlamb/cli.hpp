#pragma once

// Option parsing shared by the sky tools (docs/spec/cli-conventions.md). Errors are
// std::invalid_argument with the option name first.

#include "dlamb/vec3.hpp"

#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace dlamb::cli {

// True when strtod consumes the whole text (so "4e2" and "-5" are numbers, "4e2x" is not).
inline bool is_number(const std::string &text) {
  if (text.empty())
    return false;
  char *end = nullptr;
  (void)std::strtod(text.c_str(), &end);
  return end != text.c_str() && *end == '\0';
}

// A finite real number; the whole text must be a number.
inline double parse_number(const std::string &text, const std::string &option) {
  if (!is_number(text))
    throw std::invalid_argument(option + ": not a number: " + text);
  errno = 0;
  const double value = std::strtod(text.c_str(), nullptr);
  if (errno == ERANGE || !std::isfinite(value))
    throw std::invalid_argument(option + ": number out of range: " + text);
  return value;
}

// A signed decimal integer (digits only).
inline long long parse_integer(const std::string &text, const std::string &option) {
  std::size_t i = (!text.empty() && (text[0] == '+' || text[0] == '-')) ? 1 : 0;
  bool digits = i < text.size();
  for (std::size_t k = i; k < text.size(); ++k)
    digits = digits && std::isdigit(static_cast<unsigned char>(text[k]));
  if (!digits)
    throw std::invalid_argument(option + ": not an integer: " + text);
  errno = 0;
  const long long value = std::strtoll(text.c_str(), nullptr, 10);
  if (errno == ERANGE)
    throw std::invalid_argument(option + ": integer out of range: " + text);
  return value;
}

inline double to_double(const std::string &text, const std::string &option) {
  return parse_number(text, option);
}

// An event count: integer in 0..4294967295.
inline std::uint32_t to_count(const std::string &text, const std::string &option) {
  bool digits = !text.empty();
  for (const char c : text)
    digits = digits && std::isdigit(static_cast<unsigned char>(c));
  errno = 0;
  const unsigned long long value = digits ? std::strtoull(text.c_str(), nullptr, 10) : 0;
  if (!digits || errno == ERANGE || value > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument(option + ": must be an integer in 0..4294967295: " + text);
  return static_cast<std::uint32_t>(value);
}

// Comma separated items (an empty item is an error).
inline std::vector<std::string> split_items(const std::string &text, const std::string &option) {
  std::vector<std::string> items;
  std::size_t start = 0;
  while (true) {
    const std::size_t comma = text.find(',', start);
    const std::string item = text.substr(start, comma == std::string::npos ? comma : comma - start);
    if (item.empty())
      throw std::invalid_argument(option + ": empty item in list: " + text);
    items.push_back(item);
    if (comma == std::string::npos)
      break;
    start = comma + 1;
  }
  return items;
}

inline std::vector<double> parse_list(const std::string &text, const std::string &option) {
  std::vector<double> values;
  for (const auto &item : split_items(text, option))
    values.push_back(parse_number(item, option));
  return values;
}

// "x,y,z" with exactly three items.
inline Vec3 parse_vector(const std::string &text, const std::string &option) {
  const auto items = split_items(text, option);
  if (items.size() != 3)
    throw std::invalid_argument(option + ": needs three comma-separated values");
  return {parse_number(items[0], option), parse_number(items[1], option),
          parse_number(items[2], option)};
}

inline std::string lower_case(std::string text) {
  for (char &c : text)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return text;
}

// True when `text` starts with `prefix`, ignoring case.
inline bool starts_with_nocase(const std::string &text, const std::string &prefix) {
  return text.size() >= prefix.size() && lower_case(text.substr(0, prefix.size())) == prefix;
}

} // namespace dlamb::cli
