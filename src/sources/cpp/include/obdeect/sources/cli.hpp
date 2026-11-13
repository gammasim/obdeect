#pragma once

// Option parsing shared by the sky tools (docs/spec/cli-conventions.md). Errors are
// UsageError with the option name first.

#include "obdeect/sources/vec3.hpp"
#include "obdeect/stream/cli.hpp"
#include "obdeect/stream/parse.hpp"

#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::sources::cli {

// True for a finite number filling the whole text ("4e2" and "-5" are numbers, "4e2x" is not).
inline bool is_number(const std::string &text) {
  double value = 0.0;
  return stream::parse::finite_number(text, value);
}

// Option parsing is stream::cli (errors are UsageError with the option name first).
using stream::cli::UsageError;

inline double parse_number(const std::string &text, const std::string &option) {
  return stream::cli::number(text, option);
}

inline long long parse_integer(const std::string &text, const std::string &option) {
  return stream::cli::integer(text, option);
}

// An event count: integer in 0..4294967295.
inline std::uint32_t to_count(const std::string &text, const std::string &option) {
  long long value = 0;
  if (!stream::parse::integer(text, value) || value < 0 ||
      value > std::numeric_limits<std::uint32_t>::max())
    throw UsageError(option + ": must be an integer in 0..4294967295: " + text);
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
      throw UsageError(option + ": empty item in list: " + text);
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
    throw UsageError(option + ": needs three comma-separated values");
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

} // namespace obdeect::sources::cli
