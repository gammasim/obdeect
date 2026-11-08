#pragma once

// Command-line helpers shared by the sky-level tools: strict number parsing and the
// spectrum, pulse and angular syntax (case-insensitive
// names, parameters after a colon, optional ",delay" for pulses).

#include "dlamb/density.hpp"
#include "dlamb/source.hpp"

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dlamb::cli {

inline double to_double(const std::string &text, const std::string &option) {
  std::size_t used = 0;
  double value = 0.0;
  try {
    value = std::stod(text, &used);
  } catch (const std::exception &) {
    used = 0;
  }
  if (used == 0 || used != text.size() || !std::isfinite(value))
    throw std::invalid_argument(option + ": not a number: " + text);
  return value;
}

inline std::int64_t to_int(const std::string &text, const std::string &option) {
  std::size_t used = 0;
  long long value = 0;
  try {
    value = std::stoll(text, &used);
  } catch (const std::exception &) {
    used = 0;
  }
  if (used == 0 || used != text.size())
    throw std::invalid_argument(option + ": not an integer: " + text);
  return value;
}

// Event count: an integer in 0..2^32-1 (negative values would wrap when cast).
inline std::uint32_t to_count(const std::string &text, const std::string &option) {
  const std::int64_t value = to_int(text, option);
  if (value < 0 || value > static_cast<std::int64_t>(std::numeric_limits<std::uint32_t>::max()))
    throw std::invalid_argument(option + ": must be an integer in 0..4294967295: " + text);
  return static_cast<std::uint32_t>(value);
}

inline std::vector<double> to_list(const std::string &text, const std::string &option) {
  std::vector<double> values;
  std::size_t start = 0;
  while (start <= text.size()) {
    const std::size_t comma = text.find(',', start);
    const std::string item =
        text.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
    values.push_back(to_double(item, option));
    if (comma == std::string::npos)
      break;
    start = comma + 1;
  }
  return values;
}

// Components of "x,y,z" (or "cx,cy,cz").
inline Vec3 to_vec3(const std::string &text, const std::string &option) {
  const auto v = to_list(text, option);
  if (v.size() != 3)
    throw std::invalid_argument(option + ": needs three comma-separated values");
  return {v[0], v[1], v[2]};
}

// Case-insensitive comparison of two ASCII strings.
inline bool iequals(const std::string &a, const std::string &b) {
  if (a.size() != b.size())
    return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (std::tolower(static_cast<unsigned char>(a[i])) !=
        std::tolower(static_cast<unsigned char>(b[i])))
      return false;
  return true;
}

inline bool starts_with_ci(const std::string &text, const std::string &prefix) {
  return text.size() >= prefix.size() && iequals(text.substr(0, prefix.size()), prefix);
}

// Fields separated by `sep`; empty fields are kept.
inline std::vector<std::string> split_fields(const std::string &text, char sep) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (true) {
    const std::size_t pos = text.find(sep, start);
    fields.push_back(
        text.substr(start, pos == std::string::npos ? std::string::npos : pos - start));
    if (pos == std::string::npos)
      return fields;
    start = pos + 1;
  }
}

// True when strtod consumes the whole text (so "4e2", "-5" and "1.5E-3" are numbers, not
// file names). Range and finiteness are checked by to_double.
inline bool is_number(const std::string &text) {
  if (text.empty())
    return false;
  const char *begin = text.c_str();
  char *end = nullptr;
  std::strtod(begin, &end);
  return end == begin + text.size();
}

// Two-column table (x, y); '#' starts a comment. Throws if no rows.
inline std::vector<std::pair<double, double>> read_table(const std::string &path,
                                                         const std::string &option) {
  std::ifstream in(path);
  if (!in)
    throw std::invalid_argument(option + ": cannot open table " + path);
  std::vector<std::pair<double, double>> rows;
  std::string line;
  while (std::getline(in, line)) {
    const auto hash = line.find('#');
    if (hash != std::string::npos)
      line.resize(hash);
    std::istringstream fields(line);
    double x = 0.0, y = 0.0;
    if (fields >> x >> y)
      rows.emplace_back(x, y);
  }
  if (rows.empty())
    throw std::invalid_argument(option + ": no two-column rows in " + path);
  return rows;
}

// Consecutive rows with the same abscissa (the NSB table repeats some wavelengths) are
// merged into one node with the mean value, so the abscissae increase strictly.
inline std::vector<std::pair<double, double>>
merge_equal_nodes(const std::vector<std::pair<double, double>> &rows) {
  std::vector<std::pair<double, double>> merged;
  std::vector<double> count;
  for (const auto &[x, y] : rows) {
    if (!merged.empty() && x == merged.back().first) {
      merged.back().second += y;
      count.back() += 1.0;
    } else {
      merged.emplace_back(x, y);
      count.push_back(1.0);
    }
  }
  for (std::size_t i = 0; i < merged.size(); ++i)
    merged[i].second /= count[i];
  return merged;
}

// Spectrum: "400" (nm, 200..2000), "Cherenkov:a:b" (1/lambda^2, fixed if a == b, both
// ends in 200..2000 nm with a <= b), or a table file (nm, value).
inline Spectrum parse_spectrum(const std::string &text, const std::string &option) {
  if (is_number(text)) {
    const double wl = to_double(text, option);
    if (wl < 200.0 || wl > 2000.0)
      throw std::invalid_argument(option + ": fixed wavelength must lie in 200..2000 nm");
    return Fixed{wl};
  }
  if (starts_with_ci(text, "Cherenkov:")) {
    const auto ends = split_fields(text.substr(10), ':');
    if (ends.size() != 2)
      throw std::invalid_argument(option + ": use Cherenkov:<from_nm>:<to_nm>");
    const double a = to_double(ends[0], option), b = to_double(ends[1], option);
    if (a < 200.0 || a > 2000.0 || b < 200.0 || b > 2000.0)
      throw std::invalid_argument(option + ": Cherenkov range must lie in 200..2000 nm");
    if (a > b)
      throw std::invalid_argument(option + ": Cherenkov range needs from <= to");
    if (a == b)
      return Fixed{a};
    return InverseSquare{a, b};
  }
  const auto rows = merge_equal_nodes(read_table(text, option));
  std::vector<double> x, y;
  for (const auto &[wl, value] : rows) {
    x.push_back(wl);
    y.push_back(value);
  }
  return Tabulated{LinearDensity(x, y)};
}

// Pulse: a bare number (rectangular FWHM, or Gauss sigma when `bare_is_gauss`), or
// "Simple:w[,d]" (top hat of width w), "Gauss:s[,d]" / "Gauss-rms:s[,d]" (sigma s),
// "Gauss-fwhm:f[,d]", "Exp:tau[,d]" (mean d), "Complex:wr,wg,we[,d]", or a table file
// (time ns, value).
inline Pulse parse_pulse(const std::string &text, const std::string &option, bool bare_is_gauss) {
  if (is_number(text)) {
    if (bare_is_gauss)
      return GaussianShape{0.0, to_double(text, option)};
    const double w = to_double(text, option);
    return w > 0.0 ? Pulse{UniformRange{-0.5 * w, 0.5 * w}} : Pulse{Fixed{0.0}};
  }
  const auto colon = text.find(':');
  if (colon == std::string::npos) {
    const auto rows = merge_equal_nodes(read_table(text, option));
    std::vector<double> x, y;
    for (const auto &[time, value] : rows) {
      x.push_back(time);
      y.push_back(value);
    }
    return Tabulated{LinearDensity(x, y)};
  }
  const std::string name = text.substr(0, colon);
  const auto params = to_list(text.substr(colon + 1), option);
  // Wrong parameter counts are errors, not silently ignored extras.
  const auto arity = [&](std::size_t fewest, std::size_t most) {
    if (params.size() < fewest || params.size() > most)
      throw std::invalid_argument(option + ": " + text + " takes " + std::to_string(fewest) +
                                  " to " + std::to_string(most) + " parameters");
  };
  const auto delay_param = [&]() { return params.size() > 1 ? params[1] : 0.0; };
  if (iequals(name, "Simple")) {
    arity(1, 2);
    const double w = params[0], delay = delay_param();
    return w > 0.0 ? Pulse{UniformRange{delay - 0.5 * w, delay + 0.5 * w}} : Pulse{Fixed{delay}};
  }
  if (starts_with_ci(name, "Gauss-fwhm")) {
    arity(1, 2);
    return GaussianShape{delay_param(), params[0] / 2.35482004503095};
  }
  if (starts_with_ci(name, "Gauss")) {
    arity(1, 2);
    return GaussianShape{delay_param(), params[0]};
  }
  if (starts_with_ci(name, "Exp")) {
    arity(1, 2);
    return Exponential{delay_param(), params[0]};
  }
  if (starts_with_ci(name, "Complex")) {
    arity(3, 4);
    const double d = params.size() > 3 ? params[3] : 0.0;
    return ComplexPulse{params[0], params[1], params[2], d};
  }
  throw std::invalid_argument(option + ": unknown pulse " + text);
}

// Angular profile about the pointing direction: "isotropic", "Gauss:s" / "Gauss-rms:s"
// (deg, 2-D rms) / "Gauss-fwhm:f" (deg), "Rayleigh", "Cone:a" (deg), "Parallel",
// "FilledCone:a" (deg), or a table file (angle deg, value).
inline Angular parse_angular(const std::string &text, const std::string &option) {
  const double deg = pi / 180.0;
  if (iequals(text, "isotropic"))
    return Isotropic{};
  if (iequals(text, "Parallel"))
    return FixedAngle{0.0};
  if (starts_with_ci(text, "Gauss-fwhm:"))
    return TabulatedAngular{
        gaussian_angular_table(to_double(text.substr(11), option) / 2.35482004503095 * deg)};
  if (starts_with_ci(text, "Gauss-rms:"))
    return TabulatedAngular{gaussian_angular_table(to_double(text.substr(10), option) * deg)};
  if (starts_with_ci(text, "Gauss:"))
    return TabulatedAngular{gaussian_angular_table(to_double(text.substr(6), option) * deg)};
  if (iequals(text, "Rayleigh"))
    return TabulatedAngular{rayleigh_angular_table()};
  if (starts_with_ci(text, "Cone:"))
    return FixedAngle{1.0 - std::cos(to_double(text.substr(5), option) * deg)};
  if (starts_with_ci(text, "FilledCone:"))
    return FilledCone{1.0 - std::cos(to_double(text.substr(11), option) * deg)};
  const auto rows = merge_equal_nodes(read_table(text, option));
  std::vector<double> u, value;
  for (const auto &[angle, v] : rows) {
    u.push_back(1.0 - std::cos(angle * deg));
    value.push_back(v);
  }
  return TabulatedAngular{LinearDensity(u, value)};
}

} // namespace dlamb::cli
