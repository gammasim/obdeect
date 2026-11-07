#pragma once

// Command-line helpers shared by the sky-level tools: strict number parsing and the
// spectrum, pulse and angular syntax of sim_telarray LightEmission (case-insensitive
// names, parameters after a colon, optional ",delay" for pulses).

#include "dlamb/density.hpp"
#include "dlamb/source.hpp"

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
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

inline bool starts_with_ci(const std::string &text, const std::string &prefix) {
  if (text.size() < prefix.size())
    return false;
  for (std::size_t i = 0; i < prefix.size(); ++i)
    if (std::tolower(static_cast<unsigned char>(text[i])) !=
        std::tolower(static_cast<unsigned char>(prefix[i])))
      return false;
  return true;
}

inline bool is_number(const std::string &text) {
  if (text.empty())
    return false;
  for (const char c : text)
    if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '.'))
      return false;
  return true;
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

// Spectrum: "400" (nm, 200..2000), "Cherenkov:a:b" (1/lambda^2, fixed if a == b), or a
// table file (nm, value).
inline Spectrum parse_spectrum(const std::string &text, const std::string &option) {
  if (is_number(text)) {
    const double wl = to_double(text, option);
    if (wl < 200.0 || wl > 2000.0)
      throw std::invalid_argument(option + ": fixed wavelength must lie in 200..2000 nm");
    return Fixed{wl};
  }
  if (starts_with_ci(text, "Cherenkov:")) {
    double a = 0.0, b = 0.0;
    if (std::sscanf(text.c_str() + 10, "%lf:%lf", &a, &b) != 2)
      throw std::invalid_argument(option + ": use Cherenkov:<from_nm>:<to_nm>");
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
  const double delay = params.size() > 1 ? params[1] : 0.0;
  if (name == "Simple" || name == "simple" || name == "SIMPLE") {
    const double w = params.at(0);
    return w > 0.0 ? Pulse{UniformRange{delay - 0.5 * w, delay + 0.5 * w}} : Pulse{Fixed{delay}};
  }
  if (starts_with_ci(name, "Gauss-fwhm"))
    return GaussianShape{delay, params.at(0) / 2.35482004503095};
  if (starts_with_ci(name, "Gauss"))
    return GaussianShape{delay, params.at(0)};
  if (starts_with_ci(name, "Exp"))
    return Exponential{delay, params.at(0)};
  if (starts_with_ci(name, "Complex")) {
    const double d = params.size() > 3 ? params[3] : 0.0;
    return ComplexPulse{params.at(0), params.at(1), params.at(2), d};
  }
  throw std::invalid_argument(option + ": unknown pulse " + text);
}

// Angular profile about the pointing direction: "isotropic", "Gauss:s" / "Gauss-rms:s"
// (deg, 2-D rms) / "Gauss-fwhm:f" (deg), "Rayleigh", "Cone:a" (deg), "Parallel",
// "FilledCone:a" (deg), or a table file (angle deg, value).
inline Angular parse_angular(const std::string &text, const std::string &option) {
  const double deg = pi / 180.0;
  if (starts_with_ci(text, "isotropic") && text.size() == 9)
    return Isotropic{};
  if (text == "Parallel" || text == "parallel")
    return FixedAngle{0.0};
  if (starts_with_ci(text, "Gauss-fwhm:"))
    return TabulatedAngular{
        gaussian_angular_table(to_double(text.substr(11), option) / 2.35482004503095 * deg)};
  if (starts_with_ci(text, "Gauss-rms:"))
    return TabulatedAngular{gaussian_angular_table(to_double(text.substr(10), option) * deg)};
  if (starts_with_ci(text, "Gauss:"))
    return TabulatedAngular{gaussian_angular_table(to_double(text.substr(6), option) * deg)};
  if (text == "Rayleigh" || text == "rayleigh")
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
