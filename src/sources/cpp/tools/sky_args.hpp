#pragma once

// Spectrum, pulse and angular option texts of the sky tools, including table files
// (docs/spec/cli-conventions.md). Namespace dlamb::cli.

#include "dlamb/cli.hpp"
#include "dlamb/density.hpp"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dlamb::cli {

// Two-column table file: '#' starts a comment, lines that are not two numbers are skipped,
// consecutive equal abscissae merge into one node with the mean value.
inline std::pair<std::vector<double>, std::vector<double>> read_table(const std::string &path,
                                                                      const std::string &option) {
  std::ifstream in(path);
  if (!in)
    throw std::invalid_argument(option + ": cannot open table " + path);
  std::vector<double> x, y;
  std::vector<int> repeats;
  std::string line;
  while (std::getline(in, line)) {
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos)
      line.erase(hash);
    std::istringstream words(line);
    std::string a, b;
    if (!(words >> a >> b))
      continue;
    char *end_a = nullptr, *end_b = nullptr;
    const double va = std::strtod(a.c_str(), &end_a);
    const double vb = std::strtod(b.c_str(), &end_b);
    if (*end_a != '\0' || *end_b != '\0')
      continue;
    if (!x.empty() && x.back() == va) {
      y.back() += vb;
      ++repeats.back();
    } else {
      x.push_back(va);
      y.push_back(vb);
      repeats.push_back(1);
    }
  }
  if (x.empty())
    throw std::invalid_argument(option + ": no two-column rows in " + path);
  for (std::size_t i = 0; i < y.size(); ++i)
    y[i] /= repeats[i];
  return {x, y};
}

inline double wavelength_nm(double value, const std::string &option, const std::string &text) {
  if (!(value >= 200.0 && value <= 2000.0))
    throw std::invalid_argument(option + ": wavelength must lie in 200..2000 nm: " + text);
  return value;
}

// A wavelength in nm, "Cherenkov:a:b", or a table file (wavelength nm, value).
inline Spectrum parse_spectrum(const std::string &text, const std::string &option) {
  if (is_number(text))
    return Fixed{wavelength_nm(parse_number(text, option), option, text)};
  if (starts_with_nocase(text, "cherenkov:")) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    const std::string rest = text.substr(10);
    while (true) {
      const std::size_t colon = rest.find(':', start);
      fields.push_back(rest.substr(start, colon == std::string::npos ? colon : colon - start));
      if (colon == std::string::npos)
        break;
      start = colon + 1;
    }
    if (fields.size() != 2)
      throw std::invalid_argument(option + ": Cherenkov needs two fields a:b: " + text);
    const double a = wavelength_nm(parse_number(fields[0], option), option, text);
    const double b = wavelength_nm(parse_number(fields[1], option), option, text);
    if (a > b)
      throw std::invalid_argument(option + ": Cherenkov range must have a <= b: " + text);
    if (a == b)
      return Fixed{a};
    return InverseSquare{a, b};
  }
  auto [x, y] = read_table(text, option);
  return Tabulated{LinearDensity(std::move(x), std::move(y))};
}

// Pulse text: Simple, Gauss, Gauss-rms, Gauss-fwhm, Exp, Complex, a number, or a table file.
inline Pulse parse_pulse(const std::string &text, const std::string &option, bool bare_is_gauss) {
  if (is_number(text)) {
    const double w = parse_number(text, option);
    if (bare_is_gauss)
      return GaussianShape{0.0, w};
    return w > 0.0 ? Pulse{UniformRange{-0.5 * w, 0.5 * w}} : Pulse{Fixed{0.0}};
  }
  const std::size_t colon = text.find(':');
  if (colon == std::string::npos) {
    auto [x, y] = read_table(text, option);
    return Tabulated{LinearDensity(std::move(x), std::move(y))};
  }
  const std::string name = lower_case(text.substr(0, colon));
  const auto values = parse_list(text.substr(colon + 1), option);
  const auto arity = [&](std::size_t low, std::size_t high) {
    if (values.size() < low || values.size() > high)
      throw std::invalid_argument(option + ": wrong number of parameters: " + text);
  };
  const auto delay = [&](std::size_t index) { return values.size() > index ? values[index] : 0.0; };
  if (name == "simple") {
    arity(1, 2);
    const double w = values[0], d = delay(1);
    return w > 0.0 ? Pulse{UniformRange{d - 0.5 * w, d + 0.5 * w}} : Pulse{Fixed{d}};
  }
  if (name == "gauss" || name == "gauss-rms") {
    arity(1, 2);
    return GaussianShape{delay(1), values[0]};
  }
  if (name == "gauss-fwhm") {
    arity(1, 2);
    return GaussianShape{delay(1), values[0] / fwhm_per_sigma};
  }
  if (name == "exp") {
    arity(1, 2);
    return Exponential{delay(1), values[0]};
  }
  if (name == "complex") {
    arity(3, 4);
    return ComplexPulse{values[0], values[1], values[2], delay(3)};
  }
  throw std::invalid_argument(option + ": unknown pulse " + text);
}

inline double angle_degrees(const std::string &text, const std::string &option) {
  const double a = parse_number(text, option);
  if (!(a >= 0.0 && a <= 180.0))
    throw std::invalid_argument(option + ": angle must be finite and in [0, 180] deg: " + text);
  return a;
}

// Angular profile text; table files hold angle [deg] and value.
inline Angular parse_angular(const std::string &text, const std::string &option) {
  const std::string lower = lower_case(text);
  if (lower == "isotropic")
    return Isotropic{};
  if (lower == "parallel")
    return FixedAngle{0.0};
  if (lower == "rayleigh")
    return TabulatedAngular{rayleigh_angular_table()};
  const std::size_t colon = text.find(':');
  if (colon != std::string::npos) {
    const std::string name = lower_case(text.substr(0, colon));
    const std::string value = text.substr(colon + 1);
    if (name == "cone")
      return FixedAngle{1.0 - std::cos(angle_degrees(value, option) * pi / 180.0)};
    if (name == "filledcone") {
      const double a = angle_degrees(value, option);
      if (!(a > 0.0))
        throw std::invalid_argument(option + ": FilledCone angle must be > 0: " + text);
      return FilledCone{1.0 - std::cos(a * pi / 180.0)};
    }
    if (name == "gauss" || name == "gauss-rms")
      return TabulatedAngular{gaussian_angular_table(parse_number(value, option) * pi / 180.0)};
    if (name == "gauss-fwhm")
      return TabulatedAngular{
          gaussian_angular_table(parse_number(value, option) / fwhm_per_sigma * pi / 180.0)};
  }
  auto [angle, intensity] = read_table(text, option);
  for (double &a : angle) {
    if (!(a >= 0.0 && a <= 180.0))
      throw std::invalid_argument(option + ": table angles must lie in [0, 180] deg");
    a = 1.0 - std::cos(a * pi / 180.0);
  }
  return TabulatedAngular{LinearDensity(std::move(angle), std::move(intensity))};
}

} // namespace dlamb::cli
