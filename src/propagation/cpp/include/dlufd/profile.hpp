#pragma once

// Atmospheric profile table (atmprof format) and its interpolation. See docs/refraction.md.
//
// Columns used: altitude [km], density [g/cm^3], thickness [g/cm^2], n-1. Internally heights
// are cm. Density and n-1 are interpolated with a natural cubic spline of their logarithm
// through the levels and tabulated on a fine grid (cells of at most 3 m, every level is a grid
// point); inside a cell they are exponential, so the integrals of n-1 over height are exact
// for the tabulated profile. The thickness between two levels follows the mass of the
// interpolated density: t(h) = t_upper + (t_lower - t_upper) * (mass above h in the interval)
// / (mass of the interval).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dlufd {

namespace detail {

// (1 - exp(-y)) / y, finite at y = 0.
inline double phi(double y) noexcept {
  return std::fabs(y) < 1e-8 ? 1.0 - 0.5 * y : -std::expm1(-y) / y;
}

inline double parse_profile_number(const std::string &token, std::size_t line) {
  char *end = nullptr;
  const double value = std::strtod(token.c_str(), &end);
  if (token.empty() || end != token.c_str() + token.size() || !std::isfinite(value)) {
    throw std::invalid_argument("atmospheric profile line " + std::to_string(line) + ": '" + token +
                                "' is not a number");
  }
  return value;
}

// Natural cubic spline y(x) (zero second derivative at both ends), x strictly ascending.
class NaturalSpline {
public:
  NaturalSpline(const std::vector<double> &x, const std::vector<double> &y)
      : x_(x), y_(y), second_(x.size(), 0.) {
    const std::size_t n = x.size();
    std::vector<double> upper(n, 0.), rhs(n, 0.);
    for (std::size_t i = 1; i + 1 < n; ++i) {
      const double h0 = x[i] - x[i - 1], h1 = x[i + 1] - x[i];
      const double b = 6. * ((y[i + 1] - y[i]) / h1 - (y[i] - y[i - 1]) / h0);
      const double pivot = 2. * (h0 + h1) - h0 * upper[i - 1];
      upper[i] = h1 / pivot;
      rhs[i] = (b - h0 * rhs[i - 1]) / pivot;
    }
    for (std::size_t i = n - 1; i-- > 1;) {
      second_[i] = rhs[i] - upper[i] * second_[i + 1];
    }
  }

  [[nodiscard]] double operator()(double at) const {
    const std::size_t n = x_.size();
    const auto it = std::upper_bound(x_.begin(), x_.end(), at);
    const std::size_t i = std::min(static_cast<std::size_t>(it - x_.begin()), n - 1);
    const std::size_t k = i == 0 ? 0 : i - 1;
    const std::size_t j = std::min(k, n - 2);
    const double h = x_[j + 1] - x_[j];
    const double a = (x_[j + 1] - at) / h, b = (at - x_[j]) / h;
    return a * y_[j] + b * y_[j + 1] +
           ((a * a * a - a) * second_[j] + (b * b * b - b) * second_[j + 1]) * h * h / 6.;
  }

private:
  std::vector<double> x_, y_, second_;
};

} // namespace detail

class AtmosphereProfile {
public:
  // Gladstone-Dale ratio (n-1)/density [cm^3/g] of air; the value used to build synthetic
  // profiles.
  static constexpr double default_etadsn = 0.23;
  static constexpr std::size_t min_levels = 5;
  static constexpr std::size_t max_levels = 120;
  static constexpr int max_moment = 6; // highest power of n-1 for index_moment_cm
  static constexpr double cell_cm = 300.;

  explicit AtmosphereProfile(const std::string &text) {
    std::istringstream stream(text);
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line)) {
      ++line_number;
      std::istringstream row(line);
      std::string token;
      if (!(row >> token) || token[0] == '#') {
        continue;
      }
      double value[4];
      value[0] = detail::parse_profile_number(token, line_number);
      for (int column = 1; column < 4; ++column) {
        if (!(row >> token)) {
          throw std::invalid_argument("atmospheric profile line " + std::to_string(line_number) +
                                      ": needs 4 columns (altitude, density, thickness, n-1)");
        }
        value[column] = detail::parse_profile_number(token, line_number);
      }
      z_.push_back(value[0] * 1e5);
      rho_.push_back(value[1]);
      thick_.push_back(value[2]);
      nu_.push_back(value[3]);
    }
    validate();
    build();
  }

  [[nodiscard]] static AtmosphereProfile from_file(const std::string &path) {
    std::ifstream file(path);
    if (!file) {
      throw std::invalid_argument("cannot open atmospheric profile " + path);
    }
    std::ostringstream text;
    text << file.rdbuf();
    return AtmosphereProfile(text.str());
  }

  std::size_t levels() const noexcept { return z_.size(); }
  double bottom_m() const noexcept { return z_.front() / 100.; }
  double top_m() const noexcept { return z_.back() / 100.; }
  const std::vector<double> &table_altitude_cm() const noexcept { return z_; }
  const std::vector<double> &table_density() const noexcept { return rho_; }
  const std::vector<double> &table_thickness() const noexcept { return thick_; }
  const std::vector<double> &table_index_minus_one() const noexcept { return nu_; }

  // Lookups with heights in cm. Below the table the bottom values are returned, above it 0.
  double density_cm(double h) const {
    check(h);
    if (h > zf_.back()) {
      return 0.;
    }
    const std::size_t k = cell(h);
    return rhof_[k] * std::exp(-rate_rho_[k] * (std::max(h, zf_.front()) - zf_[k]));
  }

  double index_minus_one_cm(double h) const {
    check(h);
    if (h > zf_.back()) {
      return 0.;
    }
    const std::size_t k = cell(h);
    return nuf_[k] * std::exp(-rate_nu_[k] * (std::max(h, zf_.front()) - zf_[k]));
  }

  double thickness_cm(double h) const {
    check(h);
    if (h > z_.back()) {
      return 0.;
    }
    h = std::max(h, z_.front());
    const std::size_t j = std::min(
        static_cast<std::size_t>(std::upper_bound(z_.begin(), z_.end(), h) - z_.begin()) - 1,
        z_.size() - 2);
    const double above_level = mass_above_[node_[j + 1]];
    const double fraction = (mass_above(h) - above_level) / (mass_above_[node_[j]] - above_level);
    return thick_[j + 1] + (thick_[j] - thick_[j + 1]) * std::clamp(fraction, 0., 1.);
  }

  // Integral of (n-1) from height h to the top of the table [cm]. Below the table (n-1)
  // keeps its bottom value.
  double index_integral_cm(double h) const {
    check(h);
    const double total = moment_nodes_[1].back();
    if (h >= zf_.back()) {
      return 0.;
    }
    if (h < zf_.front()) {
      return total + nuf_[0] * (zf_.front() - h);
    }
    return total - index_moment_cm(h, 1);
  }

  // Integral of (n-1)^power from the bottom of the table to h [cm]; h above the top gives the
  // value at the top.
  double index_moment_cm(double h, int power) const {
    check(h);
    if (power < 0 || power > max_moment) {
      throw std::invalid_argument("moment order must be in [0, max_moment]");
    }
    if (h <= zf_.front()) {
      return 0.;
    }
    if (h >= zf_.back()) {
      return moment_nodes_[static_cast<std::size_t>(power)].back();
    }
    const std::size_t k = cell(h);
    const double d = h - zf_[k];
    return moment_nodes_[static_cast<std::size_t>(power)][k] +
           std::pow(nuf_[k], power) * d * detail::phi(power * rate_nu_[k] * d);
  }

  double density(double height_m) const { return density_cm(to_cm(height_m)); }
  double thickness(double height_m) const { return thickness_cm(to_cm(height_m)); }
  double refractive_index_minus_one(double height_m) const {
    return index_minus_one_cm(to_cm(height_m));
  }

private:
  static void check(double h) {
    if (!std::isfinite(h)) {
      throw std::invalid_argument("profile height must be finite");
    }
  }
  static double to_cm(double height_m) {
    const double h = 100. * height_m;
    check(h);
    return h;
  }

  // Fine-grid cell [zf_k, zf_k+1] containing h (h <= top).
  std::size_t cell(double h) const {
    const std::size_t above =
        static_cast<std::size_t>(std::upper_bound(zf_.begin(), zf_.end(), h) - zf_.begin());
    return std::min(above == 0 ? std::size_t{0} : above - 1, zf_.size() - 2);
  }

  // Mass above h (g/cm^2) of the interpolated density, h inside the table.
  double mass_above(double h) const {
    const std::size_t k = cell(h);
    const double d = h - zf_[k];
    return mass_above_[k] - rhof_[k] * d * detail::phi(rate_rho_[k] * d);
  }

  void validate() const {
    const std::size_t n = z_.size();
    if (n < min_levels || n > max_levels) {
      throw std::invalid_argument("atmospheric profile needs 5 to 120 levels, found " +
                                  std::to_string(n));
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (i > 0 && !(z_[i] > z_[i - 1])) {
        throw std::invalid_argument("profile altitudes must be strictly ascending");
      }
      if (!(rho_[i] > 0.) || !(nu_[i] > 0.)) {
        throw std::invalid_argument("profile density and n-1 must be positive");
      }
      if (!(thick_[i] >= 0.)) {
        throw std::invalid_argument("profile thickness must not be negative");
      }
      if (i > 0 && !(rho_[i] < rho_[i - 1] && thick_[i] < thick_[i - 1])) {
        throw std::invalid_argument("profile density and thickness must decrease with height");
      }
    }
  }

  void build() {
    const std::size_t n = z_.size();
    std::vector<double> log_rho(n), log_nu(n);
    for (std::size_t i = 0; i < n; ++i) {
      log_rho[i] = std::log(rho_[i]);
      log_nu[i] = std::log(nu_[i]);
    }
    const detail::NaturalSpline rho_spline(z_, log_rho), nu_spline(z_, log_nu);
    node_.assign(n, 0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
      node_[i] = zf_.size();
      const auto cells = static_cast<std::size_t>(std::ceil((z_[i + 1] - z_[i]) / cell_cm));
      for (std::size_t c = 0; c < std::max<std::size_t>(cells, 1); ++c) {
        const double z =
            z_[i] + (z_[i + 1] - z_[i]) * double(c) / double(std::max<std::size_t>(cells, 1));
        zf_.push_back(z);
        rhof_.push_back(c == 0 ? rho_[i] : std::exp(rho_spline(z)));
        nuf_.push_back(c == 0 ? nu_[i] : std::exp(nu_spline(z)));
      }
    }
    node_[n - 1] = zf_.size();
    zf_.push_back(z_.back());
    rhof_.push_back(rho_.back());
    nuf_.push_back(nu_.back());

    const std::size_t m = zf_.size();
    rate_rho_.resize(m - 1);
    rate_nu_.resize(m - 1);
    for (std::size_t k = 0; k + 1 < m; ++k) {
      const double dz = zf_[k + 1] - zf_[k];
      rate_rho_[k] = std::log(rhof_[k] / rhof_[k + 1]) / dz;
      rate_nu_[k] = std::log(nuf_[k] / nuf_[k + 1]) / dz;
    }
    mass_above_.assign(m, 0.);
    for (std::size_t k = m - 1; k-- > 0;) {
      const double dz = zf_[k + 1] - zf_[k];
      mass_above_[k] = mass_above_[k + 1] + rhof_[k] * dz * detail::phi(rate_rho_[k] * dz);
    }
    for (int power = 0; power <= max_moment; ++power) {
      auto &nodes = moment_nodes_[static_cast<std::size_t>(power)];
      nodes.assign(m, 0.);
      for (std::size_t k = 0; k + 1 < m; ++k) {
        const double dz = zf_[k + 1] - zf_[k];
        nodes[k + 1] =
            nodes[k] + std::pow(nuf_[k], power) * dz * detail::phi(power * rate_nu_[k] * dz);
      }
    }
    // The thickness column must agree with the integrated density to 5 %.
    const double column = thick_.front() - thick_.back();
    if (std::fabs(mass_above_.front() - column) > 0.05 * column) {
      throw std::invalid_argument(
          "profile thickness and density columns disagree by more than 5 %");
    }
  }

  std::vector<double> z_, rho_, thick_, nu_;            // levels
  std::vector<double> zf_, rhof_, nuf_;                 // fine grid
  std::vector<std::size_t> node_;                       // fine-grid index of each level
  std::vector<double> rate_rho_, rate_nu_, mass_above_; // per fine cell / node
  std::vector<double> moment_nodes_[max_moment + 1];
};

} // namespace dlufd
