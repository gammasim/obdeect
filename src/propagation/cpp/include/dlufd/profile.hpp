#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <cstdlib>

namespace dlufd {

// Public numeric inputs must be finite. NaN fails every boundary comparison and converting it
// to an integer is undefined behaviour; infinities would be silently clamped.
inline void require_finite(double value, const char *name) {
  if (!std::isfinite(value))
    throw std::invalid_argument(std::string(name) + " must be finite");
}

// Natural cubic spline with the coefficients of CORSIKA rpolator.c, set_1d_cubic_params()
// with clamped = 0. On interval i, y = a + b r + c r^2 + d r^3 with r = x - x[i].
class CubicSpline {
public:
  CubicSpline(std::vector<double> x, std::span<const double> y) : x_(std::move(x)) {
    const std::size_t n = x_.size();
    if (n < 4 || y.size() != n)
      throw std::invalid_argument("cubic spline needs at least 4 points and matching values");
    for (std::size_t i = 0; i + 1 < n; ++i)
      if (!(x_[i + 1] > x_[i]))
        throw std::invalid_argument("cubic spline abscissae must be strictly ascending");

    // rpolator.c set_1d_cubic_params(), natural boundary conditions.
    std::vector<double> h(n - 1), a(y.begin(), y.end()), big_a(n, 0.0), l(n), m(n - 1), z(n), c(n),
        b(n - 1), d(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i)
      h[i] = x_[i + 1] - x_[i];
    for (std::size_t i = 1; i + 1 < n; ++i)
      big_a[i] = 3. * (a[i + 1] - a[i]) / h[i] - 3. * (a[i] - a[i - 1]) / h[i - 1];
    l[0] = 1.;
    m[0] = 0.;
    z[0] = 0.;
    for (std::size_t i = 1; i + 1 < n; ++i) {
      l[i] = 2. * (x_[i + 1] - x_[i - 1]) - h[i - 1] * m[i - 1];
      m[i] = h[i] / l[i];
      z[i] = (big_a[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    l[n - 1] = 1.;
    z[n - 1] = 0.;
    c[n - 1] = 0.;
    for (std::size_t j = n - 1; j-- > 0;) {
      c[j] = z[j] - m[j] * c[j + 1];
      b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2. * c[j]) / 3.;
      d[j] = (c[j + 1] - c[j]) / (3. * h[j]);
    }
    coefficients_.resize(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i)
      coefficients_[i] = {a[i], b[i], c[i], d[i]};
    first_value_ = a.front();
    last_value_ = a.back();
  }

  // rpol_cspline(..., clip = 0): the edge values outside [x[0], x[n-1]].
  double value(double xp) const {
    if (xp < x_.front())
      return first_value_;
    if (xp > x_.back())
      return last_value_;
    return polynomial(xp);
  }

  // The cubic of the end interval also outside the range (obdeect ambient_group_index).
  double polynomial(double xp) const {
    const std::size_t i = interval(xp);
    const double r = xp - x_[i];
    const auto &k = coefficients_[i];
    return ((k[3] * r + k[2]) * r + k[1]) * r + k[0];
  }

private:
  std::size_t interval(double xp) const {
    const auto upper = std::upper_bound(x_.begin(), x_.end(), xp);
    const std::ptrdiff_t index = (upper - x_.begin()) - 1;
    const auto last = static_cast<std::ptrdiff_t>(x_.size()) - 2;
    return static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(index, 0, last));
  }

  std::vector<double> x_;
  std::vector<std::array<double, 4>> coefficients_;
  double first_value_ = 0., last_value_ = 0.;
};

// CORSIKA atmo.c atmospheric table and its interpolation.
//
// Input: CORSIKA / sim_telarray atmprof text. Lines starting with '#' and blank lines are
// skipped. The first four columns are used: altitude [km a.s.l.], density [g/cm^3],
// thickness [g/cm^2], n-1; further columns are ignored (as in the CORSIKA short format).
//
// Public heights are in m a.s.l. The *_cm methods use the units of atmo.c (cm).
// Interpolation: density, thickness and n-1 are linear in the 40000-point equidistant table
// that atmo.c derives from cubic splines of the input levels (init_fast_interpolation).
class AtmosphereProfile {
public:
  // atmo.c MAX_PROFILE: most table levels accepted.
  static constexpr std::size_t max_levels = 120;
  // atmo.c MAX_FAST_PROFILE: equidistant lookup points.
  static constexpr std::size_t fast_points = 40000;
  // atmo.c etadsn default (n-1 per density) used while the fast table of n-1 is built.
  static constexpr double default_etadsn = 0.000283 * 994186.38 / 1222.656;

  explicit AtmosphereProfile(std::string_view text) {
    std::vector<double> alt_km, rho, thick, n_1;
    parse(text, alt_km, rho, thick, n_1);
    build(alt_km, rho, thick, n_1);
  }

  static AtmosphereProfile from_file(const std::string &path) {
    std::ifstream file(path);
    if (!file)
      throw std::invalid_argument("cannot open atmospheric profile: " + path);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return AtmosphereProfile(buffer.str());
  }

  std::size_t levels() const { return alt_cm_.size(); }
  double bottom_m() const { return bottom_cm_ / 100.; }
  double top_m() const { return top_cm_ / 100.; }
  // Table columns as read, in m a.s.l. and dimensionless n-1.
  std::span<const double> table_altitude_m() const { return table_alt_m_; }
  std::span<const double> table_index_minus_one() const { return n_1_; }

  // Public height-in-metre lookups (same interpolation as the _cm methods).
  // Throw std::invalid_argument for non-finite heights (also after the conversion to cm).
  double density(double height_m) const { return density_cm(100. * height_m); }
  double thickness(double height_m) const { return thickness_cm(100. * height_m); }
  double refractive_index_minus_one(double height_m) const {
    return index_minus_one_cm(100. * height_m);
  }

  // atmo.c rhofx_: density [g/cm^3] at height [cm].
  double density_cm(double h) const {
    require_finite(h, "height");
    if (h <= bottom_cm_)
      return rho_.front();
    if (h > top_cm_ + 0.001)
      return 0.;
    const auto i = index_of(h);
    if (i >= fast_points - 1)
      return fast_rho_.back();
    const double r = fast_h_fac_ * (h - fast_alt_[i]);
    return r * (fast_rho_[i + 1] - fast_rho_[i]) + fast_rho_[i];
  }

  // atmo.c thickx_: thickness [g/cm^2] at height [cm], with the exponential top layer.
  double thickness_cm(double h) const {
    require_finite(h, "height");
    if (h <= bottom_cm_)
      return thick_.front();
    if (h >= top_cm_)
      return 0.;
    if (h > alt_cm_[alt_cm_.size() - 2])
      return top_layer_hscale_rho0_cfac_ * (std::exp(-h / top_layer_hscale_) - top_layer_exp_top_);
    const auto i = index_of(h);
    if (i >= fast_points - 1)
      return 0.;
    const double r = fast_h_fac_ * (h - fast_alt_[i]);
    return r * (fast_thick_[i + 1] - fast_thick_[i]) + fast_thick_[i];
  }

  // atmo.c refim1x_: n-1 at height [cm].
  double index_minus_one_cm(double h) const {
    require_finite(h, "height");
    if (h <= bottom_cm_)
      return n_1_.front();
    if (h >= top_cm_)
      return 0.;
    const auto i = index_of(h);
    if (i >= fast_points - 1)
      return 0.;
    const double r = fast_h_fac_ * (h - fast_alt_[i]);
    return r * (fast_n1_[i + 1] - fast_n1_[i]) + fast_n1_[i];
  }

  // Integral of n-1 over height [cm] from the bottom of the table to h (exact for the
  // piecewise-linear table). Heights above the top count to the top; below the bottom throw.
  double index_integral_cm(double h) const {
    require_finite(h, "height");
    if (h < bottom_cm_)
      throw std::invalid_argument("height below the atmospheric table");
    if (h >= top_cm_)
      return index_cumulative_cm_.back();
    const auto i = index_of(h);
    const double partial = 0.5 * (fast_n1_[i] + index_minus_one_cm(h)) * (h - fast_alt_[i]);
    return index_cumulative_cm_[i] + partial;
  }

  double bottom_cm() const { return bottom_cm_; }
  double top_cm() const { return top_cm_; }
  double fast_height_factor() const { return fast_h_fac_; }
  std::span<const double> fast_altitude_cm() const { return fast_alt_; }
  // Raw table columns in atmo.c units.

private:
  static void parse(std::string_view text, std::vector<double> &alt_km, std::vector<double> &rho,
                    std::vector<double> &thick, std::vector<double> &n_1) {
    std::size_t line_start = 0, line_number = 0;
    while (line_start <= text.size()) {
      std::size_t line_end = text.find('\n', line_start);
      if (line_end == std::string_view::npos)
        line_end = text.size();
      const std::string_view line = text.substr(line_start, line_end - line_start);
      line_start = line_end + 1;
      ++line_number;
      const auto first = line.find_first_not_of(" \t\r");
      if (first == std::string_view::npos || line[first] == '#')
        continue;
      std::array<double, 4> row{};
      std::size_t count = 0, position = first;
      while (count < row.size() && position < line.size()) {
        const auto token_end = std::min(line.find_first_of(" \t\r", position), line.size());
        const std::string token(line.substr(position, token_end - position));
        char *end = nullptr;
        row[count] = std::strtod(token.c_str(), &end);
        if (end != token.c_str() + token.size() || token.empty())
          throw std::invalid_argument("atmospheric profile line " + std::to_string(line_number) +
                                      ": not a number: " + token);
        ++count;
        position = line.find_first_not_of(" \t\r", token_end);
        if (position == std::string_view::npos)
          break;
      }
      if (count < row.size())
        throw std::invalid_argument("atmospheric profile line " + std::to_string(line_number) +
                                    ": expected 4 columns");
      alt_km.push_back(row[0]);
      rho.push_back(row[1]);
      thick.push_back(row[2]);
      n_1.push_back(row[3]);
      if (line_end >= text.size())
        break;
    }
  }

  // Table validation and the atmo.c set-up: init_common_atmosphere(), init_fast_interpolation().
  void build(const std::vector<double> &alt_km, const std::vector<double> &rho,
             const std::vector<double> &thick, const std::vector<double> &n_1) {
    const std::size_t n = alt_km.size();
    if (n < 5 || n > max_levels)
      throw std::invalid_argument("atmospheric profile needs 5 to 120 levels");
    alt_cm_.resize(n);
    rho_ = rho;
    thick_ = thick;
    n_1_ = n_1;
    table_alt_m_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(alt_km[i]) || !std::isfinite(rho[i]) || !std::isfinite(thick[i]) ||
          !std::isfinite(n_1[i]))
        throw std::invalid_argument("atmospheric profile contains non-finite values");
      if (!(rho[i] > 0.) || !(n_1[i] > 0.) || thick[i] < 0.)
        throw std::invalid_argument("atmospheric profile needs positive density and n-1");
      alt_cm_[i] = alt_km[i] * 1e5;
      table_alt_m_[i] = alt_km[i] * 1e3;
      if (i > 0 && !(alt_cm_[i] > alt_cm_[i - 1]))
        throw std::invalid_argument("atmospheric altitudes must be strictly ascending");
      if (i > 0 && !(rho[i] < rho[i - 1]))
        throw std::invalid_argument("atmospheric density must decrease strictly with height");
      if (i > 0 && !(thick[i] < thick[i - 1]))
        throw std::invalid_argument("atmospheric thickness must decrease strictly with height");
    }
    if (!(rho.back() > 0.))
      throw std::invalid_argument("top density must be positive");
    bottom_cm_ = alt_cm_.front();
    top_cm_ = alt_cm_.back();

    // Density and thickness columns must agree (init_common_atmosphere).
    double integrated = thick.back();
    for (std::size_t il = n - 1; il > 0; --il) {
      const double hs =
          (alt_cm_[il] - alt_cm_[il - 1]) / (std::log(rho[il - 1]) - std::log(rho[il]));
      const double rho0 = rho[il] / std::exp(-alt_cm_[il] / hs);
      integrated += hs * rho0 * (std::exp(-alt_cm_[il - 1] / hs) - std::exp(-alt_cm_[il] / hs));
    }
    if (std::fabs(std::log(integrated / thick[0])) > 0.05)
      throw std::invalid_argument("atmospheric density and thickness columns are inconsistent");

    // Exponential top layer (init_common_atmosphere).
    const double h2 = alt_cm_[n - 2], h1 = alt_cm_[n - 1];
    const double rho2 = rho[n - 2], rho1 = rho[n - 1];
    const double hscale = (h1 - h2) / std::log(rho2 / rho1);
    const double rho0 = rho2 / std::exp(-h2 / hscale);
    const double t2 = hscale * rho0 * (std::exp(-h2 / hscale) - std::exp(-h1 / hscale));
    const double cfac = thick[n - 2] / t2;
    top_layer_hscale_ = hscale;
    top_layer_hscale_rho0_cfac_ = rho0 * cfac * hscale;
    top_layer_exp_top_ = std::exp(-top_cm_ / hscale);

    // Cubic splines of the levels; log(thickness) excludes the top level, which has zero
    // thickness, and is extended by three points in the top layer (init_fast_interpolation).
    std::vector<double> log_rho(n), log_n1(n), log_thick(n);
    for (std::size_t i = 0; i < n; ++i) {
      log_rho[i] = std::log(rho[i]);
      log_n1[i] = std::log(n_1[i]);
      log_thick[i] = thick[i] > 0. ? std::log(thick[i]) : -1000.;
    }
    const CubicSpline spline_log_rho(alt_cm_, log_rho);
    const CubicSpline spline_log_n1(alt_cm_, log_n1);
    const double h12 = 0.5 * h1 + 0.5 * h2;
    const double h14 = 0.25 * h1 + 0.75 * h2;
    const double h18 = 0.125 * h1 + 0.875 * h2;
    const double t12 = hscale * cfac * rho0 * (std::exp(-h12 / hscale) - std::exp(-h1 / hscale));
    const double t14 = hscale * cfac * rho0 * (std::exp(-h14 / hscale) - std::exp(-h1 / hscale));
    const double t18 = hscale * cfac * rho0 * (std::exp(-h18 / hscale) - std::exp(-h1 / hscale));
    std::vector<double> ext_alt(alt_cm_.begin(), alt_cm_.end() - 1),
        ext_log_thick(log_thick.begin(), log_thick.end() - 1);
    ext_alt.insert(ext_alt.end(), {h18, h14, h12});
    ext_log_thick.insert(ext_log_thick.end(), {std::log(t18), std::log(t14), std::log(t12)});
    const CubicSpline spline_log_thick(ext_alt, ext_log_thick);

    // Equidistant table; in the top layer the exponential fit is used (linear in height).
    fast_h_fac_ = static_cast<double>(fast_points - 1) / (top_cm_ - bottom_cm_);
    std::vector<double> fast_log_rho(fast_points), fast_log_thick(fast_points),
        fast_log_n1(fast_points);
    fast_alt_.resize(fast_points);
    fast_rho_.resize(fast_points);
    fast_thick_.resize(fast_points);
    fast_n1_.resize(fast_points);
    for (std::size_t i = 0; i < fast_points; ++i) {
      const double h =
          i < fast_points - 1 ? bottom_cm_ + static_cast<double>(i) / fast_h_fac_ : top_cm_;
      fast_alt_[i] = h;
      if (n >= 2 && h > alt_cm_[n - 2]) {
        const double rho_top = cfac * rho0 * std::exp(-h / hscale);
        const double t_top =
            hscale * cfac * rho0 * (std::exp(-h / hscale) - std::exp(-top_cm_ / hscale));
        // Extrapolate the logarithms where the top-layer value is zero (the top itself).
        fast_log_rho[i] = rho_top > 0.
                              ? std::log(rho_top)
                              : fast_log_rho[i - 1] + (fast_log_rho[i - 1] - fast_log_rho[i - 2]);
        fast_log_thick[i] =
            t_top > 0. ? std::log(t_top)
                       : fast_log_thick[i - 1] + (fast_log_thick[i - 1] - fast_log_thick[i - 2]);
        // atmo.c uses etadsn here; it still holds its default value at this point.
        fast_log_n1[i] = rho_top > 0.
                             ? std::log(default_etadsn * rho_top)
                             : fast_log_n1[i - 1] + (fast_log_n1[i - 1] - fast_log_n1[i - 2]);
      } else {
        fast_log_rho[i] = spline_log_rho.value(h);
        fast_log_n1[i] = spline_log_n1.value(h);
        fast_log_thick[i] = spline_log_thick.value(h);
      }
      fast_rho_[i] = std::exp(fast_log_rho[i]);
      fast_thick_[i] = std::exp(fast_log_thick[i]);
      fast_n1_[i] = std::exp(fast_log_n1[i]);
      if (i == fast_points - 1)
        fast_thick_[i] = fast_thick_[i - 1] * (fast_thick_[i - 1] / fast_thick_[i - 2]);
      if (i >= 1 && !(fast_thick_[i] < fast_thick_[i - 1]))
        throw std::invalid_argument(
            "atmospheric thickness is not strictly decreasing after interpolation");
    }

    index_cumulative_cm_.assign(fast_points, 0.);
    for (std::size_t i = 1; i < fast_points; ++i)
      index_cumulative_cm_[i] = index_cumulative_cm_[i - 1] + 0.5 *
                                                                  (fast_n1_[i - 1] + fast_n1_[i]) *
                                                                  (fast_alt_[i] - fast_alt_[i - 1]);
  }

  // Callers have checked that h is finite; the conversion is then well defined.
  std::size_t index_of(double h) const {
    const auto i = static_cast<long long>(fast_h_fac_ * (h - bottom_cm_));
    return static_cast<std::size_t>(
        std::clamp<long long>(i, 0, static_cast<long long>(fast_points) - 1));
  }

  friend class RefractionTables;

  std::vector<double> table_alt_m_, rho_, thick_, n_1_, alt_cm_;
  double bottom_cm_ = 0., top_cm_ = 0.;
  double top_layer_hscale_ = 0., top_layer_hscale_rho0_cfac_ = 0., top_layer_exp_top_ = 0.;
  double fast_h_fac_ = 0.;
  std::vector<double> fast_alt_, fast_rho_, fast_thick_, fast_n1_, index_cumulative_cm_;
};

} // namespace dlufd
