#pragma once

// Optical-depth table of the atmosphere (extinction-table format) and the transmission of a
// photon from its emission height to the telescope. See docs/propagation/extinction.md.

#include "obdeect/propagation/quadrature.hpp"
#include "obdeect/stream/parse.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace obdeect::propagation {

// Tables are indexed by integer nm in [0, max_lambda). See WavelengthInterpolation.
inline constexpr int max_lambda = 1000;
inline constexpr double unlimited_airmass = 1e10;
inline constexpr double altitude_tolerance_m = 25.0;
inline constexpr double max_telescope_depth = 1000.0;
inline constexpr double max_total_depth = 100.0;
inline constexpr std::size_t max_table_heights = 100;

// Throws unless the wavelength is finite and non-negative.
inline void check_wavelength(double wavelength_nm) {
  if (!std::isfinite(wavelength_nm) || wavelength_nm < 0.0) {
    throw std::invalid_argument("wavelength must be finite and non-negative (nm)");
  }
}

// How a wavelength selects the table rows.
//   linear (default): optical depths linear in lambda between the rows at floor(lambda) and
//             floor(lambda)+1; the telescope depths and the efficiency bound are interpolated the
//             same way, so every lookup is continuous in lambda.
//   floor:  the row at floor(lambda), as sim_telarray does (for comparisons with sim_telarray).
// Both clamp to the first and last row.
enum class WavelengthInterpolation { floor, linear };

// "linear" or "floor"; throws std::invalid_argument for any other text.
inline WavelengthInterpolation parse_wavelength_interpolation(std::string_view text) {
  if (text == "linear") {
    return WavelengthInterpolation::linear;
  }
  if (text == "floor") {
    return WavelengthInterpolation::floor;
  }
  throw std::invalid_argument("wavelength interpolation must be linear or floor: '" +
                              std::string(text) + "'");
}

// Floor of the wavelength for table lookups, or -1 at or above max_lambda (no transmission).
// Used by the floor mode; linear-mode lookups go through TransmissionTable.
[[nodiscard]] inline int table_wavelength(double wavelength_nm) {
  check_wavelength(wavelength_nm);
  if (wavelength_nm >= max_lambda) {
    return -1;
  }
  return static_cast<int>(wavelength_nm);
}

namespace detail {

// Text of a line before a '#'.
[[nodiscard]] inline std::string_view trim_comment(std::string_view line) {
  return line.substr(0, line.find('#'));
}

[[nodiscard]] inline std::vector<std::string_view> split_whitespace(std::string_view text) {
  std::vector<std::string_view> tokens;
  std::size_t i = 0;
  while (i < text.size()) {
    while (i < text.size() && std::isspace(static_cast<unsigned char>(text[i]))) {
      ++i;
    }
    const std::size_t start = i;
    while (i < text.size() && !std::isspace(static_cast<unsigned char>(text[i]))) {
      ++i;
    }
    if (i > start) {
      tokens.push_back(text.substr(start, i - start));
    }
  }
  return tokens;
}

// A finite number filling the whole token.
[[nodiscard]] inline double parse_number(std::string_view token, const char *what) {
  double value = 0.0;
  if (!obdeect::stream::parse::finite_number(token, value)) {
    throw std::invalid_argument(std::string(what) + " '" + std::string(token) +
                                "' is not a finite number");
  }
  return value;
}

// Linear interpolation of y(x) at `at` on ascending x, clamped to the first and last value.
[[nodiscard]] inline double rpol(const std::vector<double> &x, const double *y, double at) {
  if (at <= x.front()) {
    return y[0];
  }
  if (at >= x.back()) {
    return y[x.size() - 1];
  }
  const auto upper = std::upper_bound(x.begin(), x.end(), at);
  const std::size_t i = static_cast<std::size_t>(upper - x.begin());
  return y[i - 1] + (at - x[i - 1]) / (x[i] - x[i - 1]) * (y[i] - y[i - 1]);
}

// Curved atmosphere (slant_optical_depth) uses the Gauss-Legendre table of quadrature.hpp.
// Cosines below this use h = h_t + u^2 in the first segment (the integrand 1/cos has a
// 1/sqrt singularity at h_t when cos_t -> 0).
inline constexpr double curved_substitution_cos = 0.2;
// Tolerance for cosines of the zenith angle, which may exceed 1 by rounding.
inline constexpr double curved_cos_tolerance = 1e-9;

} // namespace detail

// Optical depth inside the telescope for every integer wavelength.
struct TelescopeDepths {
  std::array<double, max_lambda> od_tel{};   // vertical depth observation level to telescope
  std::array<double, max_lambda> od_focus{}; // vertical depth over the path inside the telescope
};

class TransmissionTable {
public:
  // Default mode: linear (see WavelengthInterpolation).
  TransmissionTable(const std::string &text, double site_altitude_m,
                    WavelengthInterpolation mode = WavelengthInterpolation::linear)
      : mode_(mode) {
    if (mode != WavelengthInterpolation::floor && mode != WavelengthInterpolation::linear) {
      throw std::invalid_argument("unknown wavelength interpolation mode");
    }
    if (!std::isfinite(site_altitude_m)) {
      throw std::invalid_argument("site altitude must be finite");
    }
    parse(text);
    if (std::fabs(1000.0 * h2_km_ - site_altitude_m) > altitude_tolerance_m) {
      throw std::invalid_argument("transmission table H2 does not match the site altitude");
    }
  }

  [[nodiscard]] static TransmissionTable
  from_file(const std::string &path, double site_altitude_m,
            WavelengthInterpolation mode = WavelengthInterpolation::linear) {
    std::ifstream file(path);
    if (!file) {
      throw std::invalid_argument("cannot open transmission table " + path);
    }
    std::ostringstream text;
    text << file.rdbuf();
    return TransmissionTable(text.str(), site_altitude_m, mode);
  }

  [[nodiscard]] double observation_level_m() const noexcept { return 1000.0 * h2_km_; }
  [[nodiscard]] int first_wavelength_nm() const noexcept { return first_; }
  [[nodiscard]] int last_wavelength_nm() const noexcept { return first_ + rows_ - 1; }
  [[nodiscard]] WavelengthInterpolation wavelength_interpolation() const noexcept { return mode_; }
  // Table heights [m a.s.l.]: the observation level and the heights H1 (ascending). The extinction
  // coefficient has a kink at each of them (chord integrals split there).
  [[nodiscard]] const std::vector<double> &node_heights_m() const noexcept { return height_m_; }

  // Vertical optical depth from the observation level to height_m; linear in log10(height)
  // between the table heights, zero at and below the observation level, constant above the
  // last table height. The wavelength selects the rows as set by wavelength_interpolation()
  // (floor: truncated to integer nm; linear: interpolated) and is clamped to the table.
  template <class T> [[nodiscard]] T optical_depth(T wavelength_nm, T height_m) const {
    const Lookup at = lookup_of(static_cast<double>(wavelength_nm));
    const double h = static_cast<double>(height_m);
    if (!std::isfinite(h)) {
      throw std::invalid_argument("height must be finite");
    }
    if (h <= 1000.0 * h2_km_) {
      return T(0);
    }
    return static_cast<T>(tau_at(at, std::log10(h)));
  }

  // Extinction coefficient |d tau / dh| [1/m] at height_m: the exact derivative of optical_depth
  // (its log10-linear interpolation, per table interval; the linear blend in wavelength blends the
  // derivatives too). Zero at and below the observation level and above the last table height,
  // where tau is constant. Used by the chord integrals of near-horizontal paths (chord.hpp).
  template <class T> [[nodiscard]] T extinction_coefficient(T wavelength_nm, T height_m) const {
    const Lookup at = lookup_of(static_cast<double>(wavelength_nm));
    const double h = static_cast<double>(height_m);
    if (!std::isfinite(h)) {
      throw std::invalid_argument("height must be finite");
    }
    if (h <= 1000.0 * h2_km_) {
      return T(0);
    }
    const double log_h = std::log10(h);
    if (log_h >= grid_.back()) {
      return T(0);
    }
    return static_cast<T>(std::fabs(slope_at(at, log_h)) / (h * std::log(10.0)));
  }

  // Transmission exp(-|tau(zem) - od_tel| * airmass - od_focus). Zero below the observation
  // level, for od_tel above 1000 and for a total depth above 100.
  template <class T>
  [[nodiscard]] T transmission2(T wavelength_nm, T zem_m, T airmass, T od_tel, T od_focus) const {
    if (!std::isfinite(static_cast<double>(zem_m))) {
      throw std::invalid_argument("height must be finite");
    }
    if (!std::isfinite(static_cast<double>(airmass)) || airmass < T(0) ||
        !std::isfinite(static_cast<double>(od_tel)) || od_tel < T(0) ||
        !std::isfinite(static_cast<double>(od_focus)) || od_focus < T(0)) {
      throw std::invalid_argument("airmass and optical depths must be finite and non-negative");
    }
    if (static_cast<double>(zem_m) < 1000.0 * h2_km_ ||
        static_cast<double>(od_tel) > max_telescope_depth) {
      (void)lookup_of(static_cast<double>(wavelength_nm)); // validates the wavelength
      return T(0);
    }
    const T depth = std::fabs(optical_depth<T>(wavelength_nm, zem_m) - od_tel) * airmass + od_focus;
    if (depth > T(max_total_depth)) {
      return T(0);
    }
    return std::exp(-depth);
  }

  template <class T> [[nodiscard]] T transmission(T wavelength_nm, T zem_m, T airmass) const {
    return transmission2<T>(wavelength_nm, zem_m, airmass, T(0), T(0));
  }

  // Slant optical depth in a spherical atmosphere of Earth radius earth_radius_m:
  //   integral from h_t to zem of (d tau/dh) / cos theta(h) dh,
  // with cos theta(h) = sqrt(1 - p^2 / (R + h)^2), p = (R + h_t) sin theta_t and cos theta_t =
  // cos_zenith. Zero if zem <= h_t. The plane-parallel limit R -> infinity is
  // |tau(zem) - tau(h_t)| / cos_zenith. Computed in double precision; see
  // docs/propagation/extinction.md.
  template <class T>
  [[nodiscard]] T slant_optical_depth(T wavelength_nm, T telescope_altitude_m,
                                      T emission_altitude_m, T cos_zenith,
                                      double earth_radius_m) const {
    return static_cast<T>(slant_depth(
        static_cast<double>(wavelength_nm), static_cast<double>(telescope_altitude_m),
        static_cast<double>(emission_altitude_m), static_cast<double>(cos_zenith), earth_radius_m));
  }

  // Slant optical depth of an upward photon (emission below the telescope, zem < h_t): the same
  // integral over the heights between zem and h_t, with cos_zenith the cosine of the photon's
  // upward direction to the local vertical at the telescope. Zero if zem >= h_t. Only for upward
  // photons: slant_optical_depth keeps its zero below the telescope for downward photons. The
  // plane-parallel limit R -> infinity is |tau(h_t) - tau(zem)| / cos_zenith.
  template <class T>
  [[nodiscard]] T slant_optical_depth_upward(T wavelength_nm, T telescope_altitude_m,
                                             T emission_altitude_m, T cos_zenith,
                                             double earth_radius_m) const {
    return static_cast<T>(slant_depth(static_cast<double>(wavelength_nm),
                                      static_cast<double>(telescope_altitude_m),
                                      static_cast<double>(emission_altitude_m),
                                      static_cast<double>(cos_zenith), earth_radius_m, true));
  }

  // Transmission exp(-(slant optical depth + od_focus)) in a spherical atmosphere (see
  // slant_optical_depth). The telescope depth od_tel is tau(h_t). Zero below the observation
  // level, for od_tel above max_telescope_depth and for a total depth above max_total_depth.
  template <class T>
  [[nodiscard]] T transmission_curved(T wavelength_nm, T zem_m, T cos_zenith,
                                      T telescope_altitude_m, T od_focus,
                                      double earth_radius_m) const {
    return transmission_curved_impl<T>(wavelength_nm, zem_m, cos_zenith, telescope_altitude_m,
                                       od_focus, earth_radius_m, false);
  }

  // transmission_curved for an upward photon (zem < h_t): the slant depth is the integral over
  // [zem, h_t] (slant_optical_depth_upward). Same zero rules as transmission_curved.
  template <class T>
  [[nodiscard]] T transmission_curved_upward(T wavelength_nm, T zem_m, T cos_zenith,
                                             T telescope_altitude_m, T od_focus,
                                             double earth_radius_m) const {
    return transmission_curved_impl<T>(wavelength_nm, zem_m, cos_zenith, telescope_altitude_m,
                                       od_focus, earth_radius_m, true);
  }

  // Body of transmission_curved and transmission_curved_upward (upward selects the slant depth).
  template <class T>
  [[nodiscard]] T transmission_curved_impl(T wavelength_nm, T zem_m, T cos_zenith,
                                           T telescope_altitude_m, T od_focus,
                                           double earth_radius_m, bool upward) const {
    const double lambda = static_cast<double>(wavelength_nm);
    const double zem = static_cast<double>(zem_m);
    const double h_t = static_cast<double>(telescope_altitude_m);
    const double focus = static_cast<double>(od_focus);
    const Lookup at =
        validate_curved(lambda, h_t, zem, static_cast<double>(cos_zenith), earth_radius_m);
    if (!std::isfinite(focus) || focus < 0.0) {
      throw std::invalid_argument("od_focus must be finite and non-negative");
    }
    if (zem < observation_level_m()) {
      return T(0);
    }
    if (optical_depth<double>(lambda, h_t) > max_telescope_depth) {
      return T(0);
    }
    const double cos_t = std::min(static_cast<double>(cos_zenith), 1.0);
    const double slant = upward ? slant_depth_upward_row(at, h_t, zem, cos_t, earth_radius_m)
                                : slant_depth_row(at, h_t, zem, cos_t, earth_radius_m);
    const double depth = slant + focus;
    if (depth > max_total_depth) {
      return T(0);
    }
    return static_cast<T>(std::exp(-depth));
  }

  // Depths inside the telescope at altitude telescope_altitude_m with an optical path
  // pathlen_m: od_tel = tau(h_T), od_focus = tau(h_T + l) - tau(h_T). If tau(400 nm, h_T)
  // exceeds 1000 the table does not belong to the site and both are zero.
  [[nodiscard]] TelescopeDepths telescope_depths(double telescope_altitude_m,
                                                 double pathlen_m) const {
    if (!std::isfinite(telescope_altitude_m) || !std::isfinite(pathlen_m) || pathlen_m < 0.0) {
      throw std::invalid_argument("telescope altitude and path length must be finite, path >= 0");
    }
    TelescopeDepths depths;
    if (optical_depth<double>(400.0, telescope_altitude_m) > max_telescope_depth) {
      return depths;
    }
    for (int wavelength = 0; wavelength < max_lambda; ++wavelength) {
      const double w = static_cast<double>(wavelength);
      const double at_telescope = optical_depth<double>(w, telescope_altitude_m);
      depths.od_tel[static_cast<std::size_t>(wavelength)] = at_telescope;
      depths.od_focus[static_cast<std::size_t>(wavelength)] =
          optical_depth<double>(w, telescope_altitude_m + pathlen_m) - at_telescope;
    }
    return depths;
  }

  struct TelescopeDepthValues {
    double od_tel{};
    double od_focus{};
  };

  // Telescope depths at a wavelength [nm] from the per-integer-nm arrays of telescope_depths():
  // the array entry in floor mode, the linear blend of the entries at floor(lambda) and
  // floor(lambda)+1 in linear mode (clamped to the table). The same mode as optical_depth, so
  // od_tel equals optical_depth(lambda, h_T) and od_focus the in-telescope part, in both modes.
  // Floor mode returns the array entries exactly.
  [[nodiscard]] TelescopeDepthValues telescope_depths_at(const TelescopeDepths &depths,
                                                         double wavelength_nm) const {
    const Lookup at = lookup_of(wavelength_nm);
    const auto lower = static_cast<std::size_t>(at.lower);
    const auto upper = static_cast<std::size_t>(at.upper);
    if (upper == lower) {
      return {depths.od_tel[lower], depths.od_focus[lower]};
    }
    return {(1.0 - at.frac) * depths.od_tel[lower] + at.frac * depths.od_tel[upper],
            (1.0 - at.frac) * depths.od_focus[lower] + at.frac * depths.od_focus[upper]};
  }

private:
  // Table rows used by a wavelength: weight (1 - frac) on the row of lower [nm] and frac on the
  // row of upper. upper == lower (frac 0) in floor mode and at the last row.
  struct Lookup {
    int lower;
    int upper;
    double frac;
  };

  // Validates the wavelength and selects its rows (absolute nm, clamped to the table).
  [[nodiscard]] Lookup lookup_of(double wavelength_nm) const {
    if (!std::isfinite(wavelength_nm)) {
      throw std::invalid_argument("wavelength must be finite");
    }
    const double clamped = std::clamp(wavelength_nm, static_cast<double>(first_),
                                      static_cast<double>(last_wavelength_nm()));
    const int lower = static_cast<int>(clamped);
    if (mode_ == WavelengthInterpolation::floor || lower >= last_wavelength_nm()) {
      return {lower, lower, 0.0};
    }
    return {lower, lower + 1, clamped - static_cast<double>(lower)};
  }

  // Optical depth at log10(height) [m] of the rows of `at`, blended by frac.
  [[nodiscard]] double tau_at(const Lookup &at, double log_height) const {
    const std::size_t n = grid_.size();
    const double *lower = &tau_[static_cast<std::size_t>(at.lower - first_) * n];
    const double a = detail::rpol(grid_, lower, log_height);
    if (at.upper == at.lower) {
      return a;
    }
    const double b = detail::rpol(grid_, lower + n, log_height);
    return (1.0 - at.frac) * a + at.frac * b;
  }

  // d tau / d log10(h) of the rows of `at` at log_height (the slope of the interval containing it,
  // blended as tau_at). Requires grid_[0] <= log_height < grid_.back().
  [[nodiscard]] double slope_at(const Lookup &at, double log_height) const {
    const std::size_t n = grid_.size();
    const auto upper = std::upper_bound(grid_.begin(), grid_.end(), log_height);
    const auto i = static_cast<std::size_t>(upper - grid_.begin()); // grid_[i-1] <= x < grid_[i]
    const double *lower = &tau_[static_cast<std::size_t>(at.lower - first_) * n];
    const auto slope = [&](const double *y) {
      return (y[i] - y[i - 1]) / (grid_[i] - grid_[i - 1]);
    };
    double value = slope(lower);
    if (at.upper != at.lower) {
      value = (1.0 - at.frac) * value + at.frac * slope(lower + n);
    }
    return value;
  }

  // Validates the arguments of the curved-atmosphere lookups and returns the table rows.
  [[nodiscard]] Lookup validate_curved(double wavelength_nm, double telescope_altitude_m,
                                       double emission_altitude_m, double cos_zenith,
                                       double earth_radius_m) const {
    const Lookup at = lookup_of(wavelength_nm);
    if (!std::isfinite(telescope_altitude_m) || !std::isfinite(emission_altitude_m)) {
      throw std::invalid_argument("altitudes must be finite");
    }
    if (!std::isfinite(cos_zenith) || cos_zenith < 0.0 ||
        cos_zenith > 1.0 + detail::curved_cos_tolerance) {
      throw std::invalid_argument("cos_zenith must be finite and within [0, 1]");
    }
    if (!std::isfinite(earth_radius_m) || !(earth_radius_m > 0.0)) {
      throw std::invalid_argument("earth radius must be finite and positive");
    }
    return at;
  }

  [[nodiscard]] double slant_depth(double wavelength_nm, double telescope_altitude_m,
                                   double emission_altitude_m, double cos_zenith,
                                   double earth_radius_m, bool upward = false) const {
    const Lookup at = validate_curved(wavelength_nm, telescope_altitude_m, emission_altitude_m,
                                      cos_zenith, earth_radius_m);
    const double cos_t = std::min(cos_zenith, 1.0);
    return upward ? slant_depth_upward_row(at, telescope_altitude_m, emission_altitude_m, cos_t,
                                           earth_radius_m)
                  : slant_depth_row(at, telescope_altitude_m, emission_altitude_m, cos_t,
                                    earth_radius_m);
  }

  // Slant optical depth of validated rows, cos_zenith in [0, 1]. Segments are the table
  // intervals between consecutive heights (tau linear in log10 h), clipped to [h_t, zem]; d tau/dh
  // = a / (h ln 10) with a the slope in log10 h. The slope of a linear-mode lookup is the blend of
  // the slopes of the two rows (tau is linear in lambda, so this is the depth of the blended tau).
  // Each segment is integrated with 8-point Gauss-Legendre; the first segment at h_t uses
  // h = h_t + u^2 when cos_zenith is small.
  // cos theta(h)^2 = [d (R + h + R + h_t) + (R + h_t)^2 cos_t^2] / (R + h)^2 with d = h - h_t,
  // which has no cancellation at small angles.
  [[nodiscard]] double slant_depth_row(const Lookup &at, double telescope_altitude_m,
                                       double emission_altitude_m, double cos_zenith,
                                       double earth_radius_m) const {
    if (!(emission_altitude_m > telescope_altitude_m)) {
      return 0.0;
    }
    return slant_integral(at, telescope_altitude_m, emission_altitude_m, telescope_altitude_m,
                          cos_zenith, earth_radius_m);
  }

  // Upward photon: the integral over [emission, h_t] (zero unless emission < h_t).
  [[nodiscard]] double slant_depth_upward_row(const Lookup &at, double telescope_altitude_m,
                                              double emission_altitude_m, double cos_zenith,
                                              double earth_radius_m) const {
    if (!(emission_altitude_m < telescope_altitude_m)) {
      return 0.0;
    }
    return slant_integral(at, emission_altitude_m, telescope_altitude_m, telescope_altitude_m,
                          cos_zenith, earth_radius_m);
  }

  // Integral over the heights [lo_altitude, hi_altitude] (lo < hi) of the slant depth for the ray
  // of the telescope altitude telescope_altitude_m (d = h - h_t in the expressions below).
  [[nodiscard]] double slant_integral(const Lookup &at, double lo_altitude, double hi_altitude,
                                      double telescope_altitude_m, double cos_zenith,
                                      double earth_radius_m) const {
    if (!(hi_altitude > lo_altitude)) {
      return 0.0;
    }
    const std::size_t n = grid_.size();
    const double *tau = &tau_[static_cast<std::size_t>(at.lower - first_) * n];
    const bool blend = at.upper != at.lower;
    const double *tau_up = tau + n;
    const double q = earth_radius_m + telescope_altitude_m;
    const double q2c2 = q * q * cos_zenith * cos_zenith;
    const bool substitute = cos_zenith < detail::curved_substitution_cos;
    const double ln10 = std::log(10.0);
    double sum = 0.0;
    for (std::size_t i = 0; i + 1 < grid_.size(); ++i) {
      const double seg_lo = std::max(height_m_[i], lo_altitude);
      const double seg_hi = std::min(height_m_[i + 1], hi_altitude);
      if (!(seg_hi > seg_lo)) {
        continue;
      }
      const double width = (grid_[i + 1] - grid_[i]) * ln10;
      double slope = (tau[i + 1] - tau[i]) / width;
      if (blend) {
        slope += at.frac * ((tau_up[i + 1] - tau_up[i]) / width - slope);
      }
      if (slope == 0.0) {
        continue;
      }
      // Integrand of dh at height h = h_t + d (d >= 0): slope / (h cos theta).
      auto integrand = [&](double d) {
        const double rho = q + d;
        const double cos2 = (d * (rho + q) + q2c2) / (rho * rho);
        return slope / ((telescope_altitude_m + d) * std::sqrt(std::max(cos2, 0.0)));
      };
      double part = 0.0;
      if (substitute && height_m_[i] <= lo_altitude) {
        // Segment starting at lo: h = lo + u^2, dh = 2 u du, u in [0, sqrt(seg_hi - lo)]. For a
        // downward photon lo = h_t, so seg_lo - h_t = 0 and the integrand is that of h = h_t + u^2.
        const double top = std::sqrt(seg_hi - seg_lo);
        const double d0 = seg_lo - telescope_altitude_m;
        for (std::size_t j = 0; j < detail::gauss_points; ++j) {
          const double u = 0.5 * top * (1.0 + detail::gauss_nodes[j]);
          part += detail::gauss_weights[j] * 2.0 * u * integrand(d0 + u * u);
        }
        sum += 0.5 * top * part;
      } else {
        const double mid = 0.5 * (seg_lo + seg_hi);
        const double half = 0.5 * (seg_hi - seg_lo);
        for (std::size_t j = 0; j < detail::gauss_points; ++j) {
          const double h = mid + half * detail::gauss_nodes[j];
          part += detail::gauss_weights[j] * integrand(h - telescope_altitude_m);
        }
        sum += half * part;
      }
    }
    return sum;
  }

  static double number(const std::string &token, std::size_t line, const char *what) {
    double value = 0.0;
    if (!obdeect::stream::parse::finite_number(token, value)) {
      throw std::invalid_argument("transmission table line " + std::to_string(line) + ": " + what +
                                  " '" + token + "' is not a finite number");
    }
    return value;
  }

  void parse(const std::string &text) {
    std::istringstream stream(text);
    std::string line;
    std::size_t line_number = 0;
    std::vector<double> heights_km;
    bool have_header = false;
    std::vector<double> values;
    while (std::getline(stream, line)) {
      ++line_number;
      const std::size_t h2 = line.find("H2=");
      const std::size_t h1 = line.find("H1=");
      if (!have_header && h2 != std::string::npos && h1 != std::string::npos && h2 < h1) {
        std::string head = line.substr(h2 + 3, h1 - h2 - 3);
        std::replace(head.begin(), head.end(), ',', ' ');
        std::istringstream head_stream(head);
        std::string token;
        if (!(head_stream >> token)) {
          throw std::invalid_argument("transmission table: empty H2");
        }
        h2_km_ = number(token, line_number, "H2");
        const std::string_view after_h1 =
            detail::trim_comment(std::string_view(line).substr(h1 + 3));
        std::istringstream list{std::string(after_h1)};
        while (list >> token) {
          heights_km.push_back(number(token, line_number, "height"));
        }
        have_header = true;
        continue;
      }
      std::istringstream row{std::string(detail::trim_comment(line))};
      std::string token;
      if (!(row >> token)) {
        continue;
      }
      if (!have_header) {
        throw std::invalid_argument(
            "transmission table: data before the '# H2= ... H1= ...' header");
      }
      const double wavelength = number(token, line_number, "wavelength");
      if (wavelength != std::floor(wavelength) || wavelength < 0.0 || wavelength >= max_lambda) {
        throw std::invalid_argument("transmission table line " + std::to_string(line_number) +
                                    ": wavelength must be an integer in [0, 1000) nm");
      }
      const int w = static_cast<int>(wavelength);
      if (rows_ == 0) {
        first_ = w;
      } else if (w != first_ + rows_) {
        throw std::invalid_argument("transmission table line " + std::to_string(line_number) +
                                    ": wavelengths must rise in steps of 1 nm");
      }
      const std::size_t before = values.size();
      while (row >> token) {
        const double depth = number(token, line_number, "optical depth");
        if (depth < 0.0) {
          throw std::invalid_argument("transmission table line " + std::to_string(line_number) +
                                      ": optical depth must not be negative");
        }
        values.push_back(depth);
      }
      if (values.size() - before != heights_km.size()) {
        throw std::invalid_argument("transmission table line " + std::to_string(line_number) +
                                    ": number of values differs from the number of heights");
      }
      ++rows_;
    }
    if (heights_km.size() + 1 > max_table_heights) { // H2 counts as a height
      throw std::invalid_argument(
          "transmission table: at most 100 heights (H2 and H1) are supported");
    }
    if (!have_header || heights_km.empty() || rows_ == 0) {
      throw std::invalid_argument("transmission table needs a header with H2 and H1 and data rows");
    }
    if (!(h2_km_ > 0.0)) {
      throw std::invalid_argument("transmission table: H2 must be positive");
    }
    double previous = h2_km_;
    for (double height : heights_km) {
      if (!(height > previous)) {
        throw std::invalid_argument("transmission table: heights must rise above H2");
      }
      previous = height;
    }
    height_m_.push_back(1000.0 * h2_km_);
    height_m_.insert(height_m_.end(), heights_km.size(), 0.0);
    for (std::size_t k = 0; k < heights_km.size(); ++k) {
      height_m_[k + 1] = 1000.0 * heights_km[k];
    }
    for (double height_m : height_m_) {
      grid_.push_back(std::log10(height_m));
    }
    const std::size_t columns = heights_km.size();
    tau_.assign(static_cast<std::size_t>(rows_) * grid_.size(), 0.0);
    for (std::size_t r = 0; r < static_cast<std::size_t>(rows_); ++r) {
      for (std::size_t c = 0; c < columns; ++c) {
        tau_[r * grid_.size() + c + 1] = values[r * columns + c];
      }
    }
  }

  WavelengthInterpolation mode_{WavelengthInterpolation::floor};
  double h2_km_{};
  int first_{};
  int rows_{};
  std::vector<double> grid_;     // log10(height / m), starting at the observation level
  std::vector<double> height_m_; // heights of grid_ [m]
  std::vector<double> tau_;      // rows x grid_.size(); first column 0
};

} // namespace obdeect::propagation
