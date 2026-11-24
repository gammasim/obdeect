#pragma once

// Optional single scattering of light-source photons (docs/propagation/scattering.md).
//
// Each primary photon (frame sky: position, direction, wavelength, weight at emission) is
// scattered at N_s points sampled along its path (downward to the observation level, or upward
// to the top of the table), by optical depth. For every
// point and every telescope sphere one scattered photon is aimed through a uniform point
// of the sphere's projected disc, with weight
//
//   w_s = w P_scatter(segment) [(tau_mol' p_R + omega tau_aer' p_aer) / (tau_mol' + omega
//   tau_aer')]
//         Omega_tel T(point -> telescope).
//
// Straight lines, plane-parallel atmosphere, single scattering, no polarisation, no refraction.
// The primary photon itself is not changed; its extinction is applied by the extinction step.

#include "obdeect/propagation/chord.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/streams.hpp"
#include "obdeect/propagation/thinning.hpp"
#include "obdeect/propagation/timing.hpp"
#include "obdeect/propagation/transmission.hpp"
#include "obdeect/stream/parse.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/random.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <numbers>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace obdeect::propagation {

inline constexpr std::uint32_t max_points_per_photon = 1024;

// Limit of the counter fields used for the draws of one photon: N_s + 2 N_s n_telescopes.
inline constexpr std::uint64_t max_scatter_deviates = std::uint64_t{1} << 26;

// Rayleigh phase function (3 / 16 pi) (1 + cos^2 theta), normalised over 4 pi.
template <std::floating_point T> [[nodiscard]] T rayleigh_phase(T cos_theta) {
  return T{3.0 / (16.0 * std::numbers::pi)} * (T{1} + cos_theta * cos_theta);
}

// Rayleigh phase function with depolarisation gamma = rho / (2 - rho):
// 3 / (8 pi (2 + gamma)) [(1 + gamma) + (1 - gamma) cos^2 theta], normalised over 4 pi.
// gamma = 0 is rayleigh_phase(cos_theta).
template <std::floating_point T> [[nodiscard]] T rayleigh_phase(T cos_theta, T gamma) {
  return T{3.0 / (8.0 * std::numbers::pi)} / (T{2} + gamma) *
         ((T{1} + gamma) + (T{1} - gamma) * cos_theta * cos_theta);
}

// Henyey-Greenstein phase function (1 - g^2) / (4 pi (1 + g^2 - 2 g cos theta)^{3/2}),
// normalised over 4 pi; mean cosine g. cos_theta is measured from the propagation direction.
template <std::floating_point T> [[nodiscard]] T henyey_greenstein_phase(T cos_theta, T asymmetry) {
  const T g2 = asymmetry * asymmetry;
  const T denominator = T{1} + g2 - T{2} * asymmetry * cos_theta;
  return (T{1} - g2) / (T{4.0 * std::numbers::pi} * denominator * std::sqrt(denominator));
}

// Atmospheric state at one wavelength and altitude. tau_* are vertical optical depths from
// the altitude to the top of the atmosphere.
template <std::floating_point T> struct ScatteringState {
  T tau_mol{};
  T tau_aer{};
  T ssa{};
  T asymmetry{};
};

// Scattering table (plain text, docs/propagation/scattering.md): per integer wavelength the
// vertical profiles of the molecular and aerosol optical depth (from the altitude to the top), the
// aerosol single-scattering albedo and asymmetry parameter. Linear in log10(altitude) and in
// wavelength between the table nodes, clamped to the table edges. Optional keyword lines: a
// tabulated aerosol phase function (angles_deg, phase) and the Rayleigh depolarisation
// (depolarization).
class ScatteringTable {
public:
  // Parse the text of a scattering table. Throws std::invalid_argument on invalid input.
  explicit ScatteringTable(std::string_view text) { parse(text); }

  // Read a scattering table from disk. Throws std::runtime_error if it cannot be read.
  static ScatteringTable from_file(const std::string &path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
      throw std::runtime_error("cannot read scattering table: " + path);
    }
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return ScatteringTable(buffer.str());
  }

  [[nodiscard]] int first_wavelength_nm() const noexcept { return wl_first_; }
  [[nodiscard]] int last_wavelength_nm() const noexcept {
    return wl_first_ + static_cast<int>(rows_.size()) - 1;
  }
  [[nodiscard]] double min_altitude_m() const noexcept { return altitude_m_.front(); }
  [[nodiscard]] double max_altitude_m() const noexcept { return altitude_m_.back(); }

  // True if the table has a tabulated aerosol phase function (aerosol_phase()).
  [[nodiscard]] bool has_phase() const noexcept { return !phase_.empty(); }

  // Rayleigh depolarisation factor gamma = rho / (2 - rho); 0 without a depolarization line.
  [[nodiscard]] double rayleigh_gamma() const noexcept { return gamma_; }

  // State at a wavelength [nm] (finite; clamped to the table) and altitude [m a.s.l.] (finite
  // and positive). Throws std::invalid_argument otherwise.
  template <std::floating_point T>
  [[nodiscard]] ScatteringState<T> state(double lambda_nm, T altitude_m) const {
    const Bracket aw = wavelength_bracket(lambda_nm);
    const Bracket ah = altitude_bracket(static_cast<double>(altitude_m));
    const ScatteringState<T> lower = row_state<T>(aw.lower, ah);
    if (aw.upper == aw.lower) {
      return lower;
    }
    const ScatteringState<T> upper = row_state<T>(aw.upper, ah);
    const T f = static_cast<T>(aw.fraction);
    return {lower.tau_mol + f * (upper.tau_mol - lower.tau_mol),
            lower.tau_aer + f * (upper.tau_aer - lower.tau_aer),
            lower.ssa + f * (upper.ssa - lower.ssa),
            lower.asymmetry + f * (upper.asymmetry - lower.asymmetry)};
  }

  // Tabulated aerosol phase function p_aer(theta) [sr^-1] at a wavelength [nm] (finite), an
  // altitude [m a.s.l.] (finite, positive) and cos(theta) (clamped to [-1, 1]). Linear in
  // cos(theta) between the angle nodes, in log10(altitude) between the altitude rows and in
  // wavelength between the integer rows. Requires has_phase(). Does not allocate.
  template <std::floating_point T>
  [[nodiscard]] T aerosol_phase(double lambda_nm, T altitude_m, T cos_theta) const {
    if (!has_phase()) {
      throw std::invalid_argument("scattering table: no tabulated phase function");
    }
    const Bracket aw = wavelength_bracket(lambda_nm);
    const Bracket ah = altitude_bracket(static_cast<double>(altitude_m));
    const Bracket ac =
        bracket_of(cos_nodes_, std::clamp(static_cast<double>(cos_theta), -1.0, 1.0));
    const std::size_t n_angles = cos_nodes_.size();
    const std::size_t n_altitudes = altitude_m_.size();
    const auto at_pair = [&](std::size_t iw, std::size_t ih) {
      const double *p = phase_.data() + (iw * n_altitudes + ih) * n_angles;
      return p[ac.lower] + ac.fraction * (p[ac.upper] - p[ac.lower]);
    };
    const auto at_wavelength = [&](std::size_t iw) {
      const double lower = at_pair(iw, ah.lower);
      return lower + ah.fraction * (at_pair(iw, ah.upper) - lower);
    };
    const double lower = at_wavelength(aw.lower);
    if (aw.upper == aw.lower) {
      return static_cast<T>(lower);
    }
    return static_cast<T>(lower + aw.fraction * (at_wavelength(aw.upper) - lower));
  }

  // Altitude [m a.s.l.] in [lo, hi] where the vertical optical depth tau_mol + tau_aer above it
  // equals target at a wavelength [nm]. The depth is non-increasing in altitude, with
  // depth(lo) >= target >= depth(hi); outside that range the nearer end is returned (as a bisection
  // would). Within one altitude row interval the depth is linear in log10(altitude) (both rows and
  // the wavelength blend are), so the solution is found by a bisection over the table's altitude
  // nodes inside [lo, hi] (O(log n) state evaluations) and one closed-form interpolation in
  // log10(altitude). Deterministic; requires 0 < lo <= hi.
  [[nodiscard]] double altitude_at_depth(double lambda_nm, double target, double lo,
                                         double hi) const {
    const auto depth_at = [&](double altitude) {
      const auto s = state<double>(lambda_nm, altitude);
      return s.tau_mol + s.tau_aer;
    };
    const double depth_lo = depth_at(lo);
    const double depth_hi = depth_at(hi);
    if (!(depth_lo > target)) {
      return lo;
    }
    if (depth_hi >= target) {
      return hi;
    }
    // Breakpoints: node 0 = lo, nodes 1..m = the table altitudes in (lo, hi), node m + 1 = hi.
    const auto first = static_cast<std::size_t>(
        std::upper_bound(altitude_m_.begin(), altitude_m_.end(), lo) - altitude_m_.begin());
    const auto end = static_cast<std::size_t>(
        std::lower_bound(altitude_m_.begin(), altitude_m_.end(), hi) - altitude_m_.begin());
    const std::size_t interior = end > first ? end - first : 0;
    const auto node = [&](std::size_t i) {
      if (i == 0) {
        return lo;
      }
      return i == interior + 1 ? hi : altitude_m_[first + i - 1];
    };
    // depth(node(a)) > target >= depth(node(b)) with b = a + 1 after the search.
    std::size_t a = 0, b = interior + 1;
    double depth_a = depth_lo, depth_b = depth_hi;
    while (b - a > 1) {
      const std::size_t mid = a + (b - a) / 2;
      const double depth_mid = depth_at(node(mid));
      if (depth_mid > target) {
        a = mid;
        depth_a = depth_mid;
      } else {
        b = mid;
        depth_b = depth_mid;
      }
    }
    const double log_a = std::log10(node(a)), log_b = std::log10(node(b));
    const double log_h = log_a + (target - depth_a) / (depth_b - depth_a) * (log_b - log_a);
    return std::clamp(std::pow(10.0, log_h), node(a), node(b));
  }

private:
  struct Row {
    std::vector<double> tau_mol;
    std::vector<double> tau_aer;
    std::vector<double> ssa;
    std::vector<double> asymmetry;
  };

  // Position of a value between two nodes: lower <= upper, fraction in [0, 1). Equal nodes
  // (clamped values) have fraction 0 so that edge values are returned exactly.
  struct Bracket {
    std::size_t lower{};
    std::size_t upper{};
    double fraction{};
  };

  // Relative tolerance of the trapezoid normalisation of a phase row
  // (docs/propagation/scattering.md).
  static constexpr double phase_norm_tolerance = 0.05;

  int wl_first_{};
  std::vector<double> altitude_m_{};
  std::vector<double> log_h_m_{};
  std::vector<Row> rows_{};
  std::vector<double> cos_nodes_{}; // cos(theta) of the angle nodes, ascending
  std::vector<double> phase_{};     // [wavelength row][altitude row][cos node] [sr^-1]
  double gamma_{0.0};

  // Bracket of a value on ascending nodes x, clamped to the first and last node.
  [[nodiscard]] static Bracket bracket_of(const std::vector<double> &x, double at) {
    const std::size_t n = x.size();
    if (at <= x.front()) {
      return {0, 0, 0.0};
    }
    if (at >= x.back()) {
      return {n - 1, n - 1, 0.0};
    }
    const auto upper = std::upper_bound(x.begin(), x.end(), at);
    const auto i = static_cast<std::size_t>(upper - x.begin()) - 1;
    return {i, i + 1, (at - x[i]) / (x[i + 1] - x[i])};
  }

  // Bracket of a wavelength [nm] over the integer rows (finite; clamped to the table).
  [[nodiscard]] Bracket wavelength_bracket(double lambda_nm) const {
    if (!std::isfinite(lambda_nm)) {
      throw std::invalid_argument("scattering table: wavelength must be finite");
    }
    const std::size_t n = rows_.size();
    const double last = static_cast<double>(n - 1);
    const double pos =
        std::clamp(lambda_nm - static_cast<double>(wl_first_), 0.0, last); // in row units
    if (pos >= last) {
      return {n - 1, n - 1, 0.0};
    }
    const auto lower = static_cast<std::size_t>(std::floor(pos));
    return {lower, lower + 1, pos - static_cast<double>(lower)};
  }

  // Bracket of an altitude [m a.s.l.] over the altitude rows, in log10(altitude).
  [[nodiscard]] Bracket altitude_bracket(double altitude_m) const {
    if (!std::isfinite(altitude_m) || !(altitude_m > 0.0)) {
      throw std::invalid_argument("scattering table: altitude must be finite and positive");
    }
    return bracket_of(log_h_m_, std::log10(altitude_m));
  }

  // Row state at one wavelength row, interpolated in altitude.
  template <std::floating_point T>
  [[nodiscard]] ScatteringState<T> row_state(std::size_t iw, const Bracket &ah) const {
    const Row &row = rows_[iw];
    const auto at = [&ah](const std::vector<double> &y) {
      const double lower = y[ah.lower];
      return static_cast<T>(lower + ah.fraction * (y[ah.upper] - lower));
    };
    return {at(row.tau_mol), at(row.tau_aer), at(row.ssa), at(row.asymmetry)};
  }

  static long parse_wavelength(std::string_view token) {
    const double wavelength = parse_number(token, "wavelength");
    if (wavelength != std::floor(wavelength) || wavelength < 1.0 ||
        wavelength > static_cast<double>(std::numeric_limits<int>::max())) {
      throw std::invalid_argument("scattering table: wavelengths must be integer nm in "
                                  "[1, INT_MAX]");
    }
    return static_cast<long>(wavelength);
  }

  void parse(std::string_view text) {
    struct Entry {
      long wavelength;
      double altitude_km, tau_mol, tau_aer, ssa, asymmetry;
    };
    struct PhaseEntry {
      long wavelength;
      double altitude_km;
      std::vector<double> values; // in the order of angles_deg
    };
    std::vector<Entry> entries;
    std::vector<double> angles_deg;
    std::vector<PhaseEntry> phase_entries;
    bool has_depolarization = false;
    double rho = 0.0;

    std::size_t pos = 0;
    while (pos < text.size()) {
      std::size_t end = text.find('\n', pos);
      if (end == std::string_view::npos) {
        end = text.size();
      }
      std::string_view line = text.substr(pos, end - pos);
      pos = end + 1;
      const auto hash = line.find('#');
      if (hash != std::string_view::npos) {
        line = line.substr(0, hash);
      }
      const auto tokens = split_whitespace(line);
      if (tokens.empty()) {
        continue;
      }
      if (std::isalpha(static_cast<unsigned char>(tokens[0].front()))) {
        const std::string_view keyword = tokens[0];
        if (keyword == "angles_deg") {
          if (!angles_deg.empty()) {
            throw std::invalid_argument("scattering table: angles_deg given more than once");
          }
          if (tokens.size() < 3 || tokens.size() > 722) {
            throw std::invalid_argument("scattering table: angles_deg needs 2 to 721 values");
          }
          for (std::size_t i = 1; i < tokens.size(); ++i) {
            angles_deg.push_back(parse_number(tokens[i], "angle"));
          }
          if (angles_deg.front() != 0.0 || angles_deg.back() != 180.0) {
            throw std::invalid_argument("scattering table: angles_deg must start at 0 and end "
                                        "at 180");
          }
          for (std::size_t i = 0; i < angles_deg.size(); ++i) {
            if (i > 0 && !(angles_deg[i] > angles_deg[i - 1])) {
              throw std::invalid_argument("scattering table: angles_deg must be strictly "
                                          "ascending");
            }
          }
        } else if (keyword == "phase") {
          if (tokens.size() < 4) {
            throw std::invalid_argument("scattering table: phase needs wavelength, altitude "
                                        "and at least one value");
          }
          PhaseEntry entry{parse_wavelength(tokens[1]), parse_number(tokens[2], "altitude"), {}};
          for (std::size_t i = 3; i < tokens.size(); ++i) {
            entry.values.push_back(parse_number(tokens[i], "phase value"));
          }
          phase_entries.push_back(std::move(entry));
        } else if (keyword == "depolarization") {
          if (has_depolarization) {
            throw std::invalid_argument("scattering table: depolarization given more than once");
          }
          if (tokens.size() != 2) {
            throw std::invalid_argument("scattering table: depolarization needs one value rho");
          }
          rho = parse_number(tokens[1], "depolarization");
          if (!(rho >= 0.0 && rho < 0.5)) {
            throw std::invalid_argument("scattering table: depolarization rho must be in [0, 0.5)");
          }
          has_depolarization = true;
        } else {
          throw std::invalid_argument("scattering table: unknown keyword '" + std::string(keyword) +
                                      "'");
        }
        continue;
      }
      if (tokens.size() != 6) {
        throw std::invalid_argument("scattering table: each row needs 6 values "
                                    "(wavelength altitude tau_mol tau_aer ssa g)");
      }
      entries.push_back({parse_wavelength(tokens[0]), parse_number(tokens[1], "altitude"),
                         parse_number(tokens[2], "tau_mol"), parse_number(tokens[3], "tau_aer"),
                         parse_number(tokens[4], "ssa"), parse_number(tokens[5], "g")});
    }
    if (entries.empty()) {
      throw std::invalid_argument("scattering table: no data rows");
    }

    // Altitudes: those of the first wavelength; every wavelength repeats them in order.
    const long first = entries.front().wavelength;
    std::vector<double> altitudes;
    for (const Entry &entry : entries) {
      if (entry.wavelength != first) {
        break;
      }
      altitudes.push_back(entry.altitude_km);
    }
    const std::size_t n_altitudes = altitudes.size();
    if (n_altitudes < 2 || entries.size() % n_altitudes != 0) {
      throw std::invalid_argument("scattering table: rows must be grouped by wavelength, "
                                  "each with the same altitudes (at least two)");
    }
    for (std::size_t i = 0; i < n_altitudes; ++i) {
      if (!(altitudes[i] > 0.0) || (i > 0 && !(altitudes[i] > altitudes[i - 1]))) {
        throw std::invalid_argument("scattering table: altitudes must be positive and ascending");
      }
    }
    const std::size_t n_wavelengths = entries.size() / n_altitudes;
    rows_.assign(n_wavelengths, Row{});
    for (std::size_t iwl = 0; iwl < n_wavelengths; ++iwl) {
      Row &row = rows_[iwl];
      row.tau_mol.resize(n_altitudes);
      row.tau_aer.resize(n_altitudes);
      row.ssa.resize(n_altitudes);
      row.asymmetry.resize(n_altitudes);
      for (std::size_t ih = 0; ih < n_altitudes; ++ih) {
        const Entry &entry = entries[iwl * n_altitudes + ih];
        if (entry.wavelength != first + static_cast<long>(iwl) ||
            entry.altitude_km != altitudes[ih]) {
          throw std::invalid_argument("scattering table: rows must be grouped by wavelength, "
                                      "each with the same altitudes");
        }
        if (!(entry.tau_mol >= 0.0) || !(entry.tau_aer >= 0.0)) {
          throw std::invalid_argument("scattering table: optical depths must be >= 0");
        }
        if (ih > 0 &&
            (entry.tau_mol > row.tau_mol[ih - 1] || entry.tau_aer > row.tau_aer[ih - 1])) {
          throw std::invalid_argument("scattering table: optical depths must not increase "
                                      "with altitude");
        }
        if (!(entry.ssa >= 0.0 && entry.ssa <= 1.0) ||
            !(entry.asymmetry > -1.0 && entry.asymmetry < 1.0)) {
          throw std::invalid_argument("scattering table: ssa must be in [0, 1], g in (-1, 1)");
        }
        row.tau_mol[ih] = entry.tau_mol;
        row.tau_aer[ih] = entry.tau_aer;
        row.ssa[ih] = entry.ssa;
        row.asymmetry[ih] = entry.asymmetry;
      }
    }

    // Tabulated phase function: every (wavelength, altitude) pair once, K values per row.
    if (angles_deg.empty() != phase_entries.empty()) {
      throw std::invalid_argument(angles_deg.empty()
                                      ? "scattering table: phase rows need angles_deg"
                                      : "scattering table: angles_deg needs phase rows");
    }
    std::vector<double> cos_nodes;
    std::vector<double> phase;
    if (!angles_deg.empty()) {
      const std::size_t n_angles = angles_deg.size();
      // Ascending cos(theta): the angle 180 deg comes first. Phase rows are stored in this order.
      cos_nodes.resize(n_angles);
      for (std::size_t k = 0; k < n_angles; ++k) {
        cos_nodes[k] = std::cos(angles_deg[n_angles - 1 - k] * std::numbers::pi / 180.0);
      }
      cos_nodes.front() = -1.0;
      cos_nodes.back() = 1.0;
      phase.assign(n_wavelengths * n_altitudes * n_angles, 0.0);
      std::vector<bool> seen(n_wavelengths * n_altitudes, false);
      for (const PhaseEntry &entry : phase_entries) {
        if (entry.values.size() != n_angles) {
          throw std::invalid_argument("scattering table: each phase row needs one value per "
                                      "angle (K values)");
        }
        if (entry.wavelength < first ||
            entry.wavelength >= first + static_cast<long>(n_wavelengths)) {
          throw std::invalid_argument("scattering table: phase row wavelength is not in the table");
        }
        const auto altitude = std::find(altitudes.begin(), altitudes.end(), entry.altitude_km);
        if (altitude == altitudes.end()) {
          throw std::invalid_argument("scattering table: phase row altitude is not in the table");
        }
        const auto pair = static_cast<std::size_t>(entry.wavelength - first) * n_altitudes +
                          static_cast<std::size_t>(altitude - altitudes.begin());
        if (seen[pair]) {
          throw std::invalid_argument("scattering table: duplicate phase row");
        }
        seen[pair] = true;
        double *out = phase.data() + pair * n_angles;
        for (std::size_t k = 0; k < n_angles; ++k) {
          const double value = entry.values[n_angles - 1 - k];
          if (!(value > 0.0)) {
            throw std::invalid_argument("scattering table: phase values must be > 0");
          }
          out[k] = value;
        }
        // Trapezoid rule in cos(theta): 2 pi int p d(cos theta) must be 1 within the tolerance.
        double norm = 0.0;
        for (std::size_t k = 0; k + 1 < n_angles; ++k) {
          norm += 0.5 * (out[k] + out[k + 1]) * (cos_nodes[k + 1] - cos_nodes[k]);
        }
        norm *= 2.0 * std::numbers::pi;
        if (!(std::abs(norm - 1.0) <= phase_norm_tolerance)) {
          throw std::invalid_argument("scattering table: phase row normalisation must be within "
                                      "5 % of 1 (2 pi int p d(cos theta))");
        }
        for (std::size_t k = 0; k < n_angles; ++k) {
          out[k] /= norm;
        }
      }
      if (std::find(seen.begin(), seen.end(), false) != seen.end()) {
        throw std::invalid_argument("scattering table: phase rows must cover every (wavelength, "
                                    "altitude) pair");
      }
    }

    wl_first_ = static_cast<int>(first);
    altitude_m_.clear();
    log_h_m_.clear();
    for (const double h_km : altitudes) {
      altitude_m_.push_back(h_km * 1000.0);
      log_h_m_.push_back(std::log10(h_km * 1000.0));
    }
    cos_nodes_ = std::move(cos_nodes);
    phase_ = std::move(phase);
    gamma_ = rho / (2.0 - rho); // 0 without a depolarization line
  }

  static std::vector<std::string_view> split_whitespace(std::string_view text) {
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

  static double parse_number(std::string_view token, const char *what) {
    double value = 0.0;
    if (!obdeect::stream::parse::finite_number(token, value)) {
      throw std::invalid_argument(std::string("scattering table: invalid ") + what + " '" +
                                  std::string(token) + "'");
    }
    return value;
  }
};

// A telescope sphere in the sky frame: centre (x, y, z above the observation level) and
// radius [m]. id is the telescope_id of the output photons.
struct ScatterTelescope {
  std::uint32_t id{};
  double x_m{};
  double y_m{};
  double z_m{};
  double radius_m{};
};

struct ScatteringConfig {
  bool enabled{false};                // off: scatter_chunk returns no photons
  std::uint32_t points_per_photon{8}; // N_s scattering points sampled per primary photon
  // Air of the stage (not owned, must outlive the call). Null: vacuum timing. Otherwise the
  // time at a scattering point adds the air delay of the straight path from the point up to the
  // emission height at the primary's zenith angle (straight_path_air_delay_ns).
  const AtmosphereProfile *air_profile{nullptr};
  // Sphere runs (docs/propagation/scattering.md, "Sphere runs"): a scattered photon whose path to
  // the telescope centre is near horizontal (|h_s - h_c| / distance < near_horizontal_cos) takes
  // the chord integral of the extinction for its transmission (no h_s == h_c skip), and a primary
  // with |cz| < near_horizontal_cos is scattered along its straight path (scatter_horizontal).
  // Plane runs: false (unchanged behaviour).
  bool sphere_geometry{false};
};

namespace detail {

// Orthonormal vectors e1, e2 perpendicular to the unit vector u.
inline void perpendicular_basis(const std::array<double, 3> &u, std::array<double, 3> &e1,
                                std::array<double, 3> &e2) {
  const std::array<double, 3> a = std::abs(u[0]) < 0.9 ? std::array<double, 3>{1.0, 0.0, 0.0}
                                                       : std::array<double, 3>{0.0, 1.0, 0.0};
  e1 = {a[1] * u[2] - a[2] * u[1], a[2] * u[0] - a[0] * u[2], a[0] * u[1] - a[1] * u[0]};
  const double norm = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  for (double &component : e1) {
    component /= norm;
  }
  e2 = {u[1] * e1[2] - u[2] * e1[1], u[2] * e1[0] - u[0] * e1[2], u[0] * e1[1] - u[1] * e1[0]};
}

} // namespace detail

// Mixed phase bracket (tau_mol p_R + omega tau_aer p_aer) / (tau_mol + omega tau_aer) of one
// segment, with the aerosol phase value p_aer given (tabulated or Henyey-Greenstein) and the
// Rayleigh depolarisation gamma (gamma = 0 for the unpolarised Rayleigh phase). Requires tau_mol +
// omega tau_aer > 0.
template <std::floating_point T>
[[nodiscard]] T mixed_phase_with_aerosol(T cos_theta, T tau_mol, T tau_aer, T ssa, T gamma,
                                         T p_aer) {
  const T scattering = tau_mol + ssa * tau_aer;
  return (tau_mol * rayleigh_phase(cos_theta, gamma) + ssa * tau_aer * p_aer) / scattering;
}

inline constexpr std::size_t n_scatter_columns = 10;

// Output columns of the scattered photons, in this order.
inline constexpr std::array<obdeect::stream::Column, n_scatter_columns> scatter_output_columns = {
    obdeect::stream::Column::x,         obdeect::stream::Column::y,   obdeect::stream::Column::z,
    obdeect::stream::Column::cx,        obdeect::stream::Column::cy,  obdeect::stream::Column::cz,
    obdeect::stream::Column::t,         obdeect::stream::Column::zem, obdeect::stream::Column::w,
    obdeect::stream::Column::wavelength};

namespace detail {

// One primary photon of the kernel: position and direction (cz signed: negative downward,
// positive upward, |cz| = abs_cz > 0), emission altitude h0 and the observation level (m a.s.l.).
struct Primary {
  double x0{};
  double y0{};
  double cx{};
  double cy{};
  double cz{};
  double abs_cz{};
  double t0{};
  double weight{};
  double lambda{};
  double h0{};
  double observation_level{};
};

// A scattering point of a primary: position (z above the observation level), altitude a.s.l.,
// time at the point, the probability of scattering in its segment, the segment's vertical optical
// depths and the albedo and asymmetry at the point.
struct ScatterPoint {
  double x{};
  double y{};
  double z{};
  double h{};
  double time{};
  double p_scatter{};
  double tau_mol{};
  double tau_aer{};
  double ssa{};
  double asymmetry{};
};

using ScatteredLists = std::vector<std::vector<std::array<float, n_scatter_columns>>>;

// Scattered photons of one scattering point towards every telescope sphere, appended to
// out[telescope]. The point lies in segment j of n points per photon; the aim point of telescope
// t uses the draws n + 2 (j n_telescopes + t) and the next one. A point below the telescope
// (h_s < h_c) gives an upward scattered photon, with optical depth |tau(h_s) - tau(h_c)| and
// airmass distance / |h_s - h_c|. A point at h_s == h_c would emit horizontal photons, which are
// not part of the sphere's phase function: they are skipped.
inline void scatter_to_telescopes(const ScatteringTable &table,
                                  const TransmissionTable &transmission,
                                  std::span<const ScatterTelescope> telescopes,
                                  const Primary &primary, const ScatterPoint &point, std::size_t n,
                                  std::size_t j, const obdeect::stream::DrawKey &key,
                                  bool sphere_geometry, ScatteredLists &out) {
  const std::size_t n_telescopes = telescopes.size();
  const double lambda = primary.lambda;
  const double h_s = point.h;
  for (std::size_t t = 0; t < n_telescopes; ++t) {
    const ScatterTelescope &tel = telescopes[t];
    const double h_c = primary.observation_level + tel.z_m;
    if (h_s == h_c && !sphere_geometry) {
      continue;
    }
    const double dx = tel.x_m - point.x;
    const double dy = tel.y_m - point.y;
    const double dz = tel.z_m - point.z;
    const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!(distance > tel.radius_m)) {
      continue;
    }
    // Solid angle of the sphere seen from the point and transmission along the path. The
    // transmission table counts from the observation level, so its depth grows with altitude;
    // the absolute value covers both directions.
    const double ratio = tel.radius_m / distance;
    const double solid_angle = 2.0 * std::numbers::pi * (1.0 - std::sqrt(1.0 - ratio * ratio));
    double transmitted = 0.0;
    if (sphere_geometry && std::abs(h_s - h_c) < near_horizontal_cos * distance) {
      // Near-horizontal path: the chord integral from the point to the centre (the plane form
      // has airmass distance / |h_s - h_c|, which is singular at h_s == h_c).
      transmitted = std::exp(-flat_chord_depth(transmission, lambda, h_s, h_c, distance));
    } else {
      const double optical_depth = std::abs(transmission.optical_depth<double>(lambda, h_s) -
                                            transmission.optical_depth<double>(lambda, h_c));
      const double airmass = distance / std::abs(h_s - h_c);
      transmitted = std::exp(-optical_depth * airmass);
    }

    // Aim point: uniform over the projected disc of the sphere, perpendicular to the line of
    // sight through the centre.
    const std::array<double, 3> axis = {dx / distance, dy / distance, dz / distance};
    std::array<double, 3> e1{};
    std::array<double, 3> e2{};
    detail::perpendicular_basis(axis, e1, e2);
    const std::size_t draw = n + 2 * (j * n_telescopes + t);
    const double radial = tel.radius_m * std::sqrt(deviate(key, draw));
    const double phi = 2.0 * std::numbers::pi * deviate(key, draw + 1);
    const std::array<double, 3> offset = {radial * (std::cos(phi) * e1[0] + std::sin(phi) * e2[0]),
                                          radial * (std::cos(phi) * e1[1] + std::sin(phi) * e2[1]),
                                          radial * (std::cos(phi) * e1[2] + std::sin(phi) * e2[2])};
    const double ax = dx + offset[0];
    const double ay = dy + offset[1];
    const double az = dz + offset[2];
    const double aim = std::sqrt(ax * ax + ay * ay + az * az);
    const double nx = ax / aim;
    const double ny = ay / aim;
    const double nz = az / aim;

    // Scattering angle between the propagation direction (signed cz) and the direction to the
    // aim point.
    const double cos_theta = primary.cx * nx + primary.cy * ny + primary.cz * nz;
    const double p_aerosol = table.has_phase()
                                 ? table.aerosol_phase<double>(lambda, h_s, cos_theta)
                                 : henyey_greenstein_phase<double>(cos_theta, point.asymmetry);
    const double bracket = mixed_phase_with_aerosol<double>(
        cos_theta, point.tau_mol, point.tau_aer, point.ssa, table.rayleigh_gamma(), p_aerosol);
    const double scattered = primary.weight * point.p_scatter * bracket * solid_angle * transmitted;
    if (!(scattered > 0.0)) {
      continue;
    }
    out[t].push_back({static_cast<float>(point.x), static_cast<float>(point.y),
                      static_cast<float>(point.z), static_cast<float>(nx), static_cast<float>(ny),
                      static_cast<float>(nz), static_cast<float>(point.time),
                      static_cast<float>(h_s), static_cast<float>(scattered),
                      static_cast<float>(lambda)});
  }
}

// Segment j of a path, from its altitudes: local is the state at the point, lower and upper the
// states at the segment's lower and upper boundary. Returns false if the segment is empty.
// p_scatter = exp(-D_j) (1 - exp(-dD)) omega_mix (docs/propagation/scattering.md).
inline bool segment_point(const ScatteringState<double> &local,
                          const ScatteringState<double> &lower,
                          const ScatteringState<double> &upper, double step, double segment_weight,
                          std::size_t j, ScatterPoint &point) {
  // Vertical optical depths of the segment (ratios only enter, so vertical is enough).
  const double tau_mol = lower.tau_mol - upper.tau_mol;
  const double tau_aer = lower.tau_aer - upper.tau_aer;
  const double scattering = tau_mol + local.ssa * tau_aer;
  if (!(tau_mol + tau_aer > 0.0) || !(scattering > 0.0)) {
    return false;
  }
  point.p_scatter = std::exp(-step * static_cast<double>(j)) * segment_weight *
                    (scattering / (tau_mol + tau_aer));
  point.tau_mol = tau_mol;
  point.tau_aer = tau_aer;
  point.ssa = local.ssa;
  point.asymmetry = local.asymmetry;
  return true;
}

// Total optical depth of the scattering table from altitude h to the top [m a.s.l.].
inline double table_depth(const ScatteringTable &table, double lambda, double h) {
  const auto s = table.state<double>(lambda, h);
  return s.tau_mol + s.tau_aer;
}

// Scattering points of a vertical primary: downward (cz < 0) from the emission altitude h0 down to
// the observation level, or upward (cz > 0) from h0 up to the top of the scattering table, where
// the scattering stops. The path is cut into n segments of equal slant depth D, counted with the
// optical depth tau(h) above an altitude h (tau grows downwards: the depth from h0 is added
// downward and subtracted upward). boundary[j] is the altitude at depth j dD from h0,
// boundary[0] = h0 and boundary[n] = the end of the path. The point of segment j is drawn by the
// inverse of exp(-D) within the segment.
inline void scatter_vertical(const ScatteringConfig &config, const ScatteringTable &table,
                             const TransmissionTable &transmission,
                             std::span<const ScatterTelescope> telescopes, const Primary &primary,
                             const obdeect::stream::DrawKey &key, bool upward,
                             std::span<double> boundary, ScatteredLists &out) {
  const std::size_t n = config.points_per_photon;
  const double lambda = primary.lambda;
  const double h0 = primary.h0;
  const double abs_cz = primary.abs_cz;
  const double end = upward ? table.max_altitude_m() : primary.observation_level;
  const double tau0 = table_depth(table, lambda, h0);
  const double total = upward ? (tau0 - table_depth(table, lambda, end)) / abs_cz
                              : (table_depth(table, lambda, end) - tau0) / abs_cz;
  if (!(total > 0.0)) {
    return;
  }
  const double step = total / static_cast<double>(n);
  const double segment_weight = -std::expm1(-step); // 1 - exp(-dD)
  // The altitude bracket of the boundaries and of the scattering points: [lower, upper].
  const double lower = upward ? h0 : end;
  const double upper = upward ? end : h0;

  boundary[0] = h0;
  boundary[n] = end;
  for (std::size_t j = 1; j < n; ++j) {
    const double offset = abs_cz * step * static_cast<double>(j);
    boundary[j] =
        table.altitude_at_depth(lambda, upward ? tau0 - offset : tau0 + offset, lower, upper);
  }

  for (std::size_t j = 0; j < n; ++j) {
    // Scattering point in segment j, by inverse CDF of exp(-D) on the segment.
    const double u = deviate(key, j);
    const double depth_in_segment = -std::log1p(-u * segment_weight);
    const double depth = step * static_cast<double>(j) + depth_in_segment;
    const double seg_lower = upward ? boundary[j] : boundary[j + 1];
    const double seg_upper = upward ? boundary[j + 1] : boundary[j];
    const double h_s = table.altitude_at_depth(
        lambda, upward ? tau0 - abs_cz * depth : tau0 + abs_cz * depth, seg_lower, seg_upper);
    ScatterPoint point;
    if (!segment_point(table.state<double>(lambda, h_s), table.state<double>(lambda, seg_lower),
                       table.state<double>(lambda, seg_upper), step, segment_weight, j, point)) {
      continue;
    }
    const double s_path = (upward ? h_s - h0 : h0 - h_s) / abs_cz; // emission to the point
    point.x = primary.x0 + primary.cx * s_path;
    point.y = primary.y0 + primary.cy * s_path;
    point.z = h_s - primary.observation_level;
    point.h = h_s;
    point.time = primary.t0 + s_path / speed_of_light_m_per_ns;
    if (config.air_profile != nullptr) {
      // From |cz| directly: acos would round to pi/2 for |cz| below about 1e-16 (see
      // straight_path_air_delay_from_cos_ns). abs_cz > 0 here.
      const double air_low = upward ? h0 : h_s, air_high = upward ? h_s : h0;
      point.time +=
          straight_path_air_delay_from_cos_ns(*config.air_profile, air_low, air_high, abs_cz);
    }
    scatter_to_telescopes(table, transmission, telescopes, primary, point, n, j, key,
                          config.sphere_geometry, out);
  }
}

// Sphere runs: primaries with |cz| < near_horizontal_cos (cz = 0: horizontal) are scattered along
// their straight path (no bending of the primary; the point is at x0 + c s, altitude h0 + cz s).
// Extinction depth from the emission to path length s: D(s) = |tau(h0 + cz s) - tau(h0)| / |cz|,
// with tau the scattering table's optical depth above an altitude (exact for a straight path with
// altitude linear in s). For |cz| < horizontal_cz_exact_limit the altitude is constant to the
// path's precision and D(s) = k(h0) s, k = -d tau/dh from a central difference over 1 m; then the
// inverse is exact. Otherwise D is tabulated at horizontal_nodes + 1 path lengths (uniform in s up
// to the path's end) and inverted by linear interpolation (D is monotone). The path ends where the
// altitude leaves [observation level, table top] (ground or top), where D reaches
// horizontal_depth_cut (exp(-D) = 1e-9; the rest is not sampled) or at horizontal_path_cap_m.
// Segments of equal depth: boundaries by the inverse; the point of segment j is drawn by the
// inverse of exp(-D) within the segment, as in scatter_downward. The molecular and aerosol mix of
// a segment is the local one (central differences at its mid-altitude): along a near-horizontal
// path the segment's own vertical depths are not enough (they are zero for cz = 0).
inline constexpr double horizontal_depth_cut = 20.72326583694641; // -ln(1e-9)
inline constexpr double horizontal_path_cap_m = 1e9;
inline constexpr double horizontal_cz_exact_limit = 1e-9;
inline constexpr std::size_t horizontal_nodes = 256;

// Central difference over 1 m of the optical depths (molecular, aerosol) at altitude h [m]: the
// decrease per metre. Within 1 m of the table's ends the difference is one-sided, so the table is
// only evaluated inside [min_altitude, max_altitude] (the altitude is not extended below the
// observation level or above the top of the table).
inline std::array<double, 2> depth_slope_per_m(const ScatteringTable &table, double lambda,
                                               double h, double min_altitude, double max_altitude) {
  const double lower = std::max(min_altitude, h - 1.0);
  const double upper = std::min(max_altitude, h + 1.0);
  if (!(upper > lower)) {
    return {0.0, 0.0};
  }
  const auto s_lower = table.state<double>(lambda, lower);
  const auto s_upper = table.state<double>(lambda, upper);
  const double inverse = 1.0 / (upper - lower);
  return {(s_lower.tau_mol - s_upper.tau_mol) * inverse,
          (s_lower.tau_aer - s_upper.tau_aer) * inverse};
}

// Extinction depth [1] of a straight primary from emission altitude h0 with signed cz, as a
// function of the path length s (see the comment above), and its inverse on the path. The
// inverse is exact: the depth D(s) = |tau(h0 + cz s) - tau0| / |cz| is inverted through
// ScatteringTable::altitude_at_depth on the path's altitude range.
class HorizontalDepth {
public:
  // s_limit: the path length where the altitude leaves the atmosphere (may be large). Throws
  // nothing; the caller checks s_limit > 0.
  HorizontalDepth(const ScatteringTable &table, double lambda, double h0, double cz, double s_limit)
      : table_(table), lambda_(lambda), h0_(h0), cz_(cz) {
    if (std::abs(cz) < horizontal_cz_exact_limit) {
      flat_ = true;
      const auto slope =
          depth_slope_per_m(table, lambda, h0, table.min_altitude_m(), table.max_altitude_m());
      k_ = slope[0] + slope[1];
      s_end_ = std::min(s_limit, horizontal_depth_cut / std::max(k_, 1e-300));
      d_end_ = k_ * s_end_;
      if (!(k_ > 0.0)) {
        s_end_ = 0.0;
        d_end_ = 0.0;
      }
      return;
    }
    tau0_ = table_depth(table, lambda, h0);
    const double end_altitude = h0 + cz * s_limit;
    lo_altitude_ = std::min(h0, end_altitude);
    hi_altitude_ = std::max(h0, end_altitude);
    s_end_ = s_limit;
    if (depth(s_limit) > horizontal_depth_cut) {
      s_end_ = length_for_depth(horizontal_depth_cut);
    }
    d_end_ = depth(s_end_);
  }

  // Depth at the end of the path (the cut, or the path's end at ground or top).
  [[nodiscard]] double depth_end() const noexcept { return d_end_; }
  // Path length of the end [m].
  [[nodiscard]] double length_end() const noexcept { return s_end_; }
  // Depth from the emission to path length s [m], 0 for s <= 0.
  [[nodiscard]] double depth(double s) const {
    if (!(s > 0.0)) {
      return 0.0;
    }
    if (flat_) {
      return k_ * s;
    }
    return std::abs(table_depth(table_, lambda_, h0_ + cz_ * s) - tau0_) / std::abs(cz_);
  }
  // Path length at which the depth equals d, for 0 <= d <= depth_end() (d is clamped there).
  [[nodiscard]] double length_at(double d) const {
    const double clamped = std::clamp(d, 0.0, d_end_);
    return flat_ ? clamped / k_ : length_for_depth(clamped);
  }

private:
  // Exact inverse of the depth for a non-flat primary: tau(h) = tau0 - cz d, the altitude is found
  // in the path's altitude range [lo, hi].
  [[nodiscard]] double length_for_depth(double d) const {
    const double altitude =
        table_.altitude_at_depth(lambda_, tau0_ - cz_ * d, lo_altitude_, hi_altitude_);
    return (altitude - h0_) / cz_;
  }

  const ScatteringTable &table_;
  double lambda_;
  double h0_;
  double cz_;
  bool flat_{};
  double k_{};
  double tau0_{};
  double lo_altitude_{};
  double hi_altitude_{};
  double s_end_{};
  double d_end_{};
};

// Scattered photons of one horizontal (sphere runs, |cz| < near_horizontal_cos) primary. The
// boundary span has points_per_photon + 1 entries (path lengths of the segment boundaries).
inline void scatter_horizontal(const ScatteringConfig &config, const ScatteringTable &table,
                               const TransmissionTable &transmission,
                               std::span<const ScatterTelescope> telescopes, const Primary &primary,
                               const obdeect::stream::DrawKey &key, std::span<double> boundary,
                               ScatteredLists &out) {
  const std::size_t n = config.points_per_photon;
  const double lambda = primary.lambda;
  const double h0 = primary.h0;
  const double cz = primary.cz;
  const double top = table.max_altitude_m();
  double s_limit = horizontal_path_cap_m;
  if (cz < 0.0) {
    s_limit = std::min(s_limit, (h0 - primary.observation_level) / -cz);
  } else if (cz > 0.0) {
    s_limit = std::min(s_limit, (top - h0) / cz);
  }
  if (!(s_limit > 0.0)) {
    return;
  }
  const HorizontalDepth depth(table, lambda, h0, cz, s_limit);
  if (!(depth.depth_end() > 0.0)) {
    return;
  }
  const double d_end = depth.depth_end();
  const double step = d_end / static_cast<double>(n);
  const double segment_weight = -std::expm1(-step); // 1 - exp(-dD)
  boundary[0] = 0.0;
  boundary[n] = depth.length_end();
  for (std::size_t j = 1; j < n; ++j) {
    boundary[j] = depth.length_at(step * static_cast<double>(j));
  }

  for (std::size_t j = 0; j < n; ++j) {
    const double u = deviate(key, j);
    const double depth_in_segment = -std::log1p(-u * segment_weight);
    const double target = step * static_cast<double>(j) + depth_in_segment;
    const double s_p = std::clamp(depth.length_at(target), boundary[j], boundary[j + 1]);
    // Extinction mix at the middle of the segment (central differences of the table).
    const double h_mid = h0 + cz * 0.5 * (boundary[j] + boundary[j + 1]);
    const auto slope =
        depth_slope_per_m(table, lambda, h_mid, table.min_altitude_m(), table.max_altitude_m());
    const double k_mol = slope[0];
    const double k_aer = slope[1];
    const double k_total = k_mol + k_aer;
    if (!(k_total > 0.0)) {
      continue;
    }
    const double h_p = h0 + cz * s_p;
    const auto local = table.state<double>(lambda, h_p);
    ScatterPoint point;
    point.p_scatter = std::exp(-step * static_cast<double>(j)) * segment_weight *
                      ((k_mol + local.ssa * k_aer) / k_total);
    point.tau_mol = k_mol;
    point.tau_aer = k_aer;
    point.ssa = local.ssa;
    point.asymmetry = local.asymmetry;
    point.x = primary.x0 + primary.cx * s_p;
    point.y = primary.y0 + primary.cy * s_p;
    point.z = h_p - primary.observation_level;
    point.h = h_p;
    point.time = primary.t0 + s_p / speed_of_light_m_per_ns;
    if (config.air_profile != nullptr) {
      point.time += straight_line_air_delay_ns(*config.air_profile, h0, h_p, s_p, 1.0);
    }
    scatter_to_telescopes(table, transmission, telescopes, primary, point, n, j, key,
                          config.sphere_geometry, out);
  }
}

// Scattered photons of one primary photon (row of the chunk), appended to out[telescope].
// Horizontal primaries (cz = 0) are not scattered. Downward and upward primaries take their own
// branch (scatter_vertical, upward or downward); both use the same draws and sphere geometry.
inline void scatter_primary(const ScatteringConfig &config, const ScatteringTable &table,
                            const TransmissionTable &transmission,
                            std::span<const ScatterTelescope> telescopes,
                            const obdeect::stream::PhotonChunk &chunk, std::size_t row,
                            std::uint64_t photon, std::uint64_t seed, std::uint32_t event,
                            std::span<double> boundary, ScatteredLists &out) {
  using obdeect::stream::Column;
  Primary primary;
  primary.observation_level = transmission.observation_level_m();
  primary.x0 = chunk[Column::x][row];
  primary.y0 = chunk[Column::y][row];
  primary.cx = chunk[Column::cx][row];
  primary.cy = chunk[Column::cy][row];
  primary.t0 = chunk[Column::t][row];
  primary.weight = chunk[Column::w][row];
  primary.lambda = chunk[Column::wavelength][row];
  primary.abs_cz =
      std::sqrt(std::max(0.0, 1.0 - primary.cx * primary.cx - primary.cy * primary.cy));
  primary.h0 = primary.observation_level + chunk[Column::z][row];
  // Sign from the cz column (absent: downward). Plane runs: cz = 0 is a horizontal primary, not
  // scattered. Sphere runs: |cz| < near_horizontal_cos is scattered along its straight path.
  const float cz_column = chunk.has(Column::cz) ? chunk[Column::cz][row] : -1.0F;
  if (primary.weight <= 0.0 || primary.h0 <= primary.observation_level) {
    return;
  }
  const obdeect::stream::DrawKey key{seed, stage_atmosphere, event, 0, photon, stream_scattering};
  if (config.sphere_geometry && (cz_column == 0.0F || primary.abs_cz < near_horizontal_cos)) {
    primary.cz = cz_column == 0.0F ? 0.0 : (cz_column > 0.0F ? primary.abs_cz : -primary.abs_cz);
    scatter_horizontal(config, table, transmission, telescopes, primary, key, boundary, out);
    return;
  }
  if (cz_column == 0.0F || primary.abs_cz <= 0.0) {
    return;
  }
  const bool upward = cz_column > 0.0F;
  primary.cz = upward ? primary.abs_cz : -primary.abs_cz;
  scatter_vertical(config, table, transmission, telescopes, primary, key, upward, boundary, out);
}

} // namespace detail

// Scatter the primary photons of a sky chunk towards every telescope sphere. Appends the
// scattered photons to chunk_out[i] (one list per telescope, each chunk at most 2^16 rows,
// frame sky, telescope_id = telescopes[i].id). first_index is the index of the chunk's first
// photon among all primary photons of the event; draws are keyed on it, so the output does
// not depend on how the photons are split into chunks, as long as the chunks are passed in
// order. Output photons of each list are numbered from the rows already in that list.
//
// Requires columns x y z cx cy t w wavelength (wavelength in nm within the scattering table,
// interpolated between its integer rows: assign it before this step); cz and zem, if present, must
// agree with cx, cy and z; no p column (the acceptance is applied later). Upward and downward
// photons are scattered; horizontal ones (cz = 0) are not.
inline void scatter_chunk(const ScatteringConfig &config, const ScatteringTable &table,
                          const TransmissionTable &transmission,
                          std::span<const ScatterTelescope> telescopes,
                          const obdeect::stream::PhotonChunk &chunk, std::uint64_t seed,
                          std::uint32_t event, std::uint64_t first_index,
                          std::vector<std::vector<obdeect::stream::PhotonChunk>> &chunk_out) {
  using obdeect::stream::Column;
  if (!config.enabled) {
    return;
  }
  if (config.points_per_photon < 1 || config.points_per_photon > max_points_per_photon) {
    throw std::invalid_argument("scattering needs 1 to 1024 points per photon");
  }
  for (std::size_t i = 0; i < telescopes.size(); ++i) {
    const ScatterTelescope &tel = telescopes[i];
    if (!std::isfinite(tel.x_m) || !std::isfinite(tel.y_m) || !std::isfinite(tel.z_m) ||
        !std::isfinite(tel.radius_m) || !(tel.radius_m > 0.0) || tel.z_m < 0.0) {
      throw std::invalid_argument("telescope needs finite x, y, z >= 0 (above the observation "
                                  "level) and finite radius > 0");
    }
    for (std::size_t k = 0; k < i; ++k) {
      if (telescopes[k].id == tel.id) {
        throw std::invalid_argument("telescope ids must be unique");
      }
    }
  }
  if (chunk.frame != obdeect::stream::Frame::sky) {
    throw std::invalid_argument("scattering needs photons in frame sky");
  }
  chunk.validate();
  for (const Column required : {Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::t,
                                Column::w, Column::wavelength}) {
    if (!chunk.has(required)) {
      throw std::invalid_argument("scattering needs the columns x y z cx cy t w wavelength");
    }
  }
  if (chunk.has(Column::p)) {
    throw std::invalid_argument("scattering needs photons before acceptance (no p column)");
  }
  if (chunk_out.size() < telescopes.size()) {
    chunk_out.resize(telescopes.size());
  }
  const std::size_t rows = chunk.rows();
  if (rows == 0 || telescopes.empty()) {
    return;
  }
  const std::size_t n_points = config.points_per_photon;
  if (static_cast<std::uint64_t>(n_points) * (1 + 2 * telescopes.size()) > max_scatter_deviates) {
    throw std::invalid_argument("too many scattering points and telescopes for one photon");
  }
  const double observation_level = transmission.observation_level_m();
  const auto &x = chunk[Column::x];
  const auto &y = chunk[Column::y];
  const auto &cx = chunk[Column::cx];
  const auto &cy = chunk[Column::cy];
  const auto &z = chunk[Column::z];
  const auto &t = chunk[Column::t];
  const auto &w = chunk[Column::w];
  const auto &wavelength = chunk[Column::wavelength];
  for (std::size_t row = 0; row < rows; ++row) {
    if (!obdeect::stream::direction_in_unit_disc(cx[row], cy[row])) {
      throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
    }
    // cz_is_consistent is false for NaN; the finiteness test is kept as the first check.
    if (chunk.has(Column::cz)) {
      // Either sign (upward and horizontal primaries too; cz = 0 is not scattered).
      const double cz_column = chunk[Column::cz][row];
      if (!std::isfinite(cz_column) ||
          !obdeect::stream::cz_is_consistent(cz_column, cx[row], cy[row])) {
        throw std::invalid_argument("cz column must be finite and match cx, cy "
                                    "(cz = +-sqrt(1-cx^2-cy^2))");
      }
    }
    if (chunk.has(Column::zem)) {
      const double zem_column = chunk[Column::zem][row];
      if (!std::isfinite(zem_column) || std::abs(zem_column - (observation_level + z[row])) > 1.0) {
        throw std::invalid_argument("zem column must be finite and match z plus the observation "
                                    "level");
      }
    }
    if (!(w[row] >= 0.0F) || !std::isfinite(w[row]) || !std::isfinite(z[row]) ||
        !std::isfinite(x[row]) || !std::isfinite(y[row]) || !std::isfinite(t[row])) {
      throw std::invalid_argument("weights, positions and times must be finite, weights >= 0");
    }
    const double lambda = wavelength[row];
    if (!(lambda >= table.first_wavelength_nm() && lambda < table.last_wavelength_nm() + 1.0)) {
      throw std::invalid_argument("wavelength must lie within the scattering table (nm)");
    }
  }

  // Scattered photons per telescope, in the order of their primary photons.
  std::vector<std::vector<std::array<float, n_scatter_columns>>> out(telescopes.size());
  std::vector<double> boundary(n_points + 1);
  for (std::size_t row = 0; row < rows; ++row) {
    detail::scatter_primary(config, table, transmission, telescopes, chunk, row, first_index + row,
                            seed, event, boundary, out);
  }

  // Pack into chunks of at most 2^16 rows, numbered after the rows already in each list.
  std::uint32_t mask = 0;
  for (const Column column : scatter_output_columns) {
    mask |= obdeect::stream::column_bit(column);
  }
  for (std::size_t i = 0; i < telescopes.size(); ++i) {
    std::uint64_t offset = 0;
    for (const auto &existing : chunk_out[i]) {
      offset += existing.rows();
    }
    const auto &photons = out[i];
    for (std::size_t first = 0; first < photons.size(); first += obdeect::stream::max_chunk_rows) {
      const std::size_t n =
          std::min<std::size_t>(obdeect::stream::max_chunk_rows, photons.size() - first);
      obdeect::stream::PhotonChunk piece;
      piece.telescope_id = telescopes[i].id;
      piece.frame = obdeect::stream::Frame::sky;
      piece.first_index = offset + first;
      piece.resize(n, mask);
      for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t c = 0; c < n_scatter_columns; ++c) {
          piece[scatter_output_columns[c]][j] = photons[first + j][c];
        }
      }
      chunk_out[i].push_back(std::move(piece));
    }
  }
}

} // namespace obdeect::propagation
