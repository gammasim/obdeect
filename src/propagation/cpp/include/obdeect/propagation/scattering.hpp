#pragma once

// Optional single scattering of light-source photons (docs/scattering.md).
//
// Each primary photon (frame sky: position, direction, wavelength, weight at emission) is
// scattered at N_s points sampled along its downward path, by optical depth. For every
// point and every telescope sphere one scattered photon is aimed through a uniform point
// of the sphere's projected disc, with weight
//
//   w_s = w P_scatter(segment) [(tau_mol' p_R + omega tau_aer' p_aer) / (tau_mol' + omega
//   tau_aer')]
//         Omega_tel T(point -> telescope).
//
// Straight lines, plane-parallel atmosphere, single scattering, no polarisation, no refraction.
// The primary photon itself is not changed; its extinction is applied by the extinction step.

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

// Scattering table (plain text, docs/scattering.md): per integer wavelength the vertical
// profiles of the molecular and aerosol optical depth (from the altitude to the top), the
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

  // Relative tolerance of the trapezoid normalisation of a phase row (docs/scattering.md).
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
};

namespace detail {

// Altitude in [lo, hi] (m) where the molecular plus aerosol optical depth equals target.
// The depth is non-increasing in altitude, with depth(lo) >= target >= depth(hi). Bisection
// with a fixed number of steps, so the result does not depend on the caller.
template <std::floating_point T>
[[nodiscard]] T altitude_at_depth(const ScatteringTable &table, double lambda_nm, T target, T lo,
                                  T hi) {
  for (int step = 0; step < 64; ++step) {
    const T mid = T{0.5} * (lo + hi);
    const auto s = table.state<T>(lambda_nm, mid);
    if (s.tau_mol + s.tau_aer > target) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return T{0.5} * (lo + hi);
}

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

// Mixed phase bracket of one segment: (tau_mol p_R + omega tau_aer p_HG) / (tau_mol + omega
// tau_aer), with tau_* the segment's optical depths. Requires tau_mol + omega tau_aer > 0.
template <std::floating_point T>
[[nodiscard]] T mixed_phase(T cos_theta, T tau_mol, T tau_aer, T ssa, T asymmetry) {
  const T scattering = tau_mol + ssa * tau_aer;
  return (tau_mol * rayleigh_phase(cos_theta) +
          ssa * tau_aer * henyey_greenstein_phase(cos_theta, asymmetry)) /
         scattering;
}

// Mixed phase bracket with the aerosol phase value p_aer given (tabulated or Henyey-Greenstein)
// and the Rayleigh depolarisation gamma. Requires tau_mol + omega tau_aer > 0.
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

// Scattered photons of one primary photon (row of the chunk), appended to out[telescope].
// Sampling: the path from the emission altitude to the observation level is cut into N_s
// segments of equal slant optical depth D_j = j dD. In segment j the scattering point is
// drawn with density proportional to exp(-D), by inverse CDF, so the probability that the
// photon scatters in the segment is P = exp(-D_j) (1 - exp(-dD)) omega_mix. Each point gives
// one scattered photon per telescope.
inline void scatter_primary(const ScatteringConfig &config, const ScatteringTable &table,
                            const TransmissionTable &transmission,
                            std::span<const ScatterTelescope> telescopes,
                            const obdeect::stream::PhotonChunk &chunk, std::size_t row,
                            std::uint64_t photon, std::uint64_t seed, std::uint32_t event,
                            std::span<double> boundary,
                            std::vector<std::vector<std::array<float, n_scatter_columns>>> &out) {
  using obdeect::stream::Column;
  const double observation_level = transmission.observation_level_m();
  const double x0 = chunk[Column::x][row];
  const double y0 = chunk[Column::y][row];
  const double z0 = chunk[Column::z][row];
  const double cx = chunk[Column::cx][row];
  const double cy = chunk[Column::cy][row];
  const double t0 = chunk[Column::t][row];
  const double weight = chunk[Column::w][row];
  const double lambda = chunk[Column::wavelength][row];
  const double abs_cz = std::sqrt(std::max(0.0, 1.0 - cx * cx - cy * cy));
  const double h0 = observation_level + z0;
  if (weight <= 0.0 || abs_cz <= 0.0 || h0 <= observation_level) {
    return;
  }

  const std::size_t n = config.points_per_photon;
  const std::size_t n_telescopes = telescopes.size();
  const obdeect::stream::DrawKey key{seed, stage_atmosphere, event, 0, photon, stream_scattering};
  const auto extinction = [&](double altitude) {
    const auto s = table.state<double>(lambda, altitude);
    return s.tau_mol + s.tau_aer;
  };
  const double tau0 = extinction(h0);
  const double total = (extinction(observation_level) - tau0) / abs_cz;
  if (!(total > 0.0)) {
    return;
  }
  const double step = total / static_cast<double>(n);
  const double segment_weight = -std::expm1(-step); // 1 - exp(-dD)

  // Altitude of the segment boundaries D_j (boundary[0] = emission, boundary[n] = ground).
  boundary[0] = h0;
  boundary[n] = observation_level;
  for (std::size_t j = 1; j < n; ++j) {
    boundary[j] = altitude_at_depth<double>(
        table, lambda, tau0 + abs_cz * step * static_cast<double>(j), observation_level, h0);
  }

  for (std::size_t j = 0; j < n; ++j) {
    // Scattering point in segment j, by inverse CDF of exp(-D) on the segment.
    const double u = deviate(key, j);
    const double depth_in_segment = -std::log1p(-u * segment_weight);
    const double depth = step * static_cast<double>(j) + depth_in_segment;
    const double h_s = altitude_at_depth<double>(table, lambda, tau0 + abs_cz * depth,
                                                 boundary[j + 1], boundary[j]);
    const auto local = table.state<double>(lambda, h_s);
    const auto lower = table.state<double>(lambda, boundary[j + 1]);
    const auto upper = table.state<double>(lambda, boundary[j]);
    // Vertical optical depths of the segment (ratios only enter, so vertical is enough).
    const double tau_mol = lower.tau_mol - upper.tau_mol;
    const double tau_aer = lower.tau_aer - upper.tau_aer;
    const double scattering = tau_mol + local.ssa * tau_aer;
    if (!(tau_mol + tau_aer > 0.0) || !(scattering > 0.0)) {
      continue;
    }
    // P_scatter(segment) = exp(-D_j) (1 - exp(-dD)) omega_mix, omega_mix = (tau_mol + omega
    // tau_aer) / (tau_mol + tau_aer).
    const double p_scatter = std::exp(-step * static_cast<double>(j)) * segment_weight *
                             (scattering / (tau_mol + tau_aer));

    const double s_path = (h0 - h_s) / abs_cz; // emission to scattering point
    const double x_p = x0 + cx * s_path;
    const double y_p = y0 + cy * s_path;
    const double z_p = h_s - observation_level;

    for (std::size_t t = 0; t < n_telescopes; ++t) {
      const ScatterTelescope &tel = telescopes[t];
      const double h_c = observation_level + tel.z_m;
      if (!(h_s > h_c)) {
        continue; // the scattered photon must travel downwards
      }
      const double dx = tel.x_m - x_p;
      const double dy = tel.y_m - y_p;
      const double dz = tel.z_m - z_p;
      const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (!(distance > tel.radius_m)) {
        continue;
      }
      // Solid angle of the sphere seen from the point and transmission along the path.
      const double ratio = tel.radius_m / distance;
      const double solid_angle = 2.0 * std::numbers::pi * (1.0 - std::sqrt(1.0 - ratio * ratio));
      const double optical_depth = transmission.optical_depth<double>(lambda, h_s) -
                                   transmission.optical_depth<double>(lambda, h_c);
      const double airmass = distance / (h_s - h_c);
      const double transmitted = std::exp(-optical_depth * airmass);

      // Aim point: uniform over the projected disc of the sphere, perpendicular to the line
      // of sight through the centre.
      const std::array<double, 3> axis = {dx / distance, dy / distance, dz / distance};
      std::array<double, 3> e1{};
      std::array<double, 3> e2{};
      detail::perpendicular_basis(axis, e1, e2);
      const std::size_t draw = n + 2 * (j * n_telescopes + t);
      const double radial = tel.radius_m * std::sqrt(deviate(key, draw));
      const double phi = 2.0 * std::numbers::pi * deviate(key, draw + 1);
      const std::array<double, 3> offset = {
          radial * (std::cos(phi) * e1[0] + std::sin(phi) * e2[0]),
          radial * (std::cos(phi) * e1[1] + std::sin(phi) * e2[1]),
          radial * (std::cos(phi) * e1[2] + std::sin(phi) * e2[2])};
      const double ax = dx + offset[0];
      const double ay = dy + offset[1];
      const double az = dz + offset[2];
      const double aim = std::sqrt(ax * ax + ay * ay + az * az);
      const double nx = ax / aim;
      const double ny = ay / aim;
      const double nz = az / aim;

      // Scattering angle: between the propagation direction and the direction to the aim point.
      const double cos_theta = cx * nx + cy * ny - abs_cz * nz;
      const double p_aerosol = table.has_phase()
                                   ? table.aerosol_phase<double>(lambda, h_s, cos_theta)
                                   : henyey_greenstein_phase<double>(cos_theta, local.asymmetry);
      const double bracket = mixed_phase_with_aerosol<double>(
          cos_theta, tau_mol, tau_aer, local.ssa, table.rayleigh_gamma(), p_aerosol);
      const double scattered = weight * p_scatter * bracket * solid_angle * transmitted;
      if (!(scattered > 0.0)) {
        continue;
      }
      // Time at the scattering point; the stage adds the path from the point to the plane.
      double time_at_point = t0 + s_path / speed_of_light_m_per_ns;
      if (config.air_profile != nullptr) {
        // From |cz| directly: acos would round to pi/2 for |cz| below about 1e-16 (see
        // straight_path_air_delay_from_cos_ns). abs_cz > 0 here.
        time_at_point += straight_path_air_delay_from_cos_ns(*config.air_profile, h_s, h0, abs_cz);
      }
      out[t].push_back({static_cast<float>(x_p), static_cast<float>(y_p), static_cast<float>(z_p),
                        static_cast<float>(nx), static_cast<float>(ny), static_cast<float>(nz),
                        static_cast<float>(time_at_point), static_cast<float>(h_s),
                        static_cast<float>(scattered), static_cast<float>(lambda)});
    }
  }
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
// agree with cx, cy and z; no p column (the acceptance is applied later). Only downward photons are
// scattered.
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
    // cz_matches_direction is false for NaN; the finiteness test is kept as the first check.
    if (chunk.has(Column::cz)) {
      const double cz_column = chunk[Column::cz][row];
      if (!std::isfinite(cz_column) ||
          !obdeect::stream::cz_matches_direction(cz_column, cx[row], cy[row])) {
        throw std::invalid_argument("cz column must be finite and match cx, cy "
                                    "(cz = -sqrt(1-cx^2-cy^2))");
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
