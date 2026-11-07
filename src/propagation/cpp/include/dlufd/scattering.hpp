#pragma once

// Optional single scattering of light-source photons (docs/scattering.md).
//
// Each primary photon (frame sky: position, direction, wavelength, weight at emission) is
// scattered at N_s points sampled along its downward path, by optical depth. For every
// point and every telescope sphere one scattered photon is aimed through a uniform point
// of the sphere's projected disc, with weight
//
//   w_s = w P_scatter(segment) [(tau_mol' p_R + omega tau_aer' p_HG) / (tau_mol' + omega tau_aer')]
//         Omega_tel T(point -> telescope).
//
// Straight lines, plane-parallel atmosphere, single scattering, no polarisation, no refraction.
// The primary photon itself is not changed; its extinction is applied by the extinction step.

#include "dlufd/streams.hpp"
#include "dlufd/thinning.hpp"
#include "dlufd/timing.hpp"
#include "dlufd/transmission.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/random.hpp"

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

namespace dlufd {

inline constexpr std::uint32_t max_points_per_photon = 1024;

// Limit of the counter fields used for the draws of one photon: N_s + 2 N_s n_telescopes.
inline constexpr std::uint64_t max_scatter_deviates = std::uint64_t{1} << 26;

// Rayleigh phase function (3 / 16 pi) (1 + cos^2 theta), normalised over 4 pi.
template <std::floating_point T> [[nodiscard]] T rayleigh_phase(T cos_theta) {
  return T{3.0 / (16.0 * std::numbers::pi)} * (T{1} + cos_theta * cos_theta);
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
// aerosol single-scattering albedo and asymmetry parameter, interpolated linearly in
// log10(altitude) and clamped to the table edges, as TransmissionTable does.
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

  // State at a wavelength (truncated to integer nm, clamped to the table) and altitude
  // [m a.s.l.] (finite and positive). Throws std::invalid_argument otherwise.
  template <std::floating_point T>
  [[nodiscard]] ScatteringState<T> state(double lambda_nm, T altitude_m) const {
    if (!std::isfinite(altitude_m) || !(altitude_m > T{0})) {
      throw std::invalid_argument("scattering table: altitude must be finite and positive");
    }
    const Row &row = rows_[row_index(lambda_nm)];
    const T log_h = std::log10(altitude_m);
    return {detail::rpol(log_h_m_, row.tau_mol.data(), log_h),
            detail::rpol(log_h_m_, row.tau_aer.data(), log_h),
            detail::rpol(log_h_m_, row.ssa.data(), log_h),
            detail::rpol(log_h_m_, row.asymmetry.data(), log_h)};
  }

private:
  struct Row {
    std::vector<double> tau_mol;
    std::vector<double> tau_aer;
    std::vector<double> ssa;
    std::vector<double> asymmetry;
  };

  int wl_first_{};
  std::vector<double> altitude_m_{};
  std::vector<double> log_h_m_{};
  std::vector<Row> rows_{};

  // Row of a wavelength: truncated to integer nm and clamped to the table, in floating point
  // (a cast of an out-of-range value to an integer type is undefined).
  [[nodiscard]] std::size_t row_index(double lambda_nm) const {
    if (!std::isfinite(lambda_nm)) {
      throw std::invalid_argument("scattering table: wavelength must be finite");
    }
    const auto first = static_cast<double>(wl_first_);
    const auto last = static_cast<double>(last_wavelength_nm());
    const double clamped = std::clamp(std::trunc(lambda_nm), first, last);
    return static_cast<std::size_t>(clamped - first);
  }

  void parse(std::string_view text) {
    struct Entry {
      long wavelength;
      double altitude_km, tau_mol, tau_aer, ssa, asymmetry;
    };
    std::vector<Entry> entries;
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
      if (tokens.size() != 6) {
        throw std::invalid_argument("scattering table: each row needs 6 values "
                                    "(wavelength altitude tau_mol tau_aer ssa g)");
      }
      const double wavelength = parse_number(tokens[0], "wavelength");
      if (wavelength != std::floor(wavelength) || wavelength < 1.0) {
        throw std::invalid_argument("scattering table: wavelengths must be integer nm >= 1");
      }
      entries.push_back({static_cast<long>(wavelength), parse_number(tokens[1], "altitude"),
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
    wl_first_ = static_cast<int>(first);
    altitude_m_.clear();
    log_h_m_.clear();
    for (const double h_km : altitudes) {
      altitude_m_.push_back(h_km * 1000.0);
      log_h_m_.push_back(std::log10(h_km * 1000.0));
    }
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
    const std::string text(token);
    char *end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    if (text.empty() || end != text.c_str() + text.size() || !std::isfinite(value)) {
      throw std::invalid_argument(std::string("scattering table: invalid ") + what + " '" + text +
                                  "'");
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

inline constexpr std::size_t n_scatter_columns = 10;

// Output columns of the scattered photons, in this order.
inline constexpr std::array<gstrickd::Column, n_scatter_columns> scatter_output_columns = {
    gstrickd::Column::x,  gstrickd::Column::y,         gstrickd::Column::z, gstrickd::Column::cx,
    gstrickd::Column::cy, gstrickd::Column::cz,        gstrickd::Column::t, gstrickd::Column::zem,
    gstrickd::Column::w,  gstrickd::Column::wavelength};

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
                            const gstrickd::PhotonChunk &chunk, std::size_t row,
                            std::uint64_t photon, std::uint64_t seed, std::uint32_t event,
                            std::span<double> boundary,
                            std::vector<std::vector<std::array<float, n_scatter_columns>>> &out) {
  using gstrickd::Column;
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
  const gstrickd::DrawKey key{seed, stage_atmosphere, event, 0, photon, stream_scattering};
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

      const double cos_theta = cx * nx + cy * ny - abs_cz * nz;
      const double bracket =
          mixed_phase<double>(cos_theta, tau_mol, tau_aer, local.ssa, local.asymmetry);
      const double scattered = weight * p_scatter * bracket * solid_angle * transmitted;
      if (!(scattered > 0.0)) {
        continue;
      }
      // Time at the scattering point; the stage adds the path from the point to the plane.
      const double time_at_point = t0 + s_path / speed_of_light_m_per_ns;
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
// Requires columns x y z cx cy t w wavelength (wavelength an integer nm within the scattering
// table, > 0: assign it before this step); cz and zem, if present, must agree with cx, cy and
// z; no p column (the acceptance is applied later). Only downward photons are scattered.
inline void scatter_chunk(const ScatteringConfig &config, const ScatteringTable &table,
                          const TransmissionTable &transmission,
                          std::span<const ScatterTelescope> telescopes,
                          const gstrickd::PhotonChunk &chunk, std::uint64_t seed,
                          std::uint32_t event, std::uint64_t first_index,
                          std::vector<std::vector<gstrickd::PhotonChunk>> &chunk_out) {
  using gstrickd::Column;
  if (!config.enabled) {
    return;
  }
  if (config.points_per_photon < 1 || config.points_per_photon > max_points_per_photon) {
    throw std::invalid_argument("scattering needs 1 to 1024 points per photon");
  }
  for (std::size_t i = 0; i < telescopes.size(); ++i) {
    const ScatterTelescope &tel = telescopes[i];
    if (!std::isfinite(tel.x_m) || !std::isfinite(tel.y_m) || !std::isfinite(tel.z_m) ||
        !(tel.radius_m > 0.0) || tel.z_m < 0.0) {
      throw std::invalid_argument("telescope needs finite x, y, z >= 0 (above the observation "
                                  "level) and radius > 0");
    }
    for (std::size_t k = 0; k < i; ++k) {
      if (telescopes[k].id == tel.id) {
        throw std::invalid_argument("telescope ids must be unique");
      }
    }
  }
  if (chunk.frame != gstrickd::Frame::sky) {
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
    const double r2 =
        static_cast<double>(cx[row]) * cx[row] + static_cast<double>(cy[row]) * cy[row];
    if (!(r2 <= 1.0)) {
      throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
    }
    const double cz = -std::sqrt(1.0 - r2);
    if (chunk.has(Column::cz) &&
        std::abs(static_cast<double>(chunk[Column::cz][row]) - cz) > 1e-5) {
      throw std::invalid_argument("cz column does not match cx, cy (cz = -sqrt(1-cx^2-cy^2))");
    }
    if (chunk.has(Column::zem) && std::abs(static_cast<double>(chunk[Column::zem][row]) -
                                           (observation_level + z[row])) > 1.0) {
      throw std::invalid_argument("zem column does not match z plus the observation level");
    }
    if (!(w[row] >= 0.0F) || !std::isfinite(w[row]) || !std::isfinite(z[row]) ||
        !std::isfinite(x[row]) || !std::isfinite(y[row]) || !std::isfinite(t[row])) {
      throw std::invalid_argument("weights, positions and times must be finite, weights >= 0");
    }
    const double lambda = wavelength[row];
    if (lambda != std::floor(lambda) || lambda < table.first_wavelength_nm() ||
        lambda > table.last_wavelength_nm()) {
      throw std::invalid_argument("wavelength must be an integer nm within the scattering table");
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
    mask |= gstrickd::column_bit(column);
  }
  for (std::size_t i = 0; i < telescopes.size(); ++i) {
    std::uint64_t offset = 0;
    for (const auto &existing : chunk_out[i]) {
      offset += existing.rows();
    }
    const auto &photons = out[i];
    for (std::size_t first = 0; first < photons.size(); first += gstrickd::max_chunk_rows) {
      const std::size_t n = std::min<std::size_t>(gstrickd::max_chunk_rows, photons.size() - first);
      gstrickd::PhotonChunk piece;
      piece.telescope_id = telescopes[i].id;
      piece.frame = gstrickd::Frame::sky;
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

} // namespace dlufd
