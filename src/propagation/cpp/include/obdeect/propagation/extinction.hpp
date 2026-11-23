#pragma once

// Atmospheric acceptance of photons in frame ground: airmass and transmission per photon.

#include "obdeect/propagation/chord.hpp"
#include "obdeect/propagation/sky.hpp"
#include "obdeect/propagation/transmission.hpp"

#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace obdeect::propagation {

// Spherical atmosphere seen from one telescope (see docs/extinction.md). The local vertical of
// the telescope is the unit vector n = (nx, ny, nz) in frame ground; earth_radius_m is the Earth
// radius R and telescope_altitude_m the altitude of the telescope [m a.s.l.]. The default n is
// the vertical of a flat ground frame.
struct CurvedPath {
  double earth_radius_m{};
  double telescope_altitude_m{};
  double nx{0.0};
  double ny{0.0};
  double nz{1.0};
};

// Throws unless R and the telescope altitude are finite (R > 0) and n is a unit vector.
inline void validate(const CurvedPath &path) {
  if (!std::isfinite(path.earth_radius_m) || !(path.earth_radius_m > 0.0)) {
    throw std::invalid_argument("earth radius must be finite and positive");
  }
  if (!std::isfinite(path.telescope_altitude_m)) {
    throw std::invalid_argument("telescope altitude must be finite");
  }
  if (!std::isfinite(path.nx) || !std::isfinite(path.ny) || !std::isfinite(path.nz) ||
      !(std::abs(path.nx * path.nx + path.ny * path.ny + path.nz * path.nz - 1.0) <= 1e-9)) {
    throw std::invalid_argument("local vertical must be a unit vector");
  }
}

// Plane-parallel airmass 1 / |cz| of a downward photon; unlimited_airmass for a horizontal one.
template <class T> [[nodiscard]] T airmass(T cx, T cy) {
  const T cz2 = T(1) - cx * cx - cy * cy;
  if (!std::isfinite(cz2) || cz2 < -static_cast<T>(obdeect::stream::direction_norm_tolerance)) {
    throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
  }
  return cz2 > T(0) ? T(1) / std::sqrt(cz2) : static_cast<T>(unlimited_airmass);
}

// p[row] = atmospheric transmission of every row; 0 for wavelengths at or above max_lambda.
// The telescope depths follow the table's wavelength mode (TransmissionTable::telescope_depths_at).
// Rows [begin, end) only (default: all rows).
template <class T>
void atmospheric_acceptance(const obdeect::stream::PhotonChunk &chunk,
                            const TransmissionTable &table, const TelescopeDepths &depths,
                            std::span<T> p, std::size_t begin = 0,
                            std::size_t end = std::numeric_limits<std::size_t>::max()) {
  using obdeect::stream::Column;
  const std::size_t rows = chunk.rows();
  if (p.size() != rows) {
    throw std::invalid_argument("acceptance vector must have one entry per row");
  }
  const auto &cx = chunk[Column::cx];
  const auto &cy = chunk[Column::cy];
  const auto &zem = chunk[Column::zem];
  const auto &wavelength = chunk[Column::wavelength];
  end = std::min(end, rows);
  for (std::size_t row = begin; row < end; ++row) {
    check_wavelength(static_cast<double>(wavelength[row]));
    const int iwl = table_wavelength(static_cast<double>(wavelength[row]));
    if (iwl < 0) {
      p[row] = T(0);
      continue;
    }
    const auto od = table.telescope_depths_at(depths, static_cast<double>(wavelength[row]));
    p[row] = table.transmission2<T>(static_cast<T>(wavelength[row]), static_cast<T>(zem[row]),
                                    airmass<T>(static_cast<T>(cx[row]), static_cast<T>(cy[row])),
                                    static_cast<T>(od.od_tel), static_cast<T>(od.od_focus));
  }
}

// Spherical-atmosphere variant of atmospheric_acceptance (frame ground). Per row the direction is
// d = (cx, cy, cz) with cz = -sqrt(1 - cx^2 - cy^2), and the zenith cosine relative to the local
// vertical n is cos_zenith = -(d . n). A photon moving away from the ground (cos_zenith clearly
// negative) gets p = 0; rounding-size negatives are taken as 0. An upward photon (cz > 0 in the
// cz column, which is then used) has cos_zenith = d . n with d = (cx, cy, cz), and its slant depth
// is the integral over [zem, h_t] (transmission_curved_upward). Wavelengths are handled as in
// atmospheric_acceptance. Per-row pure: the result does not depend on row order or chunking.
template <class T>
void atmospheric_acceptance_curved(const obdeect::stream::PhotonChunk &chunk,
                                   const TransmissionTable &table, const TelescopeDepths &depths,
                                   const CurvedPath &path, std::span<T> p, std::size_t begin = 0,
                                   std::size_t end = std::numeric_limits<std::size_t>::max()) {
  using obdeect::stream::Column;
  validate(path);
  const std::size_t rows = chunk.rows();
  if (p.size() != rows) {
    throw std::invalid_argument("acceptance vector must have one entry per row");
  }
  const auto &cx = chunk[Column::cx];
  const auto &cy = chunk[Column::cy];
  const auto &zem = chunk[Column::zem];
  const auto &wavelength = chunk[Column::wavelength];
  const bool has_cz = chunk.has(Column::cz);
  const std::vector<float> no_cz;
  const auto &cz_column = has_cz ? chunk[Column::cz] : no_cz;
  end = std::min(end, rows);
  for (std::size_t row = begin; row < end; ++row) {
    check_wavelength(static_cast<double>(wavelength[row]));
    const int iwl = table_wavelength(static_cast<double>(wavelength[row]));
    if (iwl < 0) {
      p[row] = T(0);
      continue;
    }
    const double x = static_cast<double>(cx[row]);
    const double y = static_cast<double>(cy[row]);
    const double cz2 = 1.0 - x * x - y * y;
    if (!std::isfinite(cz2) || cz2 < -obdeect::stream::direction_norm_tolerance) {
      throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
    }
    const bool upward = has_cz && cz_column[row] > 0.0F;
    double cos_zenith = 0.0;
    if (upward) {
      // d . n with d = (x, y, cz). A photon moving down against the local vertical (cos clearly
      // negative) cannot arrive from below and gets p = 0.
      cos_zenith = static_cast<double>(cz_column[row]) * path.nz + x * path.nx + y * path.ny;
      if (cos_zenith < -detail::curved_cos_tolerance) {
        p[row] = T(0);
        continue;
      }
    } else {
      // -(d . n) with d = (x, y, -sqrt(1 - x^2 - y^2)).
      cos_zenith = std::sqrt(std::max(cz2, 0.0)) * path.nz - x * path.nx - y * path.ny;
      if (cos_zenith < -detail::curved_cos_tolerance) {
        p[row] = T(0);
        continue;
      }
    }
    cos_zenith = std::clamp(cos_zenith, 0.0, 1.0);
    const auto od = table.telescope_depths_at(depths, static_cast<double>(wavelength[row]));
    if (upward) {
      p[row] = table.transmission_curved_upward<T>(
          static_cast<T>(wavelength[row]), static_cast<T>(zem[row]), static_cast<T>(cos_zenith),
          static_cast<T>(path.telescope_altitude_m), static_cast<T>(od.od_focus),
          path.earth_radius_m);
    } else {
      p[row] = table.transmission_curved<T>(static_cast<T>(wavelength[row]),
                                            static_cast<T>(zem[row]), static_cast<T>(cos_zenith),
                                            static_cast<T>(path.telescope_altitude_m),
                                            static_cast<T>(od.od_focus), path.earth_radius_m);
    }
  }
}

// Geometry of the sphere rows of one telescope (frame G, sphere representation): the sky frame
// with its observation level, the Earth radius (0: flat geometry) and the telescope centre C
// (x, y, z above the observation level). A row point Q = C + r has the sky position C + r.
struct SphereGeometry {
  double observation_m{};
  double earth_radius_m{}; // 0: flat
  double telescope_x_m{};
  double telescope_y_m{};
  double telescope_z_m{};
};

// Extinction of the sphere rows (columns x y z cx cy cz zem path wavelength; docs/propagation/
// extinction.md, "Sphere rows"). For a row with path L > 0 the emission point is E = Q - L d (a
// straight chord: for a refracted row the bent path is approximated by it), h_E = zem and h_Q the
// altitude of Q (flat: observation level + Q_z; spherical: |Q - Earth centre| - R). u is the cosine
// between d and the local vertical at Q (flat: cz). Then:
//   |u| >= near_horizontal_cos: the plane formula with h_Q in place of the telescope altitude: flat
//     exp(-|tau(h_E) - tau(h_Q)| / |u|); spherical: the slant integral between h_Q and h_E
//     (transmission_curved, or _upward if h_E < h_Q).
//   |u| < near_horizontal_cos: the chord integral of the extinction coefficient from E to Q
//     (chord.hpp; flat or spherical), which stays finite for u -> 0.
// Both forms get the in-telescope depth od_focus of the run's telescope depths, and the same
// limits as the plane rows (zero below the observation level, zero for tau(h_Q) above
// max_telescope_depth and for a total depth above max_total_depth). A row with L = 0 (emitted
// inside the sphere) has no atmospheric depth: p = exp(-od_focus) with the same limits.
// Spherical rows need a cz column (the sphere representation always has one). Per row pure:
// the result does not depend on row order or chunking.
template <class T>
void atmospheric_acceptance_sphere(const obdeect::stream::PhotonChunk &chunk,
                                   const TransmissionTable &table, const TelescopeDepths &depths,
                                   const SphereGeometry &geometry, std::span<T> p,
                                   std::size_t begin = 0,
                                   std::size_t end = std::numeric_limits<std::size_t>::max()) {
  using obdeect::stream::Column;
  if (!std::isfinite(geometry.observation_m) || !std::isfinite(geometry.earth_radius_m) ||
      geometry.earth_radius_m < 0.0 || !std::isfinite(geometry.telescope_x_m) ||
      !std::isfinite(geometry.telescope_y_m) || !std::isfinite(geometry.telescope_z_m)) {
    throw std::invalid_argument("sphere geometry must be finite (Earth radius not negative)");
  }
  const std::size_t rows = chunk.rows();
  if (p.size() != rows) {
    throw std::invalid_argument("acceptance vector must have one entry per row");
  }
  if (!chunk.has(Column::cz) || !chunk.has(Column::path)) {
    throw std::invalid_argument("sphere rows need the cz and path columns");
  }
  const auto &x = chunk[Column::x];
  const auto &y = chunk[Column::y];
  const auto &z = chunk[Column::z];
  const auto &cx = chunk[Column::cx];
  const auto &cy = chunk[Column::cy];
  const auto &cz_column = chunk[Column::cz];
  const auto &zem = chunk[Column::zem];
  const auto &path = chunk[Column::path];
  const auto &wavelength = chunk[Column::wavelength];
  const double obs = geometry.observation_m;
  const double radius = geometry.earth_radius_m;
  end = std::min(end, rows);
  for (std::size_t row = begin; row < end; ++row) {
    check_wavelength(static_cast<double>(wavelength[row]));
    const double lambda = static_cast<double>(wavelength[row]);
    if (table_wavelength(lambda) < 0) {
      p[row] = T(0);
      continue;
    }
    const double emission = static_cast<double>(zem[row]);
    const double length = static_cast<double>(path[row]);
    const double dx = static_cast<double>(cx[row]);
    const double dy = static_cast<double>(cy[row]);
    const double dz = static_cast<double>(cz_column[row]);
    if (!std::isfinite(emission) || !std::isfinite(length) || length < 0.0) {
      throw std::invalid_argument("sphere rows need finite zem and path >= 0");
    }
    if (!obdeect::stream::direction_in_unit_disc(dx, dy) ||
        !obdeect::stream::cz_is_consistent(dz, dx, dy)) {
      throw std::invalid_argument("sphere row direction must be a unit vector (cz matches cx, cy)");
    }
    if (emission < obs) {
      p[row] = T(0); // below the observation level (as transmission2)
      continue;
    }
    // Row point Q: sky position (X, Y, Z) with Z above the observation level, its altitude h_q,
    // and for spherical geometry the Earth-centred vector v_q = (X, Y, R + a) with a = obs + Z.
    const double X = geometry.telescope_x_m + static_cast<double>(x[row]);
    const double Y = geometry.telescope_y_m + static_cast<double>(y[row]);
    const double a = obs + geometry.telescope_z_m + static_cast<double>(z[row]);
    const double rho2 = X * X + Y * Y;
    double h_q = a;
    std::array<double, 3> v_q{X, Y, radius + a};
    std::array<double, 3> normal{0.0, 0.0, 1.0};
    if (radius > 0.0) {
      h_q = spherical_altitude_m(rho2, a, radius);
      const double norm = std::sqrt(rho2 + (radius + a) * (radius + a));
      normal = {X / norm, Y / norm, (radius + a) / norm};
    }
    const double od_focus = table.telescope_depths_at(depths, lambda).od_focus;
    const double tau_q = table.optical_depth<double>(lambda, h_q);
    if (tau_q > max_telescope_depth) {
      p[row] = T(0);
      continue;
    }
    if (length == 0.0) {
      // Emitted inside the sphere: no atmospheric depth.
      const double depth = od_focus;
      p[row] = depth > max_total_depth ? T(0) : static_cast<T>(std::exp(-depth));
      continue;
    }
    // u: cosine between d and the local vertical at Q (signed: positive upward).
    const double u = dx * normal[0] + dy * normal[1] + dz * normal[2];
    if (radius > 0.0) {
      // Zenith cosine of the ray at Q towards the ground (upward rays: along d).
      const double cos_q = dz > 0.0 ? u : -u;
      if (cos_q < -detail::curved_cos_tolerance) {
        p[row] = T(0); // moving away from the ground at Q: no path to the ground (as curved)
        continue;
      }
      if (std::abs(u) >= near_horizontal_cos) {
        const double cos_clamped = std::clamp(cos_q, 0.0, 1.0);
        p[row] =
            emission >= h_q
                ? table.transmission_curved<T>(static_cast<T>(lambda), static_cast<T>(emission),
                                               static_cast<T>(cos_clamped), static_cast<T>(h_q),
                                               static_cast<T>(od_focus), radius)
                : table.transmission_curved_upward<T>(
                      static_cast<T>(lambda), static_cast<T>(emission), static_cast<T>(cos_clamped),
                      static_cast<T>(h_q), static_cast<T>(od_focus), radius);
        continue;
      }
      const std::array<double, 3> v_e{v_q[0] - length * dx, v_q[1] - length * dy,
                                      v_q[2] - length * dz};
      const double chord = spherical_chord_depth(table, lambda, v_e, {dx, dy, dz}, length, radius);
      const double depth = chord + od_focus;
      p[row] = depth > max_total_depth ? T(0) : static_cast<T>(std::exp(-depth));
      continue;
    }
    if (std::abs(u) >= near_horizontal_cos) {
      p[row] = table.transmission2<T>(static_cast<T>(lambda), static_cast<T>(emission),
                                      static_cast<T>(1.0 / std::abs(u)), static_cast<T>(tau_q),
                                      static_cast<T>(od_focus));
      continue;
    }
    const double chord = flat_chord_depth(table, lambda, emission, h_q, length);
    const double depth = chord + od_focus;
    p[row] = depth > max_total_depth ? T(0) : static_cast<T>(std::exp(-depth));
  }
}

} // namespace obdeect::propagation
