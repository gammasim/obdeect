#pragma once

// Atmospheric acceptance of photons in frame ground: airmass and transmission per photon.

#include "obdeect/propagation/transmission.hpp"

#include "obdeect/stream/photons.hpp"

#include <algorithm>
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
  if (!std::isfinite(cz2) || cz2 < -T(4) * std::numeric_limits<T>::epsilon()) {
    throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
  }
  return cz2 > T(0) ? T(1) / std::sqrt(cz2) : static_cast<T>(unlimited_airmass);
}

// p[row] = atmospheric transmission of every row; 0 for wavelengths at or above max_lambda.
template <class T>
void atmospheric_acceptance(const obdeect::stream::PhotonChunk &chunk,
                            const TransmissionTable &table, const TelescopeDepths &depths,
                            std::span<T> p) {
  using obdeect::stream::Column;
  const std::size_t rows = chunk.rows();
  if (p.size() != rows) {
    throw std::invalid_argument("acceptance vector must have one entry per row");
  }
  const auto &cx = chunk[Column::cx];
  const auto &cy = chunk[Column::cy];
  const auto &zem = chunk[Column::zem];
  const auto &wavelength = chunk[Column::wavelength];
  for (std::size_t row = 0; row < rows; ++row) {
    check_wavelength(static_cast<double>(wavelength[row]));
    const int iwl = table_wavelength(static_cast<double>(wavelength[row]));
    if (iwl < 0) {
      p[row] = T(0);
      continue;
    }
    const auto index = static_cast<std::size_t>(iwl);
    p[row] = table.transmission2<T>(static_cast<T>(wavelength[row]), static_cast<T>(zem[row]),
                                    airmass<T>(static_cast<T>(cx[row]), static_cast<T>(cy[row])),
                                    static_cast<T>(depths.od_tel[index]),
                                    static_cast<T>(depths.od_focus[index]));
  }
}

// Spherical-atmosphere variant of atmospheric_acceptance (frame ground). Per row the direction is
// d = (cx, cy, cz) with cz = -sqrt(1 - cx^2 - cy^2), and the zenith cosine relative to the local
// vertical n is cos_zenith = -(d . n). A photon moving away from the ground (cos_zenith clearly
// negative) gets p = 0; rounding-size negatives are taken as 0. Wavelengths are handled as in
// atmospheric_acceptance. Per-row pure: the result does not depend on row order or chunking.
template <class T>
void atmospheric_acceptance_curved(const obdeect::stream::PhotonChunk &chunk,
                                   const TransmissionTable &table, const TelescopeDepths &depths,
                                   const CurvedPath &path, std::span<T> p) {
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
  const double eps = std::numeric_limits<double>::epsilon();
  for (std::size_t row = 0; row < rows; ++row) {
    check_wavelength(static_cast<double>(wavelength[row]));
    const int iwl = table_wavelength(static_cast<double>(wavelength[row]));
    if (iwl < 0) {
      p[row] = T(0);
      continue;
    }
    const double x = static_cast<double>(cx[row]);
    const double y = static_cast<double>(cy[row]);
    const double cz2 = 1.0 - x * x - y * y;
    if (!std::isfinite(cz2) || cz2 < -4.0 * eps) {
      throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
    }
    // -(d . n) with d = (x, y, -sqrt(1 - x^2 - y^2)).
    double cos_zenith = std::sqrt(std::max(cz2, 0.0)) * path.nz - x * path.nx - y * path.ny;
    if (cos_zenith < -detail::curved_cos_tolerance) {
      p[row] = T(0);
      continue;
    }
    cos_zenith = std::clamp(cos_zenith, 0.0, 1.0);
    const auto index = static_cast<std::size_t>(iwl);
    p[row] = table.transmission_curved<T>(
        static_cast<T>(wavelength[row]), static_cast<T>(zem[row]), static_cast<T>(cos_zenith),
        static_cast<T>(path.telescope_altitude_m), static_cast<T>(depths.od_focus[index]),
        path.earth_radius_m);
  }
}

} // namespace obdeect::propagation
