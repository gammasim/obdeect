#pragma once

// Atmospheric acceptance of photons in frame ground: airmass and transmission per photon.

#include "dlufd/transmission.hpp"

#include "gstrickd/photons.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace dlufd {

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
void atmospheric_acceptance(const gstrickd::PhotonChunk &chunk, const TransmissionTable &table,
                            const TelescopeDepths &depths, std::span<T> p) {
  using gstrickd::Column;
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

} // namespace dlufd
