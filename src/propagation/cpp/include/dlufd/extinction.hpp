#pragma once

// Atmospheric acceptance of photons in frame ground (G), as sim_telarray.c:2713-2862
// does for bunches without 3D data: plane-parallel airmass 1/|cz| and
// atmospheric_transmission2 at the bunch's integer wavelength.

#include "dlufd/transmission.hpp"

#include "gstrickd/photons.hpp"

#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>

namespace dlufd {

// Airmass 1/|cz| for direction cosines cx, cy of a downward photon, cz = -sqrt(1 - cx^2 -
// cy^2). Requires cx^2 + cy^2 <= 1. cz = 0 gives unlimited_airmass (sim_telarray.c:2707).
template <std::floating_point T> [[nodiscard]] T airmass(T cx, T cy) {
  const T cz = -std::sqrt(T{1} - cx * cx - cy * cy);
  return cz == T{0} ? T{unlimited_airmass} : T{1} / -cz;
}

// Atmospheric transmission p_atm of each row, without the detector envelope and without
// the fraction of a bunch (applied in thinning). Rows need a wavelength in [0, 1000) nm
// as an integer; rows at or above max_lambda_nm get 0, as sim_telarray skips them.
// Requires the columns cx, cy, zem and wavelength; out.size() >= rows.
template <std::floating_point T>
void atmospheric_acceptance(const gstrickd::PhotonChunk &chunk, const TransmissionTable &table,
                            const TelescopeDepths &depths, std::span<T> out) {
  const auto &cx = chunk[gstrickd::Column::cx];
  const auto &cy = chunk[gstrickd::Column::cy];
  const auto &zem = chunk[gstrickd::Column::zem];
  const auto &wavelength = chunk[gstrickd::Column::wavelength];
  for (std::size_t row = 0; row < chunk.rows(); ++row) {
    const auto iwl = static_cast<int>(wavelength[row]);
    if (iwl < 0 || iwl >= max_lambda_nm) {
      out[row] = T{0};
      continue;
    }
    const T a = airmass<T>(static_cast<T>(cx[row]), static_cast<T>(cy[row]));
    out[row] = table.transmission2<T>(static_cast<double>(iwl), static_cast<T>(zem[row]), a,
                                      static_cast<T>(depths.od_tel[iwl]),
                                      static_cast<T>(depths.od_focus[iwl]));
  }
}

} // namespace dlufd
