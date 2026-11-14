#pragma once

// Light sources: photons in frame sky (absolute emission point, direction) are propagated
// along their straight line to the horizontal plane through the telescope sphere centre,
// bent by the refraction tables if given, and written in frame G relative to the telescope.
// Per photon: no allocation and no I/O.

#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/timing.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace obdeect::propagation {

// Geometry of one telescope for sky-frame photons. Heights are m; z is measured above the
// observation level, as in frame sky.
struct SkyPlane {
  double telescope_x_m{};
  double telescope_y_m{};
  double telescope_z_m{};
  double observation_m{}; // observation level a.s.l.
  // Air delay along the straight path (both null: vacuum). Used when tables is null.
  const AtmosphereProfile *profile{};
  // Bending to the plane; the tables must be built for the plane altitude (observation_m +
  // telescope_z_m). Null: no bending.
  const RefractionTables *tables{};
};

// Converts a chunk in frame sky to frame G in place and returns the number of photons that do
// not reach the plane: photons moving upwards or horizontally, emitted below the plane, or not
// bendable. Their rows stay in the chunk with w = 0, cx = cy = 0, t = 0 and zem = observation
// level, so that row indices (and thus the random draws) do not change.
//
// Input: columns x y z cx cy t w (cz optional: if present its sign selects upwards photons).
// Output: columns x y cx cy t zem w (plus the optional columns that were present); z and cz
// are removed. zem is the emission altitude a.s.l.; t includes the vacuum path, the air delay
// (profile) or the bent travel time (tables).
inline std::uint64_t sky_to_ground(obdeect::stream::PhotonChunk &chunk, const SkyPlane &plane) {
  using obdeect::stream::Column;
  if (chunk.frame != obdeect::stream::Frame::sky) {
    throw std::invalid_argument("sky_to_ground needs photons in frame sky");
  }
  chunk.validate();
  if (!chunk.has(Column::z)) {
    throw std::invalid_argument("sky-frame photons need the z column");
  }
  const std::size_t rows = chunk.rows();
  const double plane_z = plane.telescope_z_m;
  const double plane_altitude = plane.observation_m + plane_z;
  if (plane.tables != nullptr &&
      std::abs(plane.tables->observation_level_cm() - 100.0 * plane_altitude) > 1e-6) {
    throw std::invalid_argument("refraction tables must be built for the telescope plane");
  }
  if (!chunk.has(Column::zem)) {
    chunk[Column::zem].assign(rows, 0.0F);
  }
  const bool has_cz = chunk.has(Column::cz);

  auto &x = chunk[Column::x];
  auto &y = chunk[Column::y];
  auto &z = chunk[Column::z];
  auto &cx = chunk[Column::cx];
  auto &cy = chunk[Column::cy];
  auto &t = chunk[Column::t];
  auto &zem = chunk[Column::zem];
  auto &w = chunk[Column::w];
  const auto &cz_column = chunk[Column::cz];

  std::uint64_t dropped = 0;
  for (std::size_t row = 0; row < rows; ++row) {
    const double dx = cx[row];
    const double dy = cy[row];
    const double r2 = dx * dx + dy * dy;
    if (!(r2 <= 1.0)) {
      throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
    }
    const double cz_magnitude = std::sqrt(1.0 - r2);
    const double emission_z = z[row];
    if (!std::isfinite(x[row]) || !std::isfinite(y[row]) || !std::isfinite(emission_z) ||
        !std::isfinite(t[row])) {
      throw std::invalid_argument("sky photon position and time must be finite");
    }
    bool reaches = cz_magnitude > 0.0 && emission_z >= plane_z;
    if (has_cz) {
      const double cz = cz_column[row];
      if (!std::isfinite(cz)) {
        throw std::invalid_argument("sky photon cz must be finite");
      }
      if (cz >= 0.0) {
        reaches = false; // moving upwards or horizontally
      } else if (std::abs(cz + cz_magnitude) > 1e-5) {
        throw std::invalid_argument("cz column does not match cx, cy (cz = -sqrt(1-cx^2-cy^2))");
      }
    }
    const double emission_altitude = plane.observation_m + emission_z;
    double x_out = 0.0;
    double y_out = 0.0;
    double cx_out = 0.0;
    double cy_out = 0.0;
    double t_out = 0.0;
    if (reaches) {
      const double path = (emission_z - plane_z) / cz_magnitude;
      x_out = static_cast<double>(x[row]) + path * dx - plane.telescope_x_m;
      y_out = static_cast<double>(y[row]) + path * dy - plane.telescope_y_m;
      cx_out = dx;
      cy_out = dy;
      if (plane.tables != nullptr) {
        double travel_ns = 0.0;
        reaches = bend(*plane.tables, emission_altitude, cx_out, cy_out, x_out, y_out, travel_ns);
        t_out = static_cast<double>(t[row]) + travel_ns;
      } else {
        t_out = static_cast<double>(t[row]) + path / speed_of_light_m_per_ns;
        if (plane.profile != nullptr) {
          const double zenith = std::atan2(std::sqrt(r2), cz_magnitude);
          t_out +=
              straight_path_air_delay_ns(*plane.profile, plane_altitude, emission_altitude, zenith);
        }
      }
    }
    if (!reaches) {
      ++dropped;
      x_out = 0.0;
      y_out = 0.0;
      cx_out = 0.0;
      cy_out = 0.0;
      t_out = 0.0;
      w[row] = 0.0F;
      zem[row] = static_cast<float>(plane.observation_m);
    } else {
      zem[row] = static_cast<float>(emission_altitude);
    }
    x[row] = static_cast<float>(x_out);
    y[row] = static_cast<float>(y_out);
    cx[row] = static_cast<float>(cx_out);
    cy[row] = static_cast<float>(cy_out);
    t[row] = static_cast<float>(t_out);
  }
  chunk[Column::z].clear();
  chunk[Column::cz].clear();
  chunk.frame = obdeect::stream::Frame::ground;
  return dropped;
}

} // namespace obdeect::propagation
