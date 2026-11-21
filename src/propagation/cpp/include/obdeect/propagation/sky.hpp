#pragma once

// Light sources: photons in frame sky (absolute emission point, direction) are propagated
// along their straight line to the horizontal plane through the telescope sphere centre,
// bent by the refraction tables if given, and written in frame G relative to the telescope.
// Per photon: no allocation and no I/O.

#include "obdeect/propagation/curved_refraction.hpp"
#include "obdeect/propagation/dispersion.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/refraction_horizon.hpp"
#include "obdeect/propagation/timing.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace obdeect::propagation {

// Altitude a.s.l. [m] of a point at horizontal distance sqrt(rho2) [m] from the array reference
// and at flat altitude altitude_m (observation level + z), seen from the Earth centre
// C = (0, 0, -(R + observation level)) with Earth radius R: |p - C| - R. Evaluated as
// (|p - C|^2 - R^2) / (|p - C| + R) so that there is no cancellation for large R.
inline double spherical_altitude_m(double rho2, double altitude_m, double earth_radius_m) {
  const double radius = earth_radius_m + altitude_m;
  const double distance = std::sqrt(rho2 + radius * radius);
  return (rho2 + altitude_m * (2.0 * earth_radius_m + altitude_m)) / (distance + earth_radius_m);
}

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
  // Earth radius R [m] for curved emission altitudes (spherical curvature), 0: flat. The emission
  // altitude is then |p - C| - R with C = (0, 0, -(R + observation_m)) and p = (x, y, z) of the
  // emission point. Then tables must be null; the bending and air delay come from curved_tables.
  double earth_radius_m{0.0};
  // Spherical bending and air delay (earth_radius_m > 0, profile set). The telescope altitude h_t
  // lies between the levels of the two tables: curved_tables at a_lo <= h_t and curved_tables_upper
  // at a_hi > h_t (both built for earth_radius_m). The integrals are interpolated linearly in h_t
  // with curved_weight = (h_t - a_lo) / (a_hi - a_lo). curved_tables_upper is null (and the weight
  // 0) if h_t lies on the lower level or no upper level fits inside the profile. Both are null if
  // profile is null.
  const CurvedRefractionTables *curved_tables{};
  const CurvedRefractionTables *curved_tables_upper{};
  double curved_weight{};
  // True: bend the photons (refraction not applied in the input). False: air delay only.
  bool curved_bending{};
  // Spherical altitude [m a.s.l.] of the telescope plane h_t (earth_radius_m > 0): r_t = R + h_t
  // gives the telescope level of the straight line, nu_t of the direction scale. Ignored in flat
  // geometry.
  double telescope_altitude_m{};
  // Dispersion (docs/propagation/refraction.md, "Dispersion"): reference wavelength [nm] of the
  // profile's n-1. 0: no dispersion. Otherwise each row's wavelength (column wavelength, assigned)
  // scales the bending and its travel time by k and k_group (refractivity_ratio and
  // group_refractivity_ratio). Needs bending (tables, or curved tables with curved_bending).
  double dispersion_reference_nm{0.0};
};

// Impact parameter b = |(P - C) x d| [m] of the straight line from the emission point P = (x, y, z)
// (run frame, z above the observation level) with direction d = (cx, cy, -cz_magnitude), about the
// Earth centre C = (0, 0, -(R + observation_m)). The line reaches the telescope level r_t = R + h_t
// if b < r_t, at zenith angle theta_s with sin(theta_s) = b / r_t.
inline double impact_parameter_m(double x, double y, double z, double observation_m,
                                 double earth_radius_m, double cx, double cy, double cz_magnitude) {
  const double vertical = earth_radius_m + observation_m + z; // (P - C) z-component
  const double u = y * cz_magnitude + vertical * cy;
  const double v = vertical * cx + x * cz_magnitude;
  const double w = x * cy - y * cx;
  return std::sqrt(u * u + v * v + w * w);
}

// Integrals of one curved table for a photon with impact parameter b [m] (about the Earth centre)
// emitted at emission_altitude_m. The table is looked up at its own level a: the zenith angle there
// is asin(b / (R + a)). False if the line does not reach that level (b >= R + a) or the lookup
// fails (emission below the level).
inline bool curved_table_integrals(const CurvedRefractionTables &table, double impact_m,
                                   double emission_altitude_m, double k, double k_group,
                                   CurvedRayIntegrals &out) {
  const double radius = table.earth_radius_m() + table.telescope_altitude_m();
  if (!(impact_m < radius)) {
    return false;
  }
  const double theta = std::atan2(impact_m, std::sqrt((radius - impact_m) * (radius + impact_m)));
  return table.lookup_scaled(theta, emission_altitude_m, k, k_group, out);
}

// Integrals of the photon at the exact telescope altitude by direct quadrature
// (curved_ray_integrals): the reference the tables approximate, used where no table is referenced
// to the telescope level (see curved_arrival). False if the line does not reach the telescope level
// or the emission is below it.
inline bool curved_direct_integrals(const SkyPlane &plane, double emission_altitude_m,
                                    double impact_m, double k, double k_group,
                                    CurvedRayIntegrals &out) {
  const double radius = plane.earth_radius_m + plane.telescope_altitude_m;
  if (!(impact_m < radius) || !(emission_altitude_m >= plane.telescope_altitude_m)) {
    return false;
  }
  const double theta = std::atan2(impact_m, std::sqrt((radius - impact_m) * (radius + impact_m)));
  try {
    out = curved_ray_integrals(*plane.profile, plane.earth_radius_m, plane.telescope_altitude_m,
                               theta, emission_altitude_m);
    scale_ray_integrals(out, k, k_group);
    const double start = CurvedRefractionTables::accurate_dispersion_blend_start_rad();
    const double end = CurvedRefractionTables::accurate_dispersion_blend_end_rad();
    if (plane.dispersion_reference_nm != 0.0 && theta > start) {
      const CurvedRayIntegrals exact =
          curved_ray_integrals(*plane.profile, plane.earth_radius_m, plane.telescope_altitude_m,
                               theta, emission_altitude_m, k, k_group);
      const double u = std::clamp((theta - start) / (end - start), 0.0, 1.0);
      const double weight = u * u * (3.0 - 2.0 * u);
      out.shift_m += weight * (exact.shift_m - out.shift_m);
      out.excess_bent_ns += weight * (exact.excess_bent_ns - out.excess_bent_ns);
      out.excess_straight_ns += weight * (exact.excess_straight_ns - out.excess_straight_ns);
    }
  } catch (const std::invalid_argument &) {
    return false;
  }
  return true;
}

// Spherical bending and air delay of one sky-frame photon with impact parameter impact_m that
// reaches the telescope level (plane.curved_tables required). Input: the straight-line position
// x_out, y_out (already relative to the telescope), direction cx, cy of the emission, the vacuum
// path [m] to the plane and the time t_in [ns]. Output: the bent position, direction and time.
// The integrals come from the table at the telescope level if there is one, else from the blend of
// the two tables around it. A table is referenced to its own level, so where the blend is not
// possible (emission between the telescope and the upper level, a grazing line that does not reach
// the lower level, no upper level inside the profile) the integrals are computed directly at the
// telescope altitude, which is exact and rare. Returns false if the photon is dropped (no
// integrals, or the bent direction leaves the unit disc). No allocation.
inline bool curved_arrival(const SkyPlane &plane, double emission_altitude_m, double impact_m,
                           double cx, double cy, double path_m, double t_in, double k,
                           double k_group, double &x_out, double &y_out, double &cx_out,
                           double &cy_out, double &t_out) {
  CurvedRayIntegrals v{};
  const double w = plane.curved_weight;
  const bool on_level = plane.curved_tables->telescope_altitude_m() == plane.telescope_altitude_m;
  bool have = false;
  if (on_level) {
    have =
        curved_table_integrals(*plane.curved_tables, impact_m, emission_altitude_m, k, k_group, v);
  } else {
    CurvedRayIntegrals lower{};
    CurvedRayIntegrals upper{};
    if (plane.curved_tables_upper != nullptr && w > 0.0 &&
        curved_table_integrals(*plane.curved_tables, impact_m, emission_altitude_m, k, k_group,
                               lower) &&
        curved_table_integrals(*plane.curved_tables_upper, impact_m, emission_altitude_m, k,
                               k_group, upper)) {
      v.shift_m = (1.0 - w) * lower.shift_m + w * upper.shift_m;
      v.excess_bent_ns = (1.0 - w) * lower.excess_bent_ns + w * upper.excess_bent_ns;
      v.excess_straight_ns = (1.0 - w) * lower.excess_straight_ns + w * upper.excess_straight_ns;
      have = true;
    } else {
      have = curved_direct_integrals(plane, emission_altitude_m, impact_m, k, k_group, v);
    }
  }
  if (!have) {
    return false;
  }
  const double vacuum_ns = path_m / speed_of_light_m_per_ns;
  if (!plane.curved_bending) {
    t_out = t_in + vacuum_ns + v.excess_straight_ns;
    return true;
  }
  // Lateral shift of the arrival point along the horizontal direction of motion.
  const double s = std::sqrt(cx * cx + cy * cy);
  if (s > 0.0) {
    x_out += v.shift_m * cx / s;
    y_out += v.shift_m * cy / s;
  }
  // Snell's law at the ends of the ray: s_obs / s_e = (1 + nu_e) / (1 + nu_t), with the true
  // telescope altitude for nu_t; the emission altitude is clamped to the profile top as in the
  // integrals.
  const double top = plane.profile->top_m();
  const double nu_e =
      k * plane.profile->index_minus_one_cm(100.0 * std::min(emission_altitude_m, top));
  const double nu_t = k * plane.profile->index_minus_one_cm(100.0 * plane.telescope_altitude_m);
  const double scale = (1.0 + nu_e) / (1.0 + nu_t);
  cx_out = cx * scale;
  cy_out = cy * scale;
  if (!(cx_out * cx_out + cy_out * cy_out < 1.0)) {
    return false;
  }
  t_out = t_in + vacuum_ns + v.excess_bent_ns;
  return true;
}

// Converts a chunk in frame sky to frame G in place and returns the number of photons that do
// not reach the plane: photons moving upwards or horizontally, emitted below the plane, or not
// bendable. Their rows stay in the chunk with w = 0, cx = cy = 0, t = 0 and zem = observation
// level, so that row indices (and thus the random draws) do not change.
//
// Input: columns x y z cx cy t w (cz optional: if present its sign selects upwards photons).
// Output: columns x y cx cy t zem w; z is removed. If preserve_cz is true, use and retain cz. zem
// is the emission altitude a.s.l.; t includes the vacuum path, the air delay (profile) or the bent
// travel time (tables).
inline std::uint64_t sky_to_ground(obdeect::stream::PhotonChunk &chunk, const SkyPlane &plane,
                                   bool preserve_cz = false) {
  using obdeect::stream::Column;
  if (chunk.frame != obdeect::stream::Frame::sky) {
    throw std::invalid_argument("sky_to_ground needs photons in frame sky");
  }
  chunk.validate();
  if (!chunk.has(Column::z)) {
    throw std::invalid_argument("sky-frame photons need the z column");
  }
  if (!std::isfinite(plane.earth_radius_m) || plane.earth_radius_m < 0.0) {
    throw std::invalid_argument("earth radius must be finite and not negative (0: flat)");
  }
  const bool curved = plane.earth_radius_m > 0.0;
  if (curved) {
    if (plane.tables != nullptr) {
      throw std::invalid_argument(
          "curved geometry takes curved tables, not flat refraction tables");
    }
    if (!std::isfinite(plane.telescope_altitude_m)) {
      throw std::invalid_argument("telescope altitude must be finite");
    }
    if ((plane.profile != nullptr) != (plane.curved_tables != nullptr) ||
        (plane.curved_tables == nullptr && plane.curved_tables_upper != nullptr)) {
      throw std::invalid_argument("curved bending needs a profile and its lower curved table");
    }
    if (plane.curved_tables != nullptr) {
      const double lower = plane.curved_tables->telescope_altitude_m();
      bool consistent = lower <= plane.telescope_altitude_m &&
                        plane.curved_tables->earth_radius_m() == plane.earth_radius_m;
      if (plane.curved_tables_upper != nullptr) {
        const double upper = plane.curved_tables_upper->telescope_altitude_m();
        consistent = consistent && upper > lower && plane.telescope_altitude_m < upper &&
                     plane.curved_tables_upper->earth_radius_m() == plane.earth_radius_m &&
                     std::abs(plane.curved_weight -
                              (plane.telescope_altitude_m - lower) / (upper - lower)) <= 1e-9;
      } else {
        consistent = consistent && plane.curved_weight == 0.0;
      }
      if (!consistent) {
        throw std::invalid_argument(
            "curved tables must bracket the telescope altitude, with its weight and radius");
      }
    }
  } else if (plane.curved_tables != nullptr || plane.curved_tables_upper != nullptr) {
    throw std::invalid_argument("curved tables need an Earth radius above 0");
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
  const bool dispersion = plane.dispersion_reference_nm != 0.0;
  if (dispersion) {
    if (!std::isfinite(plane.dispersion_reference_nm) || plane.dispersion_reference_nm < 0.0) {
      throw std::invalid_argument("dispersion reference wavelength must be finite and positive");
    }
    const bool bending =
        curved ? (plane.curved_bending && plane.curved_tables != nullptr) : plane.tables != nullptr;
    if (!bending) {
      throw std::invalid_argument("dispersion needs bending (a profile and photons not refracted)");
    }
    if (!chunk.has(Column::wavelength)) {
      throw std::invalid_argument("dispersion needs the wavelength column (assigned)");
    }
  }

  auto &x = chunk[Column::x];
  auto &y = chunk[Column::y];
  auto &z = chunk[Column::z];
  auto &cx = chunk[Column::cx];
  auto &cy = chunk[Column::cy];
  auto &t = chunk[Column::t];
  auto &zem = chunk[Column::zem];
  auto &w = chunk[Column::w];
  auto &cz_column = chunk[Column::cz];
  const std::vector<float> no_wavelength;
  const auto &wavelength = dispersion ? chunk[Column::wavelength] : no_wavelength;

  std::uint64_t dropped = 0;
  for (std::size_t row = 0; row < rows; ++row) {
    const double dx = cx[row];
    const double dy = cy[row];
    const double r2 = dx * dx + dy * dy;
    if (!obdeect::stream::direction_in_unit_disc(dx, dy)) {
      throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
    }
    double cz_magnitude = std::sqrt(std::max(0.0, 1.0 - r2));
    if (has_cz && cz_magnitude == 0.0 && cz_column[row] < 0.0) {
      // cx^2 + cy^2 rounded above 1: the validated cz column gives the magnitude.
      cz_magnitude = -cz_column[row];
    }
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
      } else if (!obdeect::stream::cz_matches_direction(cz, dx, dy)) {
        throw std::invalid_argument("cz column does not match cx, cy (cz = -sqrt(1-cx^2-cy^2))");
      }
      if (preserve_cz && cz < 0.0) {
        cz_magnitude = -cz;
        reaches = emission_z >= plane_z;
      }
    }
    double emission_altitude = plane.observation_m + emission_z;
    if (curved) {
      const double x_emit = x[row];
      const double y_emit = y[row];
      emission_altitude = spherical_altitude_m(x_emit * x_emit + y_emit * y_emit, emission_altitude,
                                               plane.earth_radius_m);
    }
    double x_out = 0.0;
    double y_out = 0.0;
    double cx_out = 0.0;
    double cy_out = 0.0;
    double t_out = 0.0;
    // Dispersion factors of this row (1 without dispersion; see SkyPlane::dispersion_reference_nm).
    double k = 1.0;
    double k_group = 1.0;
    if (dispersion) {
      const double lambda = static_cast<double>(wavelength[row]);
      if (lambda == 0.0) {
        throw std::invalid_argument("dispersion needs assigned wavelengths (none is 0 nm)");
      }
      k = refractivity_ratio(lambda, plane.dispersion_reference_nm);
      k_group = group_refractivity_ratio(lambda, plane.dispersion_reference_nm);
    }
    if (reaches) {
      const double path = (emission_z - plane_z) / cz_magnitude;
      x_out = static_cast<double>(x[row]) + path * dx - plane.telescope_x_m;
      y_out = static_cast<double>(y[row]) + path * dy - plane.telescope_y_m;
      cx_out = dx;
      cy_out = dy;
      if (curved) {
        // The straight line must reach the telescope level r_t = R + h_t (impact parameter b).
        const double r_t = plane.earth_radius_m + plane.telescope_altitude_m;
        const double b =
            impact_parameter_m(static_cast<double>(x[row]), static_cast<double>(y[row]), emission_z,
                               plane.observation_m, plane.earth_radius_m, dx, dy, cz_magnitude);
        reaches = b < r_t;
        if (reaches) {
          const double t_in = static_cast<double>(t[row]);
          if (plane.curved_tables != nullptr) {
            reaches = curved_arrival(plane, emission_altitude, b, dx, dy, path, t_in, k, k_group,
                                     x_out, y_out, cx_out, cy_out, t_out);
          } else {
            t_out = t_in + path / speed_of_light_m_per_ns;
          }
        }
      } else if (plane.tables != nullptr) {
        // Beyond the series limit of bend() (near the horizon) the ray is bent by quadrature
        // instead of being dropped. Both leave the outputs unchanged when they fail.
        double travel_ns = 0.0;
        reaches = bend(*plane.tables, emission_altitude, cx_out, cy_out, x_out, y_out, travel_ns, k,
                       k_group);
        if (!reaches && bend_series_limit_exceeded(*plane.tables, emission_altitude, dx, dy, k)) {
          reaches = bend_quadrature(*plane.tables, emission_altitude, cx_out, cy_out, x_out, y_out,
                                    travel_ns, k, k_group);
        }
        t_out = static_cast<double>(t[row]) + travel_ns;
      } else {
        t_out = static_cast<double>(t[row]) +
                path * (preserve_cz && has_cz ? std::sqrt(r2 + cz_magnitude * cz_magnitude) : 1.0) /
                    speed_of_light_m_per_ns;
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
    if (has_cz && preserve_cz) {
      if (!reaches) {
        cz_column[row] = -1.0F;
      } else if (plane.tables != nullptr || plane.curved_tables != nullptr) {
        cz_column[row] =
            static_cast<float>(-std::sqrt(std::max(0.0, 1.0 - cx_out * cx_out - cy_out * cy_out)));
      }
    }
  }
  chunk[Column::z].clear();
  if (!preserve_cz) {
    chunk[Column::cz].clear();
  }
  chunk.frame = obdeect::stream::Frame::ground;
  return dropped;
}

} // namespace obdeect::propagation
