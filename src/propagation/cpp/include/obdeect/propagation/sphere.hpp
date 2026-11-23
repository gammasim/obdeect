#pragma once

// Sphere representation of frame G (docs/propagation/propagation.md, "Sphere representation"). A
// G row is the point where the photon's path inside the telescope sphere begins, relative to the
// sphere centre C = (telescope x, y, z) of the sky frame, with radius R (telescope_r_m):
//   emitted outside the sphere: the entry point (|r| = R);
//   emitted inside the sphere:  the emission point (|r| < R, path 0).
// Columns x y z cx cy cz t zem w path: direction at the row, time there, zem = emission altitude
// a.s.l., path = distance travelled from emission to the row point [m].
//
// sky_to_sphere() converts a chunk of sky-frame photons per row (see the header comment of
// sky_to_sphere). The rules, in order:
//   1. z < 0 (emission below the observation level): below_ground.
//   2. emission inside the sphere: the emission point, path 0, time and direction unchanged.
//   3. no refraction (no tables, no curved bending): the straight line, smallest s > 0 with
//      |P + s u - C| = R; the air delay of the straight path (flat or spherical).
//   4. refraction: the hybrid rule if the plane-mode computation reaches the telescope plane with
//      |cz| >= 0.1 at the crossing and |z_e - z_C| > R; else the 3D ray integrator (ray_march.hpp).
// Rows that do not reach the sphere are dropped as in plane mode (w = 0, zeros, cz = -1). No
// allocation and no I/O per photon.

#include "obdeect/propagation/curved_refraction.hpp"
#include "obdeect/propagation/ray_march.hpp"
#include "obdeect/propagation/sky.hpp"
#include "obdeect/propagation/timing.hpp"

#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace obdeect::propagation {

// Hybrid rule (step 4): the plane-mode crossing is used if |cz| at the crossing is at least this.
inline constexpr double hybrid_min_abs_cz = 0.1;

// Air delay [ns] of the straight path from P (sky frame, z above observation) along the unit
// direction u for length_m, in spherical geometry: the altitude along the line is taken from |r|
// (spherical_altitude_m). Gauss-Legendre on pieces of at most 250 m, as straight_line_air_delay_ns.
inline double straight_line_air_delay_sky_ns(const AtmosphereProfile &profile,
                                             const std::array<double, 3> &p,
                                             const std::array<double, 3> &u, double length_m,
                                             double earth_radius_m, double observation_m) {
  if (!(length_m > 0.0)) {
    return 0.0;
  }
  const RayMedium medium{&profile, 1.0, 1.0, observation_m, earth_radius_m};
  const double mean = detail::straight_path_mean(length_m, [&](double s) {
    const std::array<double, 3> x{p[0] + s * u[0], p[1] + s * u[1], p[2] + s * u[2]};
    return profile.refractive_index_minus_one(ray_altitude_m(medium, x));
  });
  return length_m * mean / speed_of_light_m_per_ns;
}

// Plane-mode arrival of one row (the body of sky_to_ground for one row, see sky.hpp). Inputs: the
// sky position (x, y, z), the direction (dx, dy) with |cz| = cz_magnitude and the length factor
// of the path, the emission altitude, the time t_in and the dispersion factors. Outputs: the bent
// position at the telescope plane relative to the telescope (x_out, y_out), the direction there
// and the time. Returns the drop reason (Drop::none if the plane is reached).
inline Drop plane_arrival(const SkyPlane &plane, double x, double y, double z, double dx, double dy,
                          double cz_magnitude, double length_factor, bool upward,
                          double emission_altitude, double t_in, double k, double k_group,
                          double &x_out, double &y_out, double &cx_out, double &cy_out,
                          double &t_out) {
  const bool curved = plane.earth_radius_m > 0.0;
  const double plane_z = plane.telescope_z_m;
  const double plane_altitude = plane.observation_m + plane_z;
  const double r2 = dx * dx + dy * dy;
  Drop reason = Drop::none;
  if (upward) {
    if (z < 0.0) {
      reason = Drop::below_ground;
    } else if (z > plane_z) {
      reason = Drop::moving_away;
    }
  } else if (z < plane_z) {
    reason = z < 0.0 ? Drop::below_ground : Drop::moving_away;
  }
  if (reason != Drop::none) {
    return reason;
  }
  const double path = (upward ? plane_z - z : z - plane_z) / cz_magnitude;
  x_out = x + path * dx - plane.telescope_x_m;
  y_out = y + path * dy - plane.telescope_y_m;
  cx_out = dx;
  cy_out = dy;
  const double zenith = std::atan2(std::sqrt(r2), cz_magnitude);
  t_out = 0.0;
  const double t_path = t_in;
  if (curved) {
    const double r_t = plane.earth_radius_m + plane.telescope_altitude_m;
    const double b = impact_parameter_m(x, y, z, plane.observation_m, plane.earth_radius_m, dx, dy,
                                        upward ? -cz_magnitude : cz_magnitude);
    if (!(b < r_t)) {
      return Drop::misses_level;
    }
    if (upward) {
      if (emission_altitude > plane.telescope_altitude_m) {
        return Drop::moving_away;
      }
      return upward_arrival(plane.profile, plane.curved_bending, emission_altitude,
                            plane.telescope_altitude_m, dx, dy, path, zenith, t_path, k, k_group,
                            x_out, y_out, cx_out, cy_out, t_out);
    }
    if (plane.curved_tables != nullptr) {
      return curved_arrival(plane, emission_altitude, b, dx, dy, path, t_path, k, k_group, x_out,
                            y_out, cx_out, cy_out, t_out)
                 ? Drop::none
                 : Drop::not_bendable;
    }
    t_out = t_path + path / speed_of_light_m_per_ns;
    return Drop::none;
  }
  const AtmosphereProfile *flat_profile =
      plane.profile != nullptr ? plane.profile
                               : (plane.tables != nullptr ? &plane.tables->profile() : nullptr);
  if (upward) {
    return upward_arrival(flat_profile, plane.tables != nullptr, emission_altitude, plane_altitude,
                          dx, dy, path * length_factor, zenith, t_path, k, k_group, x_out, y_out,
                          cx_out, cy_out, t_out);
  }
  if (plane.tables != nullptr) {
    double travel_ns = 0.0;
    bool bent =
        bend(*plane.tables, emission_altitude, cx_out, cy_out, x_out, y_out, travel_ns, k, k_group);
    if (!bent && bend_series_limit_exceeded(*plane.tables, emission_altitude, dx, dy, k)) {
      bent = bend_quadrature(*plane.tables, emission_altitude, cx_out, cy_out, x_out, y_out,
                             travel_ns, k, k_group);
    }
    t_out = t_path + travel_ns;
    return bent ? Drop::none : Drop::not_bendable;
  }
  t_out = t_path + path * length_factor / speed_of_light_m_per_ns;
  if (plane.profile != nullptr) {
    t_out += straight_path_air_delay_ns(*plane.profile, plane_altitude, emission_altitude, zenith);
  }
  return Drop::none;
}

// Whether the sphere rows of this plane use refraction (steps 3 and 4 of the header comment).
inline bool sphere_refracts(const SkyPlane &plane) {
  if (plane.earth_radius_m > 0.0) {
    return plane.curved_bending && plane.curved_tables != nullptr;
  }
  return plane.tables != nullptr;
}

// Integrator of one row (step 4e): from the sky position (x, y, z) with unit direction d, time
// t_in, to the sphere. Outputs: the entry point (absolute sky position), direction, time and path
// in state. Returns none, below_ground or misses_sphere.
inline Drop integrate_sphere_row(const SkyPlane &plane, double radius_m, const RayState &start,
                                 double k, double k_group, RayState &entry) {
  const RayMedium medium{plane.profile, k, k_group, plane.observation_m, plane.earth_radius_m};
  const std::array<double, 3> centre{plane.telescope_x_m, plane.telescope_y_m, plane.telescope_z_m};
  entry = start;
  switch (march_to_sphere(medium, entry, centre, radius_m)) {
  case MarchEnd::entered:
    return Drop::none;
  case MarchEnd::below_ground:
    return Drop::below_ground;
  case MarchEnd::escaped:
  case MarchEnd::closest_passed:
  case MarchEnd::step_limit:
    return Drop::misses_sphere;
  }
  return Drop::misses_sphere;
}

// Hybrid rule (step 4d). The plane-mode state at the telescope plane, Q (relative to the telescope,
// in the plane through C), with the bent direction (cx_out, cy_out, +-cz) and time t_cross, is the
// starting point of the 3D ray integrator (ray_march.hpp), not a straight tangent line: a straight
// line from Q is off by up to tens of metres for shallow rays (|cz| about 0.16, Q far from C).
// The straight line through Q decides the side of the entry: its smaller root s (entry at Q + s u,
// |Q + s u - C| = R; s < 0: before Q). Before Q (s < 0) the ray is integrated backwards from Q
// (reverse direction, time decreasing with n_g), otherwise forwards; the sphere is met with the
// straight finish of the last step (march_to_sphere / march_back_to_sphere). Output: the entry
// relative to C, the bent direction there, the time, and the path = the distance from the emission
// to the entry (the chord, as the straight rule). Misses: the straight line does not meet the
// sphere (misses_sphere), or the integration reaches the ground first (below_ground) or stops
// (misses_sphere). k and k_group: the refractivity ratios of the photon's wavelength (as
// integrate_sphere_row).
inline Drop hybrid_sphere_row(const SkyPlane &plane, double radius_m, double x_out, double y_out,
                              double cx_out, double cy_out, bool upward, double t_cross, double k,
                              double k_group, const std::array<double, 3> &emission,
                              std::array<double, 3> &entry_rel, std::array<double, 3> &u,
                              double &t_out, double &path_out) {
  const double cz_magnitude = std::sqrt(std::max(0.0, 1.0 - cx_out * cx_out - cy_out * cy_out));
  u = {cx_out, cy_out, upward ? cz_magnitude : -cz_magnitude};
  const std::array<double, 3> q{x_out, y_out, 0.0};
  const double qd = q[0] * u[0] + q[1] * u[1] + q[2] * u[2];
  const double q2 = q[0] * q[0] + q[1] * q[1];
  const double w2 = q2 - qd * qd; // |q x u|^2
  const double disc = radius_m * radius_m - w2;
  if (!(disc >= 0.0)) {
    return Drop::misses_sphere;
  }
  const double s = -qd - std::sqrt(disc);
  const std::array<double, 3> centre{plane.telescope_x_m, plane.telescope_y_m, plane.telescope_z_m};
  const RayMedium medium{plane.profile, k, k_group, plane.observation_m, plane.earth_radius_m};
  RayState state;
  state.x = {x_out + plane.telescope_x_m, y_out + plane.telescope_y_m, plane.telescope_z_m};
  state.d = u;
  state.t = t_cross;
  state.path = 0.0;
  const MarchEnd end = s < 0.0 ? march_back_to_sphere(medium, state, centre, radius_m)
                               : march_to_sphere(medium, state, centre, radius_m);
  if (end == MarchEnd::below_ground) {
    return Drop::below_ground;
  }
  if (end != MarchEnd::entered) {
    return Drop::misses_sphere;
  }
  entry_rel = {state.x[0] - centre[0], state.x[1] - centre[1], state.x[2] - centre[2]};
  u = state.d;
  t_out = state.t;
  const double dx = state.x[0] - emission[0];
  const double dy = state.x[1] - emission[1];
  const double dz = state.x[2] - emission[2];
  path_out = std::sqrt(dx * dx + dy * dy + dz * dz);
  return Drop::none;
}

// Converts a chunk of sky-frame photons to the sphere representation in place (frame G, columns
// x y z cx cy cz t zem w path; see the header comment). Input: x y z cx cy t w (cz optional, as
// sky_to_ground; use_cz_column as there: false for photons in the sky frame, true for scattered
// photons whose stored cz is the direction used). Dropped rows stay in the chunk with w = 0, zeros
// and cz = -1. Returns the drop reasons (DropCounts). Throws std::invalid_argument for invalid
// input, as sky_to_ground does; sphere_radius_m must be finite and positive.
inline DropCounts sky_to_sphere(obdeect::stream::PhotonChunk &chunk, const SkyPlane &plane,
                                double sphere_radius_m, bool use_cz_column = false) {
  using obdeect::stream::Column;
  if (chunk.frame != obdeect::stream::Frame::sky) {
    throw std::invalid_argument("sky_to_sphere needs photons in frame sky");
  }
  chunk.validate();
  if (!chunk.has(Column::z)) {
    throw std::invalid_argument("sky-frame photons need the z column");
  }
  if (!(std::isfinite(sphere_radius_m) && sphere_radius_m > 0.0)) {
    throw std::invalid_argument("telescope sphere radius must be finite and positive");
  }
  if (!std::isfinite(plane.earth_radius_m) || plane.earth_radius_m < 0.0) {
    throw std::invalid_argument("earth radius must be finite and not negative (0: flat)");
  }
  const bool curved = plane.earth_radius_m > 0.0;
  if (!curved && (plane.curved_tables != nullptr || plane.curved_tables_upper != nullptr)) {
    throw std::invalid_argument("curved tables need an Earth radius above 0");
  }
  if (curved && plane.tables != nullptr) {
    throw std::invalid_argument("curved geometry takes curved tables, not flat refraction tables");
  }
  if (!curved && plane.tables != nullptr &&
      std::abs(plane.tables->observation_level_cm() -
               100.0 * (plane.observation_m + plane.telescope_z_m)) > 1e-6) {
    throw std::invalid_argument("refraction tables must be built for the telescope plane");
  }
  if (sphere_refracts(plane) && plane.profile == nullptr) {
    throw std::invalid_argument("refraction needs an atmospheric profile");
  }
  const std::size_t rows = chunk.rows();
  if (!chunk.has(Column::zem)) {
    chunk[Column::zem].assign(rows, 0.0F);
  }
  const bool has_cz = chunk.has(Column::cz);
  const bool dispersion = plane.dispersion_reference_nm != 0.0;
  if (dispersion) {
    if (!std::isfinite(plane.dispersion_reference_nm) || plane.dispersion_reference_nm < 0.0) {
      throw std::invalid_argument("dispersion reference wavelength must be finite and positive");
    }
    if (!sphere_refracts(plane)) {
      throw std::invalid_argument("dispersion needs bending (a profile and photons not refracted)");
    }
    if (!chunk.has(Column::wavelength)) {
      throw std::invalid_argument("dispersion needs the wavelength column (assigned)");
    }
  }
  const std::array<double, 3> centre{plane.telescope_x_m, plane.telescope_y_m, plane.telescope_z_m};
  const double radius2 = sphere_radius_m * sphere_radius_m;
  const bool refracts = sphere_refracts(plane);

  auto &x = chunk[Column::x];
  auto &y = chunk[Column::y];
  auto &z = chunk[Column::z];
  auto &cx = chunk[Column::cx];
  auto &cy = chunk[Column::cy];
  auto &t = chunk[Column::t];
  auto &zem = chunk[Column::zem];
  auto &w = chunk[Column::w];
  auto &cz_column = chunk[Column::cz];
  auto &path_column = chunk[Column::path];
  path_column.assign(rows, 0.0F);
  cz_column.resize(rows, 0.0F);
  const std::vector<float> no_wavelength;
  const auto &wavelength = dispersion ? chunk[Column::wavelength] : no_wavelength;

  DropCounts dropped;
  for (std::size_t row = 0; row < rows; ++row) {
    const double dx = cx[row];
    const double dy = cy[row];
    const double r2 = dx * dx + dy * dy;
    if (!obdeect::stream::direction_in_unit_disc(dx, dy)) {
      throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
    }
    const double px = x[row], py = y[row], pz = z[row];
    if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) || !std::isfinite(t[row])) {
      throw std::invalid_argument("photon position and time in the sky frame must be finite");
    }
    double cz_magnitude = std::sqrt(std::max(0.0, 1.0 - r2));
    double length_factor = 1.0;
    double cz_signed = -cz_magnitude;
    if (has_cz) {
      const double cz = cz_column[row];
      if (!std::isfinite(cz)) {
        throw std::invalid_argument("cz of a photon in the sky frame must be finite");
      }
      if (!obdeect::stream::cz_is_consistent(cz, dx, dy)) {
        throw std::invalid_argument("cz column does not match cx, cy (cz = -sqrt(1-cx^2-cy^2))");
      }
      if (use_cz_column) {
        cz_magnitude = std::fabs(cz);
        length_factor = std::sqrt(r2 + cz_magnitude * cz_magnitude);
      } else if (cz_magnitude == 0.0) {
        cz_magnitude = std::fabs(cz);
      }
      cz_signed = cz;
    }
    // One direction for the geometry (as the plane mode): |cz| from cx, cy (the column only if
    // use_cz_column, or if cx^2 + cy^2 rounded above 1), so that the straight line, the plane
    // crossing and the integrator start agree to the last bit of cx and cy.
    const double sign = cz_signed > 0.0 ? 1.0 : (cz_signed < 0.0 ? -1.0 : 0.0);
    const double cz_geometry = sign * cz_magnitude;
    const bool upward = cz_geometry > 0.0;
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
    // Emission altitude a.s.l. (spherical geometry: from |r| about the Earth centre).
    const double emission_altitude =
        curved ? spherical_altitude_m(px * px + py * py, plane.observation_m + pz,
                                      plane.earth_radius_m)
               : plane.observation_m + pz;
    // Direction used for the straight geometry (unit: the path is a distance).
    const double norm = std::sqrt(r2 + cz_geometry * cz_geometry);
    const std::array<double, 3> u{dx / norm, dy / norm, cz_geometry / norm};
    const std::array<double, 3> p{px, py, pz};
    const std::array<double, 3> p_rel{px - centre[0], py - centre[1], pz - centre[2]};

    double cx_out = dx, cy_out = dy, cz_out = cz_signed;
    double t_out = static_cast<double>(t[row]);
    double path_out = 0.0;
    std::array<double, 3> row_rel{0.0, 0.0, 0.0};
    Drop reason = Drop::none;
    const double q2 = p_rel[0] * p_rel[0] + p_rel[1] * p_rel[1] + p_rel[2] * p_rel[2];

    if (pz < 0.0) {
      reason = Drop::below_ground; // the ground is opaque (checked before the sphere)
    } else if (q2 <= radius2) {
      // Step 2: emitted inside the sphere.
      row_rel = p_rel;
    } else if (!refracts) {
      // Step 3: straight line.
      const double qd = p_rel[0] * u[0] + p_rel[1] * u[1] + p_rel[2] * u[2];
      const double cxp = p_rel[1] * u[2] - p_rel[2] * u[1];
      const double cyp = p_rel[2] * u[0] - p_rel[0] * u[2];
      const double czp = p_rel[0] * u[1] - p_rel[1] * u[0];
      const double disc = radius2 - (cxp * cxp + cyp * cyp + czp * czp);
      const double s = -qd - std::sqrt(std::max(0.0, disc));
      if (!(disc >= 0.0) || !(s > 0.0)) {
        reason = Drop::misses_sphere;
      } else if (pz + s * u[2] < 0.0) {
        reason = Drop::below_ground;
      } else {
        row_rel = {p_rel[0] + s * u[0], p_rel[1] + s * u[1], p_rel[2] + s * u[2]};
        path_out = s;
        double delay = 0.0;
        if (plane.profile != nullptr) {
          if (curved) {
            delay = straight_line_air_delay_sky_ns(*plane.profile, p, u, s, plane.earth_radius_m,
                                                   plane.observation_m);
          } else {
            const double h_entry = plane.observation_m + pz + s * u[2];
            delay = straight_line_air_delay_ns(*plane.profile, plane.observation_m + pz, h_entry, s,
                                               1.0);
          }
        }
        t_out = static_cast<double>(t[row]) + s / speed_of_light_m_per_ns + delay;
      }
    } else {
      // Step 4: refraction. The plane-mode computation decides between the hybrid rule and the
      // integrator.
      double x_p = 0.0, y_p = 0.0, cx_p = 0.0, cy_p = 0.0, t_p = 0.0;
      const Drop plane_reason = plane_arrival(
          plane, px, py, pz, dx, dy, cz_magnitude, length_factor, upward, emission_altitude,
          static_cast<double>(t[row]), k, k_group, x_p, y_p, cx_p, cy_p, t_p);
      const double cz_at_plane = std::sqrt(std::max(0.0, 1.0 - cx_p * cx_p - cy_p * cy_p));
      if (plane_reason == Drop::none && cz_at_plane >= hybrid_min_abs_cz &&
          std::abs(pz - plane.telescope_z_m) > sphere_radius_m) {
        std::array<double, 3> u_t{};
        reason = hybrid_sphere_row(plane, sphere_radius_m, x_p, y_p, cx_p, cy_p, upward, t_p, k,
                                   k_group, p, row_rel, u_t, t_out, path_out);
        cx_out = u_t[0];
        cy_out = u_t[1];
        cz_out = u_t[2];
      } else {
        RayState start;
        start.x = p;
        start.d = u;
        start.t = static_cast<double>(t[row]);
        start.path = 0.0;
        RayState entry;
        reason = integrate_sphere_row(plane, sphere_radius_m, start, k, k_group, entry);
        if (reason == Drop::none) {
          row_rel = {entry.x[0] - centre[0], entry.x[1] - centre[1], entry.x[2] - centre[2]};
          cx_out = entry.d[0];
          cy_out = entry.d[1];
          cz_out = entry.d[2];
          t_out = entry.t;
          path_out = entry.path;
        }
      }
    }

    if (reason != Drop::none) {
      dropped.count(reason);
      x[row] = 0.0F;
      y[row] = 0.0F;
      z[row] = 0.0F;
      cx[row] = 0.0F;
      cy[row] = 0.0F;
      cz_column[row] = -1.0F;
      t[row] = 0.0F;
      w[row] = 0.0F;
      zem[row] = static_cast<float>(plane.observation_m);
      path_column[row] = 0.0F;
      continue;
    }
    x[row] = static_cast<float>(row_rel[0]);
    y[row] = static_cast<float>(row_rel[1]);
    z[row] = static_cast<float>(row_rel[2]);
    cx[row] = static_cast<float>(cx_out);
    cy[row] = static_cast<float>(cy_out);
    cz_column[row] = static_cast<float>(cz_out);
    t[row] = static_cast<float>(t_out);
    zem[row] = static_cast<float>(emission_altitude);
    path_column[row] = static_cast<float>(path_out);
  }
  chunk.frame = obdeect::stream::Frame::ground;
  return dropped;
}

} // namespace obdeect::propagation
