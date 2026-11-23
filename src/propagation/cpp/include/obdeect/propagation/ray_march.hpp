#pragma once

// Three-dimensional ray integration for the sphere representation (docs/propagation/refraction.md,
// "Integrator"). A ray with refractive index n = 1 + k nu(h) and group index n_g = 1 + k_group
// nu(h) follows
//   dX/ds = d,   dd/ds = (g - (g . d) d) / n,   dT/ds = n_g / c,   g = grad n = k nu'(h) up,
// with up = (0, 0, 1) in flat geometry and the radial unit vector from the Earth centre in
// spherical geometry. The frame is frame sky (z above the observation level). Arithmetic is in
// double; no allocation, no I/O. Inputs and outputs may be float or double.

#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/sky.hpp"
#include "obdeect/propagation/timing.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace obdeect::propagation {

// Medium of a ray: profile (null: vacuum), k and k_group as in bend(), observation level [m a.s.l.]
// and Earth radius (0: flat geometry).
struct RayMedium {
  const AtmosphereProfile *profile{};
  double k{1.0};
  double k_group{1.0};
  double observation_m{};
  double earth_radius_m{0.0};
};

// State of a ray: position (sky frame, z above the observation level) [m], unit direction, time
// [ns] and the path length [m] travelled from the start.
struct RayState {
  std::array<double, 3> x{};
  std::array<double, 3> d{};
  double t{};
  double path{};
};

// Altitude a.s.l. [m] of a sky-frame position.
inline double ray_altitude_m(const RayMedium &medium, const std::array<double, 3> &x) {
  const double h = medium.observation_m + x[2];
  if (medium.earth_radius_m > 0.0) {
    return spherical_altitude_m(x[0] * x[0] + x[1] * x[1], h, medium.earth_radius_m);
  }
  return h;
}

// Right-hand side of the ray equations at position x, direction d: (dX/ds, dd/ds, dt/ds) in one
// array of 7 values (dX, dd, dT).
inline std::array<double, 7> ray_rhs(const RayMedium &medium, const std::array<double, 3> &x,
                                     const std::array<double, 3> &d) {
  const double altitude = ray_altitude_m(medium, x);
  double nu = 0.0;
  double slope_cm = 0.0; // d(n-1)/dh [1/cm]
  if (medium.profile != nullptr) {
    medium.profile->index_and_slope_cm(100.0 * altitude, nu, slope_cm);
  }
  const double n = 1.0 + medium.k * nu;
  const double n_group = 1.0 + medium.k_group * nu;
  const double slope = 100.0 * slope_cm * medium.k; // d n / dh [1/m]
  std::array<double, 3> up{0.0, 0.0, 1.0};
  if (medium.earth_radius_m > 0.0) {
    const double ux = x[0], uy = x[1], uz = x[2] + medium.earth_radius_m + medium.observation_m;
    const double r = std::sqrt(ux * ux + uy * uy + uz * uz);
    up = {ux / r, uy / r, uz / r};
  }
  const std::array<double, 3> g{slope * up[0], slope * up[1], slope * up[2]};
  const double g_dot_d = g[0] * d[0] + g[1] * d[1] + g[2] * d[2];
  std::array<double, 7> out{};
  for (std::size_t i = 0; i < 3; ++i) {
    out[i] = d[i];
    out[3 + i] = (g[i] - g_dot_d * d[i]) / n;
  }
  out[6] = n_group / speed_of_light_m_per_ns;
  return out;
}

// One classical RK4 step of length h [m] along the ray; the direction is renormalised afterwards.
inline void ray_rk4_step(const RayMedium &medium, RayState &state, double h) {
  using Vec = std::array<double, 3>;
  const auto advance = [](const RayState &base, const std::array<double, 7> &k, double f) {
    RayState s = base;
    for (std::size_t i = 0; i < 3; ++i) {
      s.x[i] = base.x[i] + f * k[i];
      s.d[i] = base.d[i] + f * k[3 + i];
    }
    s.t = base.t + f * k[6];
    return s;
  };
  const std::array<double, 7> k1 = ray_rhs(medium, state.x, state.d);
  const RayState s2 = advance(state, k1, 0.5 * h);
  const std::array<double, 7> k2 = ray_rhs(medium, s2.x, s2.d);
  const RayState s3 = advance(state, k2, 0.5 * h);
  const std::array<double, 7> k3 = ray_rhs(medium, s3.x, s3.d);
  const RayState s4 = advance(state, k3, h);
  const std::array<double, 7> k4 = ray_rhs(medium, s4.x, s4.d);
  Vec x{};
  Vec d{};
  for (std::size_t i = 0; i < 3; ++i) {
    x[i] = state.x[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    d[i] = state.d[i] + h / 6.0 * (k1[3 + i] + 2.0 * k2[3 + i] + 2.0 * k3[3 + i] + k4[3 + i]);
  }
  const double norm = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  state.x = x;
  state.d = {d[0] / norm, d[1] / norm, d[2] / norm};
  state.t += h / 6.0 * (k1[6] + 2.0 * k2[6] + 2.0 * k3[6] + k4[6]);
  state.path += h;
}

// Intersection of the straight segment X + s d, 0 <= s <= length (d unit), with the sphere of
// centre C and radius R. On a hit returns true and the smallest such s.
inline bool segment_sphere_entry(const std::array<double, 3> &x, const std::array<double, 3> &d,
                                 const std::array<double, 3> &centre, double radius_m,
                                 double length_m, double &s) {
  const double qx = x[0] - centre[0], qy = x[1] - centre[1], qz = x[2] - centre[2];
  const double b = qx * d[0] + qy * d[1] + qz * d[2];
  const double c = qx * qx + qy * qy + qz * qz - radius_m * radius_m;
  if (c <= 0.0) {
    s = 0.0;
    return true;
  }
  // disc = R^2 - |q x d|^2 (no cancellation for distant starts)
  const double cx = qy * d[2] - qz * d[1], cy = qz * d[0] - qx * d[2], cz = qx * d[1] - qy * d[0];
  const double disc = radius_m * radius_m - (cx * cx + cy * cy + cz * cz);
  if (!(disc >= 0.0)) {
    return false;
  }
  const double root = -b - std::sqrt(disc);
  if (root >= 0.0 && root <= length_m) {
    s = root;
    return true;
  }
  return false;
}

// Moves the state onto the sphere (centre C, radius R) along the straight line through its point
// and direction d: the displacement along d is -q.d - sqrt(R^2 - |q x d|^2) with q = x - C (the
// entry of the line, the smaller root of |q + s d|^2 = R^2; a backward march takes the same root
// as its entry). Time and path follow the displacement at the state's group index. The straight
// line is accurate only for a short displacement: a curved ray leaves it by about |kappa| s^2 / 2
// (kappa the curvature of the ray). The callers take an RK4 step onto the estimated crossing first,
// so the displacement here is a fraction of a metre.
inline void straight_to_sphere(const RayMedium &medium, RayState &state,
                               const std::array<double, 3> &centre, double radius_m) {
  const std::array<double, 3> &d = state.d;
  const double qx = state.x[0] - centre[0], qy = state.x[1] - centre[1],
               qz = state.x[2] - centre[2];
  const double qd = qx * d[0] + qy * d[1] + qz * d[2];
  const double cx = qy * d[2] - qz * d[1], cy = qz * d[0] - qx * d[2], cz = qx * d[1] - qy * d[0];
  const double disc = radius_m * radius_m - (cx * cx + cy * cy + cz * cz);
  const double displacement = -qd - std::sqrt(std::max(0.0, disc));
  const double altitude = ray_altitude_m(medium, state.x);
  const double ng =
      1.0 +
      medium.k_group *
          (medium.profile != nullptr ? medium.profile->refractive_index_minus_one(altitude) : 0.0);
  for (std::size_t i = 0; i < 3; ++i) {
    state.x[i] += displacement * d[i];
  }
  state.t += displacement * ng / speed_of_light_m_per_ns;
  state.path += displacement;
}

// Result of integrating a ray to the telescope sphere.
enum class MarchEnd : std::uint8_t { entered, below_ground, escaped, closest_passed, step_limit };

// Integrates the ray from state (sky frame) until it enters the sphere (centre C, radius R). On
// entry, state holds the entry point (absolute sky position), its direction, time and path.
// Stops with the reason otherwise (state is then the last point reached). Each step is step_m
// long at most; the sphere is checked on the straight segment of each step, so a step that enters
// the sphere is taken to the straight estimate of the crossing by RK4 and finished with the
// straight residual (straight_to_sphere). max_steps bounds the work.
inline MarchEnd march_to_sphere(const RayMedium &medium, RayState &state,
                                const std::array<double, 3> &centre, double radius_m,
                                double step_m = 50.0, std::size_t max_steps = 200000) {
  const double top_m = medium.profile != nullptr ? medium.profile->top_m() : 0.0;
  for (std::size_t steps = 0;; ++steps) {
    const std::array<double, 3> q{state.x[0] - centre[0], state.x[1] - centre[1],
                                  state.x[2] - centre[2]};
    const double q2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
    if (q2 <= radius_m * radius_m) {
      return MarchEnd::entered; // reached by a step that missed the straight segment
    }
    const double altitude = ray_altitude_m(medium, state.x);
    if (altitude < medium.observation_m) {
      return MarchEnd::below_ground;
    }
    if (medium.profile != nullptr && altitude > top_m && state.d[2] > 0.0) {
      return MarchEnd::escaped;
    }
    const double away = q[0] * state.d[0] + q[1] * state.d[1] + q[2] * state.d[2];
    if (away > 0.0) {
      return MarchEnd::closest_passed;
    }
    if (steps >= max_steps) {
      return MarchEnd::step_limit;
    }
    double s_hit = 0.0;
    if (segment_sphere_entry(state.x, state.d, centre, radius_m, step_m, s_hit)) {
      // RK4 onto the straight estimate of the crossing, then the straight residual.
      ray_rk4_step(medium, state, s_hit);
      straight_to_sphere(medium, state, centre, radius_m);
      return MarchEnd::entered;
    }
    ray_rk4_step(medium, state, step_m);
  }
}

// Integrates the ray backwards from state (sky frame, forward direction d as in RayState) until it
// meets the sphere (centre C, radius R): the first point on the backward ray at which |x - C| = R,
// from either side of the sphere (the state may be inside it). Backward steps of -step_m; each step
// is checked on its straight segment for the two roots of |x - sigma d - C| = R, 0 < sigma <= step.
// On a hit, state holds the crossing (absolute sky position), the forward direction there, the time
// (decreased along the ray: dT/ds = n_g / c with s running backwards) and the path (decreased).
// Stops with the reason otherwise: below_ground, escaped (above the table top moving upwards) or
// step_limit. The crossing of a step that reaches it is the RK4 step onto the straight estimate
// plus the straight residual (straight_to_sphere).
inline MarchEnd march_back_to_sphere(const RayMedium &medium, RayState &state,
                                     const std::array<double, 3> &centre, double radius_m,
                                     double step_m = 50.0, std::size_t max_steps = 200000) {
  const double top_m = medium.profile != nullptr ? medium.profile->top_m() : 0.0;
  for (std::size_t steps = 0;; ++steps) {
    const std::array<double, 3> q{state.x[0] - centre[0], state.x[1] - centre[1],
                                  state.x[2] - centre[2]};
    const double q2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
    const double altitude = ray_altitude_m(medium, state.x);
    if (altitude < medium.observation_m) {
      return MarchEnd::below_ground;
    }
    if (medium.profile != nullptr && altitude > top_m && state.d[2] < 0.0) {
      return MarchEnd::escaped;
    }
    if (steps >= max_steps) {
      return MarchEnd::step_limit;
    }
    // Roots sigma of |q - sigma d|^2 = R^2 (distance along -d). The point q - sigma d is an entry
    // of the forward ray when its radial direction is inward for d: (q - sigma d) . d < 0, i.e.
    // sigma > q . d. That is the larger root only (the smaller one is the exit of the forward ray).
    const double qd = q[0] * state.d[0] + q[1] * state.d[1] + q[2] * state.d[2];
    const double disc = qd * qd - (q2 - radius_m * radius_m);
    if (disc >= 0.0) {
      const double sigma = qd + std::sqrt(disc);
      if (sigma > 0.0 && sigma <= step_m) {
        // RK4 back onto the straight estimate of the crossing, then the straight residual (the
        // same root: the entry is the smaller displacement along d, see straight_to_sphere).
        ray_rk4_step(medium, state, -sigma);
        straight_to_sphere(medium, state, centre, radius_m);
        return MarchEnd::entered;
      }
    }
    ray_rk4_step(medium, state, -step_m);
  }
}

} // namespace obdeect::propagation
