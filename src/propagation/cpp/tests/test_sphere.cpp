// Sphere representation of frame G (docs/propagation/propagation.md, "Sphere representation"):
// profile slope and straight-line air delay (analytic), the straight sphere entry (many directions,
// horizontal and upward), inside emission, misses_sphere, below_ground, the 3D ray integrator
// (constant index: exact; downward against bend_quadrature; upward against bend_upward; a turning
// ray against an independent Snell integration; sagitta of a horizontal ray), the hybrid rule
// against the integrator in the overlap region, the spherical geometry of a horizontal photon,
// ground_to_sky of sphere rows, and independence of chunking and row order.
#include "check.hpp"
#include "stream_fixtures.hpp"

#include "obdeect/propagation/ground_to_sky.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/ray_march.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/refraction_horizon.hpp"
#include "obdeect/propagation/sky.hpp"
#include "obdeect/propagation/sphere.hpp"
#include "obdeect/propagation/timing.hpp"

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using namespace obdeect_test;

using obdeect::propagation::AtmosphereProfile;
using obdeect::propagation::bend_quadrature;
using obdeect::propagation::bend_upward;
using obdeect::propagation::DropCounts;
using obdeect::propagation::RayMedium;
using obdeect::propagation::RayState;
using obdeect::propagation::RefractionTables;
using obdeect::propagation::SkyPlane;
using obdeect::stream::Column;
using obdeect::stream::PhotonChunk;

constexpr double kObs = 2000.0;                                 // observation level [m]
constexpr double kC = obdeect::stream::speed_of_light_m_per_ns; // [m/ns]
constexpr double kR = obdeect::stream::earth_radius_m;          // [m]

AtmosphereProfile ecmwf_profile() {
  return AtmosphereProfile::from_file(data_file("atmprof_ecmwf_north_winter_fixed.dat"));
}

// The ECMWF profile with n - 1 = nu at every level (density and thickness unchanged): no bending.
AtmosphereProfile constant_index_profile(double nu) {
  std::ifstream file(data_file("atmprof_ecmwf_north_winter_fixed.dat"));
  std::string line;
  std::string text;
  while (std::getline(file, line)) {
    std::istringstream row(line);
    std::vector<std::string> tokens;
    std::string token;
    while (row >> token) {
      tokens.push_back(token);
    }
    if (tokens.size() < 4 || tokens[0][0] == '#') {
      continue;
    }
    char buf[256];
    std::snprintf(buf, sizeof buf, "%s %s %s %.17g\n", tokens[0].c_str(), tokens[1].c_str(),
                  tokens[2].c_str(), nu);
    text += buf;
  }
  return AtmosphereProfile(text);
}

// Flat plane: telescope at (tx, ty, tz) above the observation level; no profile, no tables.
SkyPlane plane_at(double tx, double ty, double tz) {
  SkyPlane plane;
  plane.telescope_x_m = tx;
  plane.telescope_y_m = ty;
  plane.telescope_z_m = tz;
  plane.observation_m = kObs;
  return plane;
}

// One sky chunk row: x y z cx cy cz t (w = 1).
PhotonChunk sky_chunk(const std::vector<std::array<double, 7>> &rows) {
  PhotonChunk chunk;
  chunk.frame = obdeect::stream::Frame::sky;
  chunk.resize(
      rows.size(),
      obdeect::stream::column_bit(Column::x) | obdeect::stream::column_bit(Column::y) |
          obdeect::stream::column_bit(Column::z) | obdeect::stream::column_bit(Column::cx) |
          obdeect::stream::column_bit(Column::cy) | obdeect::stream::column_bit(Column::cz) |
          obdeect::stream::column_bit(Column::t) | obdeect::stream::column_bit(Column::w));
  for (std::size_t r = 0; r < rows.size(); ++r) {
    const auto &row = rows[r];
    chunk[Column::x][r] = static_cast<float>(row[0]);
    chunk[Column::y][r] = static_cast<float>(row[1]);
    chunk[Column::z][r] = static_cast<float>(row[2]);
    chunk[Column::cx][r] = static_cast<float>(row[3]);
    chunk[Column::cy][r] = static_cast<float>(row[4]);
    chunk[Column::cz][r] = static_cast<float>(row[5]);
    chunk[Column::t][r] = static_cast<float>(row[6]);
    chunk[Column::w][r] = 1.0F;
  }
  return chunk;
}

// One float32 unit in the last place of v (the storage resolution of a time or position).
double float_ulp(double v) {
  const float f = static_cast<float>(v);
  return static_cast<double>(std::nextafter(f, std::numeric_limits<float>::infinity()) - f);
}

// Signed cz from cx, cy in float32 storage (as the chunks hold it).
double signed_cz(double cx, double cy, double sign) {
  const double fx = static_cast<float>(cx), fy = static_cast<float>(cy);
  return sign * static_cast<float>(std::sqrt(std::max(0.0, 1.0 - fx * fx - fy * fy)));
}

// Dropped-row test: w = 0 and zero columns.
bool is_dropped(const PhotonChunk &chunk, std::size_t row) {
  return chunk[Column::w][row] == 0.0F && chunk[Column::x][row] == 0.0F;
}

// Integrates the ray with fixed RK4 steps until it crosses the level z_stop (flat, relative to the
// observation level) moving downward (downward = true: d_z < 0 and z <= z_stop) or upward. The
// crossing of the last step is interpolated linearly. False if the ray never crosses z_stop.
bool march_to_level(const RayMedium &medium, RayState &state, double z_stop, double step_m,
                    bool downward) {
  for (std::size_t n = 0; n < 100000000; ++n) {
    RayState previous = state;
    obdeect::propagation::ray_rk4_step(medium, state, step_m);
    const bool crossed = downward ? (state.d[2] < 0.0 && state.x[2] <= z_stop)
                                  : (state.d[2] > 0.0 && state.x[2] >= z_stop);
    if (crossed) {
      const double f = (z_stop - previous.x[2]) / (state.x[2] - previous.x[2]);
      for (std::size_t i = 0; i < 3; ++i) {
        state.x[i] = previous.x[i] + f * (state.x[i] - previous.x[i]);
        state.d[i] = previous.d[i] + f * (state.d[i] - previous.d[i]);
      }
      state.t = previous.t + f * (state.t - previous.t);
      return true;
    }
  }
  return false;
}

// Independent reference for a turning ray (flat geometry, k = 1): emitted upward at h_e with
// direction cosine s_e = sin(theta_e), turning at the apex and returning to h_e. Snell's law gives
// s(h) = s_e (1 + nu_e) / (1 + nu(h)), and the ray's horizontal displacement and time over the
// return to h_e are 2 * integral over h in [h_e, h_a] of tan(theta) dh and of (1 + nu) / cos(theta)
// dh / c. Parametrised by v = sqrt(1 - s), so that 1 - s^2 = v^2 (2 - v^2) exactly: the apex
// singularity is removed without any cancellation. The height of each node is found by bisection
// of nu(h) = s_e (1 + nu_e) / s - 1. Midpoint rule with n nodes; independent of the ODE.
struct SnellLeg {
  double dx{};
  double dt{};
};

SnellLeg snell_turning_reference(const AtmosphereProfile &profile, double h_e, double s_e,
                                 std::size_t n) {
  const double nu_e = profile.refractive_index_minus_one(h_e);
  const double target_scale = s_e * (1.0 + nu_e);
  const double v_e = std::sqrt(1.0 - s_e);
  SnellLeg leg;
  for (std::size_t i = 0; i < n; ++i) {
    const double v = (static_cast<double>(i) + 0.5) * v_e / static_cast<double>(n);
    const double s = 1.0 - v * v;
    const double target = target_scale / s - 1.0; // nu at the node
    double lo = h_e, hi = h_e + 400.0;            // nu decreases with height
    for (int it = 0; it < 80; ++it) {
      const double mid = 0.5 * (lo + hi);
      (profile.refractive_index_minus_one(mid) > target ? lo : hi) = mid;
    }
    const double h = 0.5 * (lo + hi);
    const double slope = profile.index_minus_one_slope_cm(100.0 * h) * 100.0; // d nu / dh [1/m]
    const double dh_ds = (-target_scale / (s * s)) / slope;                   // dh / ds > 0
    const double weight = 2.0 * (v_e / static_cast<double>(n)) / std::sqrt(2.0 - v * v);
    leg.dx += 2.0 * weight * s * dh_ds;
    leg.dt += 2.0 * weight * (1.0 + profile.refractive_index_minus_one(h)) * dh_ds / kC;
  }
  return leg;
}

void test_profile_slope() {
  const AtmosphereProfile profile = ecmwf_profile();
  // Heights off the tabulated levels (n - 1 has a kink at a level; the central difference would
  // average the one-sided slopes there).
  for (const double h_m : {2345.7, 5123.4, 9123.4, 15555.5, 25777.7, 45432.1, 80123.4}) {
    const double h_cm = 100.0 * h_m;
    const double step_cm = 0.001;
    const double numeric =
        (profile.index_minus_one_cm(h_cm + step_cm) - profile.index_minus_one_cm(h_cm - step_cm)) /
        (2.0 * step_cm);
    const double analytic = profile.index_minus_one_slope_cm(h_cm);
    // Central difference of an exponential cell: relative error ~ (1 cm rate)^2.
    CHECK_NEAR(analytic, numeric, 1e-6 * std::abs(numeric) + 1e-20);
  }
  // Constant n - 1: zero slope.
  CHECK_NEAR(constant_index_profile(3e-4).index_minus_one_slope_cm(5e5), 0.0, 1e-18);
}

void test_straight_line_delay() {
  const AtmosphereProfile profile = ecmwf_profile();
  // Non-horizontal line: heights linear in the path, equal to the zenith-angle integral.
  for (const double zenith : {0.2, 0.9, 1.3, 1.5}) {
    const double low = kObs + 300.0;
    const double high = 9000.0;
    const double length = (high - low) / std::cos(zenith);
    const double reference =
        obdeect::propagation::straight_path_air_delay_ns(profile, low, high, zenith);
    const double value =
        obdeect::propagation::straight_line_air_delay_ns(profile, low, high, length, 1.0);
    CHECK_NEAR(value, reference, 1e-9 * std::abs(reference));
  }
  // Horizontal line: exact, nu(h) L / c.
  const double h = 12345.0;
  const double length = 4000.0;
  CHECK_NEAR(obdeect::propagation::straight_line_air_delay_ns(profile, h, h, length, 1.0),
             profile.refractive_index_minus_one(h) * length / kC, 1e-15);
  // k_group scales linearly.
  CHECK_NEAR(obdeect::propagation::straight_line_air_delay_ns(profile, h, h, length, 2.5),
             2.5 * profile.refractive_index_minus_one(h) * length / kC, 1e-15);
}

// Straight entries: many directions (horizontal excluded), vacuum. The telescope is at z = 40 m so
// that upward photons start below the sphere (z = 5 m) outside it.
void test_straight_entry() {
  const double R = 12.0;
  const std::array<double, 3> centre{7.0, -4.0, 40.0};
  const SkyPlane plane = plane_at(centre[0], centre[1], centre[2]);
  int checked = 0;
  for (double az = 0.0; az < 6.3; az += 0.7) {
    for (const double cz_mag : {0.05, 0.3, 0.7, 0.99}) {
      for (const double sign : {-1.0, 1.0}) {
        const double cx = static_cast<float>(std::sqrt(1.0 - cz_mag * cz_mag) * std::cos(az));
        const double cy = static_cast<float>(std::sqrt(1.0 - cz_mag * cz_mag) * std::sin(az));
        const double cz = signed_cz(cx, cy, sign);
        // Line through the centre offset by 6 m (inside R): Q = C + 6 e, e perpendicular to (cx,
        // cy).
        const double e_norm = std::sqrt(cx * cx + cy * cy);
        const double qx = centre[0] - 6.0 * cy / e_norm, qy = centre[1] + 6.0 * cx / e_norm;
        const double qz = centre[2];
        // Start: downward 300 m before Q (above), upward so that the start is at z = 5 m.
        const double L = sign < 0.0 ? 300.0 : (centre[2] - 5.0) / cz;
        const double px = qx - L * cx, py = qy - L * cy, pz = qz - L * cz;
        PhotonChunk chunk = sky_chunk({{px, py, pz, cx, cy, cz, 0.0}});
        const DropCounts counts = obdeect::propagation::sky_to_sphere(chunk, plane, R);
        CHECK(counts.total() == 0);
        if (counts.total() != 0) {
          continue;
        }
        const double sx = chunk[Column::x][0], sy = chunk[Column::y][0], sz = chunk[Column::z][0];
        CHECK_NEAR(std::sqrt(sx * sx + sy * sy + sz * sz), R, 1e-4);
        // The row lies on the straight line: |(row + C - P) x u| small, and path = |row + C - P|.
        const double ax = sx + centre[0] - px, ay = sy + centre[1] - py, az2 = sz + centre[2] - pz;
        const double bx = ay * cz - az2 * cy, by = az2 * cx - ax * cz, bz = ax * cy - ay * cx;
        CHECK_NEAR(std::sqrt(bx * bx + by * by + bz * bz), 0.0, 2e-3);
        const double path = chunk[Column::path][0];
        CHECK_NEAR(path, std::sqrt(ax * ax + ay * ay + az2 * az2), 2e-3);
        // Stored as float32: one ulp of the time (about 6e-5 ns at 1000 ns) is the resolution.
        CHECK_NEAR(chunk[Column::t][0], path / kC, 2.0 * float_ulp(path / kC));
        CHECK_NEAR(chunk[Column::cz][0], cz, 1e-6);
        ++checked;
      }
    }
  }
  CHECK(checked == 72); // 9 azimuths x 4 cosines x 2 signs
}

void test_horizontal_and_upward_entry() {
  // Horizontal photon at the telescope altitude: straight in x, through the centre.
  const double R = 12.0;
  const SkyPlane plane = plane_at(0.0, 0.0, 3.0);
  PhotonChunk chunk = sky_chunk({{-200.0, 0.0, 3.0, 1.0, 0.0, 0.0, 0.0}});
  const DropCounts counts = obdeect::propagation::sky_to_sphere(chunk, plane, R);
  CHECK(counts.total() == 0);
  CHECK_NEAR(chunk[Column::x][0], -R, 1e-4);
  CHECK_NEAR(chunk[Column::y][0], 0.0, 1e-4);
  CHECK_NEAR(chunk[Column::z][0], 0.0, 1e-4);
  CHECK_NEAR(chunk[Column::path][0], 200.0 - R, 1e-4);
  CHECK_NEAR(chunk[Column::cz][0], 0.0, 1e-7);
  // Upward photon from below the sphere (telescope at z = 40 m, sphere bottom at 28 m): enters at
  // the bottom, z = 28 m, which is z = -R relative to the centre.
  const SkyPlane high = plane_at(0.0, 0.0, 40.0);
  PhotonChunk up = sky_chunk({{0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 0.0}});
  CHECK(obdeect::propagation::sky_to_sphere(up, high, R).total() == 0);
  CHECK_NEAR(up[Column::z][0], -R, 1e-4);
  CHECK_NEAR(up[Column::path][0], 23.0, 1e-4);
  CHECK_NEAR(up[Column::cz][0], 1.0, 1e-7);
}

void test_inside_below_and_missed() {
  const double R = 12.0;
  const SkyPlane plane = plane_at(0.0, 0.0, 3.0);
  // Inside: emission at 5 m from the centre.
  PhotonChunk inside = sky_chunk({{3.0, 4.0, 3.0, 0.0, 0.0, -1.0, 2.5}});
  CHECK(obdeect::propagation::sky_to_sphere(inside, plane, R).total() == 0);
  CHECK_NEAR(inside[Column::x][0], 3.0, 1e-5);
  CHECK_NEAR(inside[Column::y][0], 4.0, 1e-5);
  CHECK_NEAR(inside[Column::z][0], 0.0, 1e-5);
  CHECK_NEAR(inside[Column::path][0], 0.0, 0.0);
  CHECK_NEAR(inside[Column::t][0], 2.5, 1e-6);
  CHECK_NEAR(inside[Column::zem][0], kObs + 3.0, 1e-3);
  // Misses: a downward line at 20 m lateral offset (impact > R).
  PhotonChunk miss = sky_chunk({{20.0, 0.0, 500.0, 0.0, 0.0, -1.0, 0.0}});
  const DropCounts lost = obdeect::propagation::sky_to_sphere(miss, plane, R);
  CHECK(lost.misses_sphere == 1);
  CHECK(lost.total() == 1);
  CHECK(is_dropped(miss, 0));
  // Moving away: a line pointing away from the sphere also misses.
  PhotonChunk away = sky_chunk({{0.0, 0.0, 500.0, 0.0, 0.0, 1.0, 0.0}});
  CHECK(obdeect::propagation::sky_to_sphere(away, plane, R).misses_sphere == 1);
  // Below ground: emitted at z < 0.
  PhotonChunk below = sky_chunk({{0.0, 0.0, -5.0, 0.0, 0.0, -1.0, 0.0}});
  CHECK(obdeect::propagation::sky_to_sphere(below, plane, R).below_ground == 1);
  // Below ground: a downward line from above that reaches z < 0 before the sphere.
  PhotonChunk through = sky_chunk({{-800.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0}});
  // Horizontal at z = 1 m: enters the sphere (z range 3 +- 12 contains 1): not below ground.
  CHECK(obdeect::propagation::sky_to_sphere(through, plane, R).total() == 0);
  PhotonChunk deep = sky_chunk({{-800.0, 0.0, -30.0, 0.0, 0.0, -1.0, 0.0}});
  // Emitted below ground: counted below_ground before the geometry.
  CHECK(obdeect::propagation::sky_to_sphere(deep, plane, R).below_ground == 1);
  // Steep downward line that enters the sphere at z = -2 m, below the observation level.
  const double s = 0.8;
  const double cz = -std::sqrt(1.0 - s * s);
  const double x0 = 0.0, z0 = 3.0 + 40.0;
  // Line through (0, 0, -2) region: start at x = -40 * s / |cz|... choose start so the line hits
  // the centre's bottom region (z = 3 - 12 = -9 m is the bottom of the sphere; the centre line at
  // z = 3 - ...). Compute start so that the line passes (0, 0, 3 - 11.9).
  const double target_z = 3.0 - 11.9;
  const double L = (z0 - target_z) / -cz;
  PhotonChunk low = sky_chunk({{x0 - L * s, 0.0, z0, s, 0.0, cz, 0.0}});
  // The entry is at z below 0 (the sphere bottom is -9 m): the photon is below ground.
  const DropCounts counts = obdeect::propagation::sky_to_sphere(low, plane, R);
  CHECK(counts.below_ground == 1);
}

// Integrator: constant index (no bending): exact against the straight line with delay.
void test_integrator_constant_index() {
  const double nu = 3.0e-4;
  const AtmosphereProfile profile = constant_index_profile(nu);
  const double R = 12.0;
  const std::array<double, 3> centre{4.0, -2.0, 25.0};
  const RayMedium medium{&profile, 1.0, 1.0, kObs, 0.0};
  const std::array<double, 3> u{0.6, 0.48, -0.64};
  RayState state;
  state.x = {centre[0] - 1500.0 * u[0] + 3.0, centre[1] - 1500.0 * u[1], centre[2] - 1500.0 * u[2]};
  state.d = u;
  state.t = 0.0;
  const RayState start = state;
  const auto end = obdeect::propagation::march_to_sphere(medium, state, centre, R, 50.0);
  CHECK(end == obdeect::propagation::MarchEnd::entered);
  // Straight entry (closed form along the line).
  const std::array<double, 3> q{start.x[0] - centre[0], start.x[1] - centre[1],
                                start.x[2] - centre[2]};
  const double b = q[0] * u[0] + q[1] * u[1] + q[2] * u[2];
  const double c = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] - R * R;
  const double s_exact = -b - std::sqrt(b * b - c);
  CHECK_NEAR(state.path, s_exact, 1e-6);
  for (std::size_t i = 0; i < 3; ++i) {
    CHECK_NEAR(state.x[i], start.x[i] + s_exact * u[i], 1e-6);
  }
  // Time: path times n_group, straight.
  CHECK_NEAR(state.t, s_exact * (1.0 + nu) / kC, 1e-9);
}

// Integrator against bend_quadrature for downward rays (flat, ECMWF tables). bend_quadrature takes
// the straight-line position at the observation level and returns the bent one.
void test_integrator_matches_downward_bending() {
  const AtmosphereProfile profile = ecmwf_profile();
  const RefractionTables tables(profile, kObs);
  const RayMedium medium{&profile, 1.0, 1.0, kObs, 0.0};
  for (const double s_e : {0.3, 0.8}) {
    const double h_e = 4000.0; // emission altitude a.s.l.
    const double cz = -std::sqrt(1.0 - s_e * s_e);
    const double z0 = h_e - kObs;
    const double x_straight = z0 * s_e / -cz; // straight line reaches z = 0 here
    double bcx = s_e, bcy = 0.0, bx = x_straight, by = 0.0, bt = 0.0;
    CHECK(bend_quadrature(tables, h_e, bcx, bcy, bx, by, bt));
    RayState state;
    state.x = {0.0, 0.0, z0};
    state.d = {s_e, 0.0, cz};
    state.t = 0.0;
    CHECK(march_to_level(medium, state, 0.0, 0.25, true));
    CHECK_NEAR(state.x[0], bx, 1e-5);
    CHECK_NEAR(state.x[1], by, 1e-5);
    CHECK_NEAR(state.d[0], bcx, 1e-8);
    CHECK_NEAR(state.t, bt, 1e-5);
    std::fprintf(stderr, "downward s=%.1f: integrator x %.7f bend %.7f (dx %.2e m), dt %.2e ns\n",
                 s_e, state.x[0], bx, state.x[0] - bx, state.t - bt);
  }
}

// Integrator against bend_upward for upward rays from below the plane (flat, ECMWF). bend_upward
// takes the straight-line position at the plane and returns the bent one.
void test_integrator_matches_upward_bending() {
  const AtmosphereProfile profile = ecmwf_profile();
  const RayMedium medium{&profile, 1.0, 1.0, kObs, 0.0};
  const double h_e = kObs + 200.0;
  const double plane_m = kObs + 1500.0;
  const double s_e = 0.7;
  const double cz = std::sqrt(1.0 - s_e * s_e);
  const double x_straight = s_e * (plane_m - h_e) / cz;
  double cx = s_e, cy = 0.0, x = x_straight, y = 0.0, t = 0.0;
  CHECK(bend_upward(profile, h_e, plane_m, cx, cy, x, y, t, 1.0, 1.0));
  RayState state;
  state.x = {0.0, 0.0, h_e - kObs};
  state.d = {s_e, 0.0, cz};
  state.t = 0.0;
  CHECK(march_to_level(medium, state, plane_m - kObs, 0.25, false));
  CHECK_NEAR(state.x[0], x, 1e-5);
  CHECK_NEAR(state.d[0], cx, 1e-8);
  CHECK_NEAR(state.t, t, 1e-5);
  std::fprintf(stderr, "upward: integrator x %.7f bend_upward %.7f, dt %.2e ns\n", state.x[0], x,
               state.t - t);
}

// A turning ray: emitted upward, apex 300 m above, returns to its emission level. Independent
// Snell integration (snell_turning_reference) against the integrator.
void test_turning_ray_against_snell() {
  const AtmosphereProfile profile = ecmwf_profile();
  const RayMedium medium{&profile, 1.0, 1.0, kObs, 0.0};
  const double h_e = kObs + 10.0;
  const double nu_e = profile.refractive_index_minus_one(h_e);
  // Apex 300 m above: (1 + nu_e) s_e = 1 + nu(h_a) with h_a = h_e + 300 m (nu(h_a) ~ nu_e - 1e-5).
  const double nu_a = profile.refractive_index_minus_one(h_e + 300.0);
  const double s_e = (1.0 + nu_a) / (1.0 + nu_e);
  CHECK(s_e < 1.0);
  const double cz = std::sqrt(1.0 - s_e * s_e);
  RayState state;
  state.x = {0.0, 0.0, h_e - kObs};
  state.d = {s_e, 0.0, cz};
  state.t = 0.0;
  // Back down to the emission level (z = h_e - obs), after the apex.
  CHECK(march_to_level(medium, state, h_e - kObs, 0.5, true));
  // The reference's quadrature noise is about 2e-3 m and 5e-3 ns over this 318 km path (the profile
  // has a kink at every level and cell, so the midpoint rule does not converge monotonically: 80
  // 000, 160 000, 320 000 and 640 000 nodes differ by up to 1e-2 ns). 640 000 nodes keep that noise
  // inside the tolerances below.
  const SnellLeg reference = snell_turning_reference(profile, h_e, s_e, 640000);
  std::fprintf(stderr, "turning reference: dx %.9f m, dt %.9f ns\n", reference.dx, reference.dt);
  CHECK(state.d[2] < 0.0);
  CHECK_NEAR(state.x[0], reference.dx, 5e-3);
  CHECK_NEAR(state.t, reference.dt, 1e-2);
  std::fprintf(stderr,
               "turning ray: dx integrator %.9f reference %.9f diff %.3e m; dt diff %.3e ns\n",
               state.x[0], reference.dx, state.x[0] - reference.dx, state.t - reference.dt);
}

// Horizontal photon at the telescope altitude: the sagitta of the ray over L is L^2 |nu'| / 2 (the
// ray curves toward the denser air below).
void test_horizontal_sagitta() {
  const AtmosphereProfile profile = ecmwf_profile();
  const RayMedium medium{&profile, 1.0, 1.0, kObs, 0.0};
  const double h0 = kObs + 15.0; // altitude of the photon at the origin
  const double L = 2000.0;
  RayState state;
  state.x = {0.0, 0.0, h0 - kObs};
  state.d = {1.0, 0.0, 0.0};
  state.t = 0.0;
  // Integrate in x: fixed steps of 0.5 m along the ray, until the horizontal distance is L.
  while (state.x[0] < L) {
    obdeect::propagation::ray_rk4_step(medium, state, 0.5);
  }
  const double sag = (h0 - kObs) - state.x[2];
  const double slope =
      std::abs(profile.index_minus_one_slope_cm(100.0 * h0)) * 100.0; // |nu'| [1/m]
  const double analytic = 0.5 * slope * L * L;
  std::fprintf(stderr, "sagitta: integrator %.6f m, analytic 1/2 |nu'| L^2 = %.6f m\n", sag,
               analytic);
  CHECK(sag > 0.0);
  CHECK(std::abs(sag - analytic) <= 1e-3 * analytic); // nu' varies over L: 3e-4 relative here
}

// Spherical geometry, vacuum: a horizontal photon along the flat frame's x axis at z = 0 passes a
// telescope at (L, 0, 0) at the spherical altitude obs + L^2 / (2 R) (the flat plane is curved).
void test_spherical_horizontal_vacuum() {
  const double L = 20000.0;
  const double h_t = obdeect::propagation::spherical_altitude_m(L * L, kObs, kR);
  // Exact: sqrt(L^2 + (R + obs)^2) - R; to second order L^2 / (2 (R + obs)).
  CHECK_NEAR(h_t - kObs, std::sqrt(L * L + (kR + kObs) * (kR + kObs)) - kR - kObs, 1e-6);
  CHECK_NEAR(h_t - kObs, L * L / (2.0 * (kR + kObs)), 1e-3);
  SkyPlane plane = plane_at(L, 0.0, 0.0);
  plane.earth_radius_m = kR;
  plane.telescope_altitude_m = h_t;
  PhotonChunk chunk = sky_chunk({{0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}});
  const DropCounts counts = obdeect::propagation::sky_to_sphere(chunk, plane, 12.0);
  CHECK(counts.total() == 0);
  CHECK_NEAR(chunk[Column::x][0], -12.0, 1e-4);
  CHECK_NEAR(chunk[Column::path][0], L - 12.0, 1e-3);
  CHECK_NEAR(chunk[Column::zem][0], kObs, 1e-3);
}

// Hybrid rule against the integrator in the overlap region |cz| = 0.1 to 0.3 (downward rays from
// well above the telescope, ECMWF tables). The row of sky_to_sphere (hybrid rule) and the forced
// integrator must agree within 1e-4 m and 1e-4 ns.
void test_hybrid_against_integrator() {
  const AtmosphereProfile profile = ecmwf_profile();
  const double R = 12.0;
  const double tz = 10.0;
  const RefractionTables tables(profile, kObs + tz);
  SkyPlane plane = plane_at(0.0, 0.0, tz);
  plane.profile = &profile;
  plane.tables = &tables;
  double worst_x = 0.0, worst_t = 0.0;
  int hybrid_rows = 0;
  for (const double cz_mag : {0.1, 0.15, 0.2, 0.25, 0.3}) {
    for (const double offset : {-6.0, -2.0, 0.0, 3.0, 8.0}) {
      const double cx = static_cast<float>(std::sqrt(1.0 - cz_mag * cz_mag));
      const double cz = signed_cz(cx, 0.0, -1.0);
      const double z_start = 25000.0; // sky z above the observation level
      // Start x such that the bent ray crosses the telescope plane at x = offset: the bending over
      // this path (kilometres at |cz| = 0.1) is removed by fixed-point iteration on plane_arrival.
      double x_start = offset - (z_start - tz) * cx / cz_mag;
      for (int iteration = 0; iteration < 4; ++iteration) {
        double xp = 0, yp = 0, cxp = 0, cyp = 0, tp = 0;
        obdeect::propagation::plane_arrival(plane, x_start, 0.0, z_start, cx, 0.0, cz_mag, 1.0,
                                            false, kObs + z_start, 0.0, 1.0, 1.0, xp, yp, cxp, cyp,
                                            tp);
        x_start += offset - xp;
      }
      PhotonChunk chunk = sky_chunk({{x_start, 0.0, z_start, cx, 0.0, cz, 0.0}});
      // Hybrid rule in double precision: plane-mode crossing and tangent line.
      const double sx = chunk[Column::x][0], sz = chunk[Column::z][0];
      const double scx = chunk[Column::cx][0];
      double x_p = 0, y_p = 0, cx_p = 0, cy_p = 0, t_p = 0;
      const double cz_plane = std::sqrt(1.0 - scx * scx);
      const auto reason =
          obdeect::propagation::plane_arrival(plane, sx, 0.0, sz, scx, 0.0, cz_plane, 1.0, false,
                                              kObs + sz, 0.0, 1.0, 1.0, x_p, y_p, cx_p, cy_p, t_p);
      if (reason != obdeect::propagation::Drop::none) {
        continue;
      }
      std::array<double, 3> entry_rel{}, u_t{};
      double t_h = 0.0, path_h = 0.0;
      const auto hr =
          obdeect::propagation::hybrid_sphere_row(plane, R, x_p, y_p, cx_p, cy_p, false, t_p, 1.0,
                                                  1.0, {sx, 0.0, sz}, entry_rel, u_t, t_h, path_h);
      if (hr != obdeect::propagation::Drop::none) {
        continue;
      }
      // Forced integrator from the same start point, in double precision.
      RayState start;
      start.x = {sx, 0.0, sz};
      start.d = {scx, 0.0, -cz_plane};
      start.t = 0.0;
      RayState entry;
      const auto ir = obdeect::propagation::integrate_sphere_row(plane, R, start, 1.0, 1.0, entry);
      if (ir != obdeect::propagation::Drop::none) {
        continue;
      }
      ++hybrid_rows;
      const double dx = entry_rel[0] - (entry.x[0] - 0.0);
      const double dy = entry_rel[1] - entry.x[1];
      const double dz = entry_rel[2] - (entry.x[2] - tz);
      worst_x = std::max(worst_x, std::sqrt(dx * dx + dy * dy + dz * dz));
      worst_t = std::max(worst_t, std::abs(t_h - entry.t));
    }
  }
  std::fprintf(
      stderr,
      "hybrid vs integrator (double): %d rows, max position diff %.3e m, max time diff %.3e ns\n",
      hybrid_rows, worst_x, worst_t);
  CHECK(hybrid_rows >= 10);
  CHECK(worst_x <= 1e-4);
  CHECK(worst_t <= 1e-4);
}

// Hybrid rule against the integrator (forced from the emission) for shallow rays: the plane-mode
// crossing Q with |cz| = 0.1 to 0.3 and impact parameters |Q| up to R (the entry before Q for |Q| <
// R, after it otherwise). The hybrid integrates from Q, so both agree to the integrator's step.
void test_hybrid_shallow_against_integrator() {
  const AtmosphereProfile profile = ecmwf_profile();
  const double R = 12.0;
  const double tz = 10.0;
  const RefractionTables tables(profile, kObs + tz);
  SkyPlane plane = plane_at(0.0, 0.0, tz);
  plane.profile = &profile;
  plane.tables = &tables;
  double worst_x = 0.0, worst_t = 0.0;
  int rows = 0;
  for (const double cz_mag : {0.1, 0.13, 0.16, 0.2, 0.25, 0.3}) {
    for (const double offset : {-12.0, -9.0, -6.0, -3.0, -1.0, 0.0, 1.0, 3.0, 6.0, 9.0, 12.0}) {
      const double cx = static_cast<float>(std::sqrt(1.0 - cz_mag * cz_mag));
      const double cz = signed_cz(cx, 0.0, -1.0);
      const double z_start = 25000.0;
      double x_start = offset - (z_start - tz) * cx / cz_mag;
      for (int iteration = 0; iteration < 4; ++iteration) {
        double xp = 0, yp = 0, cxp = 0, cyp = 0, tp = 0;
        obdeect::propagation::plane_arrival(plane, x_start, 0.0, z_start, cx, 0.0, cz_mag, 1.0,
                                            false, kObs + z_start, 0.0, 1.0, 1.0, xp, yp, cxp, cyp,
                                            tp);
        x_start += offset - xp;
      }
      PhotonChunk chunk = sky_chunk({{x_start, 0.0, z_start, cx, 0.0, cz, 0.0}});
      const double sx = chunk[Column::x][0], sz = chunk[Column::z][0];
      const double scx = chunk[Column::cx][0];
      double x_p = 0, y_p = 0, cx_p = 0, cy_p = 0, t_p = 0;
      const double cz_plane = std::sqrt(1.0 - scx * scx);
      if (obdeect::propagation::plane_arrival(plane, sx, 0.0, sz, scx, 0.0, cz_plane, 1.0, false,
                                              kObs + sz, 0.0, 1.0, 1.0, x_p, y_p, cx_p, cy_p,
                                              t_p) != obdeect::propagation::Drop::none) {
        continue;
      }
      std::array<double, 3> entry_rel{}, u_t{};
      double t_h = 0.0, path_h = 0.0;
      if (obdeect::propagation::hybrid_sphere_row(plane, R, x_p, y_p, cx_p, cy_p, false, t_p, 1.0,
                                                  1.0, {sx, 0.0, sz}, entry_rel, u_t, t_h,
                                                  path_h) != obdeect::propagation::Drop::none) {
        continue;
      }
      RayState start;
      start.x = {sx, 0.0, sz};
      start.d = {scx, 0.0, -cz_plane};
      RayState entry;
      if (obdeect::propagation::integrate_sphere_row(plane, R, start, 1.0, 1.0, entry) !=
          obdeect::propagation::Drop::none) {
        continue;
      }
      ++rows;
      const double dx = entry_rel[0] - entry.x[0];
      const double dy = entry_rel[1] - entry.x[1];
      const double dz = entry_rel[2] - (entry.x[2] - tz);
      const double rowerr = std::sqrt(dx * dx + dy * dy + dz * dz);
      worst_x = std::max(worst_x, rowerr);
      worst_t = std::max(worst_t, std::abs(t_h - entry.t));
    }
  }
  std::fprintf(stderr,
               "hybrid shallow vs integrator: %d rows, max position diff %.3e m, max time diff "
               "%.3e ns\n",
               rows, worst_x, worst_t);
  CHECK(rows >= 40);
  CHECK(worst_x <= 1e-4);
  CHECK(worst_t <= 1e-4);
}

// Photon 148 of the validation run (sphere representation, upward from 4.13 m above the observation
// level, |cz| = 0.16, telescope centre at 20 m, R = 12.5 m, ECMWF tables, no dispersion). Its plane
// crossing is at the telescope plane 97.8 m along the ray, 22.6 m after the sphere entry, so the
// hybrid rule integrates backwards from the crossing. The reference is an independent scipy
// integration of the ray equation from the emission (DOP853, rtol 1e-12; validation adapter
// sphere_geometry.py): entry (1.22508583, -2.34695243, -12.21642251) relative to the centre,
// t = 103.23170775791107 ns. The straight line from a 50 m step was off by 1.1e-4 ns (about 3e-5 m
// along the ray): the entry is taken after an RK4 step onto the straight estimate, then the
// residual.
void test_hybrid_photon_148() {
  const AtmosphereProfile profile = ecmwf_profile();
  const double R = 12.5;
  const double tz = 20.0;
  const RefractionTables tables(profile, kObs + tz);
  SkyPlane plane = plane_at(0.0, 0.0, tz);
  plane.profile = &profile;
  plane.tables = &tables;
  const std::array<double, 3> ref_entry{1.22508583, -2.34695243, -12.21642251};
  const double ref_t = 103.23170775791107;
  // Double precision: plane-mode crossing and the hybrid rule (no float32 storage).
  {
    const double x = -14.034931182861328, y = -18.51064109802246, z = 4.129739284515381;
    const double cx = 0.6773982048034668, cy = 0.7175125479698181;
    const double cz = std::sqrt(1.0 - cx * cx - cy * cy);
    const double t0 = 28.0711669921875;
    double x_p = 0, y_p = 0, cx_p = 0, cy_p = 0, t_p = 0;
    CHECK(obdeect::propagation::plane_arrival(plane, x, y, z, cx, cy, cz, 1.0, true, kObs + z, t0,
                                              1.0, 1.0, x_p, y_p, cx_p, cy_p,
                                              t_p) == obdeect::propagation::Drop::none);
    std::array<double, 3> entry_rel{}, u_out{};
    double t_out = 0.0, path_out = 0.0;
    CHECK(obdeect::propagation::hybrid_sphere_row(plane, R, x_p, y_p, cx_p, cy_p, true, t_p, 1.0,
                                                  1.0, {x, y, z}, entry_rel, u_out, t_out,
                                                  path_out) == obdeect::propagation::Drop::none);
    std::fprintf(stderr, "photon 148 (double): dt %.3e ns, dpos %.3e m\n", t_out - ref_t,
                 std::sqrt((entry_rel[0] - ref_entry[0]) * (entry_rel[0] - ref_entry[0]) +
                           (entry_rel[1] - ref_entry[1]) * (entry_rel[1] - ref_entry[1]) +
                           (entry_rel[2] - ref_entry[2]) * (entry_rel[2] - ref_entry[2])));
    CHECK_NEAR(t_out, ref_t, 1e-6);
    for (int i = 0; i < 3; ++i) {
      CHECK_NEAR(entry_rel[i], ref_entry[i], 1e-6);
    }
  }
  // The kernel as run (sky_to_sphere, float32 columns): the stored time has a float32 spacing of
  // 7.6e-6 ns at 103 ns.
  {
    PhotonChunk chunk =
        sky_chunk({{-14.034931182861328, -18.51064109802246, 4.129739284515381, 0.6773982048034668,
                    0.7175125479698181, 0.1621956080198288, 28.0711669921875}});
    const DropCounts counts = obdeect::propagation::sky_to_sphere(chunk, plane, R);
    CHECK(counts.total() == 0);
    CHECK_NEAR(static_cast<double>(chunk[Column::t][0]), ref_t, 1e-5);
    CHECK_NEAR(static_cast<double>(chunk[Column::x][0]), ref_entry[0], 1e-5);
    CHECK_NEAR(static_cast<double>(chunk[Column::y][0]), ref_entry[1], 1e-5);
    CHECK_NEAR(static_cast<double>(chunk[Column::z][0]), ref_entry[2], 1e-5);
  }
}

// Chunking and row order do not change the rows (the kernel is per row).
void test_chunk_and_order_invariance() {
  const AtmosphereProfile profile = ecmwf_profile();
  const double R = 12.0;
  const RefractionTables tables(profile, kObs + 10.0);
  SkyPlane plane = plane_at(3.0, 1.0, 10.0);
  plane.profile = &profile;
  plane.tables = &tables;
  std::vector<std::array<double, 7>> rows;
  for (int i = 0; i < 24; ++i) {
    const double cx = static_cast<float>(0.5 * std::cos(0.3 * i));
    const double cy = static_cast<float>(0.4 * std::sin(0.3 * i));
    const double sign = i % 3 == 0 ? 1.0 : -1.0; // a few upward rows (from below the plane)
    const double cz = signed_cz(cx, cy, sign);
    const double z = sign > 0 ? 3.0 : 4000.0 + 900.0 * i;
    rows.push_back({3.0 * i - 40.0, 2.0 * i, z, cx, cy, cz, 0.5 * i});
  }
  PhotonChunk whole = sky_chunk(rows);
  const DropCounts all = obdeect::propagation::sky_to_sphere(whole, plane, R);
  std::vector<std::array<double, 7>> first(rows.begin(), rows.begin() + 9);
  std::vector<std::array<double, 7>> second(rows.begin() + 9, rows.end());
  PhotonChunk a = sky_chunk(second);
  PhotonChunk b = sky_chunk(first);
  const DropCounts ca = obdeect::propagation::sky_to_sphere(a, plane, R);
  const DropCounts cb = obdeect::propagation::sky_to_sphere(b, plane, R);
  CHECK(all.total() == ca.total() + cb.total());
  CHECK(all.total() > 0);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const PhotonChunk &part = i < 9 ? b : a;
    const std::size_t local = i < 9 ? i : i - 9;
    for (const Column c : {Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::cz,
                           Column::t, Column::zem, Column::w, Column::path}) {
      CHECK(whole[c][i] == part[c][local]);
    }
  }
}

// The plane-mode arrival of the hybrid rule (plane_arrival) reproduces sky_to_ground for a row.
void test_plane_arrival_matches_plane_mode() {
  const AtmosphereProfile profile = ecmwf_profile();
  const RefractionTables tables(profile, kObs + 10.0);
  SkyPlane plane = plane_at(3.0, 1.0, 10.0);
  plane.profile = &profile;
  plane.tables = &tables;
  const double cx = static_cast<float>(0.2), cy = static_cast<float>(-0.1);
  const double cz = signed_cz(cx, cy, -1.0);
  const double x = -2000.0, y = 150.0, z = 14000.0;
  PhotonChunk chunk = sky_chunk({{x, y, z, cx, cy, cz, 0.0}});
  chunk.resize(1, chunk.mask() & ~obdeect::stream::column_bit(Column::cz));
  chunk[Column::x][0] = static_cast<float>(x);
  chunk[Column::y][0] = static_cast<float>(y);
  chunk[Column::z][0] = static_cast<float>(z);
  chunk[Column::cx][0] = static_cast<float>(cx);
  chunk[Column::cy][0] = static_cast<float>(cy);
  chunk[Column::t][0] = 0.0F;
  chunk[Column::w][0] = 1.0F;
  PhotonChunk plane_chunk = chunk;
  (void)obdeect::propagation::sky_to_ground(plane_chunk, plane);
  double x_out = 0, y_out = 0, cx_out = 0, cy_out = 0, t_out = 0;
  const double fz = static_cast<float>(z);
  const double cz_mag = std::sqrt(1.0 - cx * cx - cy * cy);
  const auto reason = obdeect::propagation::plane_arrival(
      plane, static_cast<float>(x), static_cast<float>(y), fz, cx, cy, cz_mag, 1.0, false,
      kObs + fz, 0.0, 1.0, 1.0, x_out, y_out, cx_out, cy_out, t_out);
  CHECK(reason == obdeect::propagation::Drop::none);
  CHECK_NEAR(x_out, plane_chunk[Column::x][0], 1e-3);
  CHECK_NEAR(y_out, plane_chunk[Column::y][0], 1e-3);
  CHECK_NEAR(t_out, plane_chunk[Column::t][0], 2.0 * float_ulp(t_out));
}

// Round trip: a straight sphere row (vacuum) back to the sky frame gives the emission point.
void test_ground_to_sky_round_trip() {
  const double R = 12.0;
  const SkyPlane plane = plane_at(7.0, -4.0, 3.0);
  // Start on a line 2 m from the centre (through the sphere), 300 m along the path from it.
  const double ux = static_cast<float>(0.3), uy = static_cast<float>(-0.2);
  const double uz = signed_cz(ux, uy, -1.0);
  const double px = 7.0 - 300.0 * ux + 1.0, py = -4.0 - 300.0 * uy + 2.0, pz = 3.0 - 300.0 * uz;
  PhotonChunk chunk = sky_chunk({{px, py, pz, ux, uy, uz, 4.0}});
  // Carry the stored (float32) start point into the comparison below.
  const DropCounts counts = obdeect::propagation::sky_to_sphere(chunk, plane, R);
  CHECK(counts.total() == 0);
  chunk.frame = obdeect::stream::Frame::ground;
  obdeect::propagation::GroundPlane ground{7.0, -4.0, 3.0, kObs, nullptr};
  std::uint64_t dropped = 0;
  const auto pieces = obdeect::propagation::ground_to_sky_sphere(chunk, ground, dropped);
  CHECK(dropped == 0);
  CHECK(pieces.size() == 1);
  CHECK_NEAR(pieces[0][Column::x][0], static_cast<float>(px), 2e-3);
  CHECK_NEAR(pieces[0][Column::y][0], static_cast<float>(py), 2e-3);
  CHECK_NEAR(pieces[0][Column::z][0], static_cast<float>(pz), 2e-3);
  // The path column is float32: its ulp (1.5e-5 m at 250 m) is 5e-5 ns of time, the time ulp adds.
  CHECK_NEAR(pieces[0][Column::t][0], 4.0,
             float_ulp(chunk[Column::path][0]) / kC + 2.0 * float_ulp(4.0));
}

void test_validation() {
  const SkyPlane plane = plane_at(0.0, 0.0, 3.0);
  PhotonChunk chunk = sky_chunk({{0.0, 0.0, 100.0, 0.0, 0.0, -1.0, 0.0}});
  CHECK(throws_with([&] { (void)obdeect::propagation::sky_to_sphere(chunk, plane, 0.0); },
                    "sphere radius must be finite and positive"));
  PhotonChunk ground = chunk;
  ground.frame = obdeect::stream::Frame::ground;
  CHECK(throws_with([&] { (void)obdeect::propagation::sky_to_sphere(ground, plane, 12.0); },
                    "needs photons in frame sky"));
}

} // namespace

// Spherical straight rule (sky_to_sphere, no refraction): a straight segment from an emission with
// z >= 0 to its entry point stays above the observation level. z is linear along the segment and
// the spherical altitude is at least obs + z (sqrt(rho^2 + (R + a)^2) >= R + a), so the segment
// cannot dip below the ground between its ends. Random grazing downward chords: every accepted row
// is checked by sampling the altitude along its segment from the emission to the entry point.
void test_straight_chords_stay_above_ground() {
  SkyPlane plane;
  plane.telescope_x_m = 0.0;
  plane.telescope_y_m = 0.0;
  plane.telescope_z_m = 5.0;
  plane.observation_m = kObs;
  plane.earth_radius_m = kR;
  plane.telescope_altitude_m = kObs + 5.0;
  const double radius = 12.0;
  std::mt19937_64 rng(20261123);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  const double pi = std::acos(-1.0);
  std::vector<std::array<double, 7>> rows;
  for (std::size_t i = 0; i < 6000; ++i) {
    const double phi = 2.0 * pi * uniform(rng);
    const double rho = 20000.0 * uniform(rng);
    const double cz = 0.02 + 0.18 * uniform(rng);
    const double c_xy = std::sqrt(1.0 - cz * cz);
    // Emission above the line to the telescope height at its horizontal distance: the line passes
    // near the sphere and stays at z >= 0 for the emission.
    const double z = 5.0 + rho * cz / c_xy + 8.0 * (uniform(rng) - 0.5);
    // Towards the telescope axis, with a random offset in azimuth.
    const double azimuth =
        std::atan2(-std::sin(phi), -std::cos(phi)) + 1.0e-3 * (uniform(rng) - 0.5);
    rows.push_back({static_cast<float>(rho * std::cos(phi)),
                    static_cast<float>(rho * std::sin(phi)), static_cast<float>(z),
                    static_cast<float>(c_xy * std::cos(azimuth)),
                    static_cast<float>(c_xy * std::sin(azimuth)), static_cast<float>(-cz), 0.0});
  }
  PhotonChunk chunk = sky_chunk(rows);
  (void)sky_to_sphere(chunk, plane, radius, false);
  std::size_t accepted = 0;
  for (std::size_t r = 0; r < rows.size(); ++r) {
    if (chunk[Column::w][r] == 0.0F) {
      continue;
    }
    ++accepted;
    const double ex = plane.telescope_x_m + chunk[Column::x][r];
    const double ey = plane.telescope_y_m + chunk[Column::y][r];
    const double ez = plane.telescope_z_m + chunk[Column::z][r];
    double min_altitude = 1e300;
    for (int k = 0; k <= 4000; ++k) {
      const double t = static_cast<double>(k) / 4000.0;
      const double x = rows[r][0] + t * (ex - rows[r][0]);
      const double y = rows[r][1] + t * (ey - rows[r][1]);
      const double zz = rows[r][2] + t * (ez - rows[r][2]);
      const double a = kObs + zz;
      min_altitude = std::min(min_altitude, std::sqrt(x * x + y * y + (kR + a) * (kR + a)) - kR);
    }
    CHECK(min_altitude >= kObs - 1e-2);
  }
  CHECK(accepted > 500);
}

int main() {
  return run_checks([] {
    test_profile_slope();
    test_straight_line_delay();
    test_straight_entry();
    test_horizontal_and_upward_entry();
    test_inside_below_and_missed();
    test_validation();
    test_straight_chords_stay_above_ground();
    test_integrator_constant_index();
    test_integrator_matches_downward_bending();
    test_integrator_matches_upward_bending();
    test_turning_ray_against_snell();
    test_horizontal_sagitta();
    test_spherical_horizontal_vacuum();
    test_hybrid_against_integrator();
    test_hybrid_shallow_against_integrator();
    test_hybrid_photon_148();
    test_chunk_and_order_invariance();
    test_plane_arrival_matches_plane_mode();
    test_ground_to_sky_round_trip();
  });
}
