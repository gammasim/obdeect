// Chord integrals and the extinction of sphere rows (docs/propagation/extinction.md, "Chord
// integral" and "Sphere rows"): closed forms of the chord quadrature (constant and exponential
// extinction, split at table heights), the table's exact extinction coefficient against finite
// differences, the flat chord against the plane formula, the spherical chord against an independent
// integral and its dip in closed form, continuity at the switch |u| = near_horizontal_cos (100 km
// chords), sphere rows against the plane rows (flat and spherical, same photon: Q at the telescope
// centre), zero-path rows, input checks, and independence of row order and ranges.
#include "check.hpp"

#include "obdeect/propagation/chord.hpp"
#include "obdeect/propagation/extinction.hpp"
#include "obdeect/propagation/transmission.hpp"

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using obdeect::propagation::TelescopeDepths;
using obdeect::propagation::TransmissionTable;
using obdeect::propagation::WavelengthInterpolation;
using obdeect::stream::Column;
using obdeect::stream::PhotonChunk;

constexpr double kObs = 2000.0;                        // observation level [m a.s.l.]
constexpr double kR = obdeect::stream::earth_radius_m; // Earth radius [m]
constexpr double kLambda = 355.0;                      // [nm]
constexpr double kTelescopeZ = 30.0;                   // telescope centre above H2 [m]

std::string data_file(const std::string &name) {
  return std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/" + name;
}

TransmissionTable site_table(WavelengthInterpolation mode) {
  return TransmissionTable::from_file(data_file("atm_trans_2000_1_10_0_0_2000.dat"), kObs, mode);
}

// Depth of a flat chord from the table's own optical depth: |tau(h_1) - tau(h_2)| / |u| (the plane
// formula with the chord's cosine u, |u| > 0).
double plane_depth(const TransmissionTable &table, double h_1, double h_2, double abs_u) {
  return std::abs(table.optical_depth<double>(kLambda, h_1) -
                  table.optical_depth<double>(kLambda, h_2)) /
         abs_u;
}

// Column mask for the sphere rows.
std::uint32_t sphere_mask() {
  std::uint32_t mask = 0;
  for (const Column c : {Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::cz,
                         Column::zem, Column::path, Column::wavelength}) {
    mask |= obdeect::stream::column_bit(c);
  }
  return mask;
}

// One sphere row (frame G) as stored: position relative to C, direction, emission altitude,
// chord length, wavelength.
struct Row {
  double x, y, z, cx, cy, cz, zem, path, lambda;
};

PhotonChunk sphere_chunk(const std::vector<Row> &rows) {
  PhotonChunk chunk;
  chunk.frame = obdeect::stream::Frame::ground;
  chunk.resize(rows.size(), sphere_mask());
  for (std::size_t r = 0; r < rows.size(); ++r) {
    const Row &row = rows[r];
    chunk[Column::x][r] = static_cast<float>(row.x);
    chunk[Column::y][r] = static_cast<float>(row.y);
    chunk[Column::z][r] = static_cast<float>(row.z);
    chunk[Column::cx][r] = static_cast<float>(row.cx);
    chunk[Column::cy][r] = static_cast<float>(row.cy);
    chunk[Column::cz][r] = static_cast<float>(row.cz);
    chunk[Column::zem][r] = static_cast<float>(row.zem);
    chunk[Column::path][r] = static_cast<float>(row.path);
    chunk[Column::wavelength][r] = static_cast<float>(row.lambda);
  }
  return chunk;
}

// Direction of a row with d.n = cos_n on the local vertical n = (0, 0, 1) (flat or at the
// telescope centre): cz = -cos_n (downward) or +cos_n (upward), cx = sqrt(1 - cos^2), cy = 0.
// cx is stored in float32 and cz is derived from that value, so that the plane kernels (airmass
// from cx, cy) and the sphere kernel (cz) see the same direction.
Row direction_row(double cos_n, bool upward, double x, double y, double z, double path,
                  double lambda) {
  const double cx = static_cast<float>(std::sqrt(1.0 - cos_n * cos_n));
  const double cz = std::sqrt(std::max(0.0, 1.0 - cx * cx));
  Row row{x, y, z, cx, 0.0, upward ? cz : -cz, 0.0, path, lambda};
  return row;
}

// Expected emission altitude of a row relative to its end Q (flat): Q altitude plus the vertical
// step of the chord, path * cos (downward: E above Q).
double flat_emission(double h_q, double path, double cos_n, bool upward) {
  return upward ? h_q - path * cos_n : h_q + path * cos_n;
}

template <class F> void expect_invalid(F &&f, const char *what) {
  try {
    f();
  } catch (const std::invalid_argument &) {
    return;
  } catch (...) {
    std::fprintf(stderr, "%s: unexpected exception type\n", what);
    ++check_failures;
    return;
  }
  std::fprintf(stderr, "%s: no exception\n", what);
  ++check_failures;
}

// Closed forms of the chord quadrature.
void test_chord_integral_closed_forms() {
  using obdeect::propagation::chord_integral;
  using obdeect::propagation::chord_integral_split;
  using obdeect::propagation::max_chord_breaks;
  // Constant extinction: k L.
  for (const double length : {1234.5, 1.0e5}) {
    CHECK_NEAR(chord_integral(length, 0.0, [](double) { return 2.5e-5; }), 2.5e-5 * length,
               1e-12 * 2.5e-5 * length);
  }
  // Exponential extinction k0 exp(-(h0 + g s) / H) along a chord with a constant slope g (m per m):
  // closed form k0 exp(-h0/H) (H/g) (1 - exp(-g L / H)); horizontal (g = 0): k0 exp(-h0/H) L.
  const double k0 = 3.0e-4, scale = 8000.0, h0 = 5000.0, length = 20000.0;
  const auto exponential = [&](double slope) {
    return [=](double s) { return k0 * std::exp(-(h0 + slope * s) / scale); };
  };
  const double g = 0.05;
  // The altitude changes by g length along the chord (monotone): the pieces follow the change.
  CHECK_NEAR(chord_integral(length, g * length, exponential(g)),
             k0 * std::exp(-h0 / scale) * (scale / g) * (1.0 - std::exp(-g * length / scale)),
             1e-12 * k0 * std::exp(-h0 / scale) * scale / g);
  CHECK_NEAR(chord_integral(length, 0.0, exponential(0.0)), k0 * std::exp(-h0 / scale) * length,
             1e-12 * k0 * length);
  // Splitting at breaks does not change the integral (the split parts are smooth).
  std::array<double, max_chord_breaks> breaks{};
  breaks[0] = 3000.0;
  breaks[1] = 11111.0;
  const auto altitude = [&](double s) { return h0 + g * s; };
  CHECK_NEAR(chord_integral_split(length, breaks, 2, altitude, exponential(g)),
             chord_integral(length, g * length, exponential(g)), 1e-13 * k0 * scale / g);
  CHECK_NEAR(chord_integral(0.0, 0.0, exponential(g)), 0.0, 0.0);
}

// The table's extinction coefficient is the derivative of its optical depth.
void test_extinction_coefficient(const TransmissionTable &table) {
  const double delta = 0.05; // m: central difference inside one table interval (no node within)
  for (const double h : {2100.3, 3333.3, 7300.0, 12345.6, 25000.7}) {
    const double fd = (table.optical_depth<double>(kLambda, h + delta) -
                       table.optical_depth<double>(kLambda, h - delta)) /
                      (2.0 * delta);
    const double k = table.extinction_coefficient<double>(kLambda, h);
    CHECK(k > 0.0);
    CHECK_NEAR(k, std::abs(fd), 1e-7 * k);
  }
  CHECK_NEAR(table.extinction_coefficient<double>(kLambda, kObs - 10.0), 0.0, 0.0);
  CHECK_NEAR(table.extinction_coefficient<double>(kLambda, kObs), 0.0, 0.0);
  CHECK_NEAR(table.extinction_coefficient<double>(kLambda, 1.0e7), 0.0, 0.0);
  // Linear blend in wavelength: the coefficient at 355.5 nm is the mean of 355 and 356 nm.
  if (table.wavelength_interpolation() == WavelengthInterpolation::linear) {
    const double h = 6100.0;
    const double mid = table.extinction_coefficient<double>(355.5, h);
    const double mean = 0.5 * (table.extinction_coefficient<double>(355.0, h) +
                               table.extinction_coefficient<double>(356.0, h));
    CHECK_NEAR(mid, mean, 1e-12 * mean);
  }
}

// A horizontal chord at constant height: k(h) L exactly.
void test_horizontal_flat_chord(const TransmissionTable &table) {
  const double h = 6000.0, length = 5000.0;
  CHECK_NEAR(obdeect::propagation::flat_chord_depth(table, kLambda, h, h, length),
             table.extinction_coefficient<double>(kLambda, h) * length,
             1e-12 * table.extinction_coefficient<double>(kLambda, h) * length);
}

// A flat chord with altitude linear in s: the integral of k over the altitude is the plane
// formula |tau(h_1) - tau(h_2)| / |u| with u = dh/ds, for any |u| (also the near-horizontal ones).
void test_flat_chord_matches_plane_formula(const TransmissionTable &table) {
  const double length = 3000.0, h_q = 4500.0;
  for (const double abs_u : {0.02, 0.05, 0.1, 0.3, 0.9}) {
    const double h_e = h_q + length * abs_u;
    const double chord = obdeect::propagation::flat_chord_depth(table, kLambda, h_e, h_q, length);
    CHECK_NEAR(chord, plane_depth(table, h_e, h_q, abs_u), 1e-9);
  }
}

// Spherical chord: a symmetric chord that dips. Start at altitude H with the direction at angle
// alpha below the horizontal, end at the same altitude after L = 2 (R + H) sin(alpha). Its midpoint
// (s = L/2) is the lowest point, at altitude (R + H) cos(alpha) - R. The integral is compared with
// an independent Simpson sum over the altitudes from the vector |v + s d| - R, split where the
// chord crosses a table height (the coefficient is discontinuous there; the crossings are found by
// bisection of the vector altitude, which is monotone on each half of the chord).
void test_spherical_chord_dip(const TransmissionTable &table) {
  // alpha small: the lowest point (R + H) cos(alpha) - R stays above the observation level. The
  // start height is not a table node (the reference evaluates the end points), and the chord
  // crosses the node at 4000 m twice (the split).
  const double h_start = 4200.0, alpha = 0.01;
  const double r0 = kR + h_start;
  const std::array<double, 3> v{0.0, 0.0, r0};
  const std::array<double, 3> d{std::cos(alpha), 0.0, -std::sin(alpha)};
  const double length = 2.0 * r0 * std::sin(alpha);
  // Altitude at the midpoint, from the vector sum (independent of the chord formulas).
  const double mid = length / 2.0;
  const double r_mid = std::sqrt((v[0] + mid * d[0]) * (v[0] + mid * d[0]) +
                                 (v[1] + mid * d[1]) * (v[1] + mid * d[1]) +
                                 (v[2] + mid * d[2]) * (v[2] + mid * d[2]));
  CHECK_NEAR(r_mid - kR, r0 * std::cos(alpha) - kR, 1e-6);
  CHECK(r_mid - kR < h_start); // the chord dips below its ends
  const auto altitude = [&](double s) {
    const double x = v[0] + s * d[0], y = v[1] + s * d[1], z = v[2] + s * d[2];
    return std::sqrt(x * x + y * y + z * z) - kR;
  };
  std::vector<double> cuts{0.0, mid, length};
  for (const double node : table.node_heights_m()) {
    for (const auto [lo0, hi0] :
         {std::array<double, 2>{0.0, mid}, std::array<double, 2>{mid, length}}) {
      double lo = lo0, hi = hi0;
      if ((altitude(lo) - node) * (altitude(hi) - node) >= 0.0) {
        continue;
      }
      for (int it = 0; it < 200; ++it) {
        const double m = 0.5 * (lo + hi);
        ((altitude(m) - node) * (altitude(lo) - node) > 0.0 ? lo : hi) = m;
      }
      cuts.push_back(0.5 * (lo + hi));
    }
  }
  std::sort(cuts.begin(), cuts.end());
  double simpson = 0.0;
  const int n_part = 200000; // the reference converges as 1/n^2 (the split points are exact)
  for (std::size_t c = 0; c + 1 < cuts.size(); ++c) {
    const double a = cuts[c], b = cuts[c + 1];
    if (!(b > a)) {
      continue;
    }
    const double step = (b - a) / n_part;
    double part = 0.0;
    for (int i = 0; i <= n_part; ++i) {
      const double s = a + i * step;
      const double weight = (i == 0 || i == n_part) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
      part += weight * table.extinction_coefficient<double>(kLambda, altitude(s));
    }
    simpson += part * step / 3.0;
  }
  const double chord =
      obdeect::propagation::spherical_chord_depth(table, kLambda, v, d, length, kR);
  CHECK(chord > 0.0);
  CHECK_NEAR(chord, simpson, 1e-8 * chord);
}

// The two forms of the extinction at |u| = near_horizontal_cos agree (chords up to 100 km): the
// plane form (flat) and the slant integral (spherical) against the chord integral, both table
// interpolations.
void test_continuity_at_switch(const TransmissionTable &table) {
  const double u = obdeect::propagation::near_horizontal_cos;
  for (const double length : {1.0e3, 1.0e4, 1.0e5}) {
    // Flat: Q at 4 km, E above Q by length * u.
    const double h_q = 4000.0;
    const double h_e = h_q + length * u;
    const double flat_plane = plane_depth(table, h_e, h_q, u);
    const double flat_chord =
        obdeect::propagation::flat_chord_depth(table, kLambda, h_e, h_q, length);
    CHECK_NEAR(flat_chord, flat_plane, 1e-9);
    // Spherical: Q at 4 km on the vertical, direction with d.n = -u (downward): slant integral
    // (transmission_curved's depth) against the spherical chord from E.
    const std::array<double, 3> v_q{0.0, 0.0, kR + h_q};
    const std::array<double, 3> d{std::sqrt(1.0 - u * u), 0.0, -u};
    const std::array<double, 3> v_e{v_q[0] - length * d[0], v_q[1] - length * d[1],
                                    v_q[2] - length * d[2]};
    const double h_e_sph = std::sqrt(v_e[0] * v_e[0] + v_e[1] * v_e[1] + v_e[2] * v_e[2]) - kR;
    const double slant = table.slant_optical_depth<double>(kLambda, h_q, h_e_sph, u, kR);
    const double sph_chord =
        obdeect::propagation::spherical_chord_depth(table, kLambda, v_e, d, length, kR);
    CHECK_NEAR(sph_chord, slant, 1e-9);
  }
}

// Sphere rows against the plane rows. The row point Q at the telescope centre (x = y = z = 0), so
// h_Q = h_T and the local vertical is the telescope's: the same photon in both representations.
// Flat: plane rows of atmospheric_acceptance (downward). Spherical: plane rows of
// atmospheric_acceptance_curved (downward and upward). Rows with |u| >= near_horizontal_cos only.
void test_sphere_rows_match_plane_rows(const TransmissionTable &table) {
  using obdeect::propagation::CurvedPath;
  const double telescope_z = 200.0;
  const double h_t = kObs + telescope_z;
  const TelescopeDepths depths = table.telescope_depths(h_t, 20.0);
  std::vector<Row> rows;
  for (const double cos_n : {0.1, 0.25, 0.6, 0.95}) {
    for (const double path : {2000.0, 30000.0, 90000.0}) {
      rows.push_back(direction_row(cos_n, false, 0.0, 0.0, 0.0, path, kLambda));
      rows.back().zem = flat_emission(h_t, path, cos_n, false);
    }
  }
  // Emitted inside the sphere: path 0, no atmospheric depth.
  rows.push_back(direction_row(0.5, false, 0.0, 0.0, 0.0, 0.0, kLambda));
  rows.back().zem = h_t;
  // Flat comparison.
  {
    PhotonChunk sphere = sphere_chunk(rows);
    PhotonChunk plane;
    plane.frame = obdeect::stream::Frame::ground;
    std::uint32_t mask = 0;
    for (const Column c : {Column::cx, Column::cy, Column::zem, Column::wavelength}) {
      mask |= obdeect::stream::column_bit(c);
    }
    plane.resize(rows.size(), mask);
    for (std::size_t r = 0; r < rows.size(); ++r) {
      plane[Column::cx][r] = static_cast<float>(rows[r].cx);
      plane[Column::cy][r] = static_cast<float>(rows[r].cy);
      plane[Column::zem][r] = static_cast<float>(rows[r].zem);
      plane[Column::wavelength][r] = static_cast<float>(rows[r].lambda);
    }
    std::vector<double> p_plane(rows.size(), 0.0), p_sphere(rows.size(), 0.0);
    obdeect::propagation::atmospheric_acceptance<double>(plane, table, depths, p_plane);
    obdeect::propagation::SphereGeometry geometry{kObs, 0.0, 0.0, 0.0, telescope_z};
    obdeect::propagation::atmospheric_acceptance_sphere<double>(sphere, table, depths, geometry,
                                                                p_sphere);
    for (std::size_t r = 0; r < rows.size(); ++r) {
      if (rows[r].path == 0.0) {
        CHECK_NEAR(p_sphere[r], std::exp(-depths.od_focus[static_cast<std::size_t>(kLambda)]),
                   1e-12);
        continue;
      }
      CHECK(p_plane[r] > 0.0);
      CHECK_NEAR(p_sphere[r], p_plane[r], 1e-6 * p_plane[r]);
    }
  }
  // Spherical comparison: downward and upward rows, Q above the telescope centre.
  {
    std::vector<Row> curved_rows;
    for (const double cos_n : {0.1, 0.25, 0.6, 0.95}) {
      for (const double path : {2000.0, 30000.0, 90000.0}) {
        curved_rows.push_back(direction_row(cos_n, false, 0.0, 0.0, 0.0, path, kLambda));
        curved_rows.back().zem = h_t + path * cos_n;
      }
      // Upward rows start below the telescope (the vertical step stays above the observation
      // level).
      for (const double path : {100.0, 150.0}) {
        curved_rows.push_back(direction_row(cos_n, true, 0.0, 0.0, 0.0, path, kLambda));
        curved_rows.back().zem = h_t - path * cos_n;
      }
    }
    PhotonChunk sphere = sphere_chunk(curved_rows);
    PhotonChunk plane;
    plane.frame = obdeect::stream::Frame::ground;
    std::uint32_t mask = 0;
    for (const Column c : {Column::cx, Column::cy, Column::cz, Column::zem, Column::wavelength}) {
      mask |= obdeect::stream::column_bit(c);
    }
    plane.resize(curved_rows.size(), mask);
    for (std::size_t r = 0; r < curved_rows.size(); ++r) {
      plane[Column::cx][r] = static_cast<float>(curved_rows[r].cx);
      plane[Column::cy][r] = static_cast<float>(curved_rows[r].cy);
      plane[Column::cz][r] = static_cast<float>(curved_rows[r].cz);
      plane[Column::zem][r] = static_cast<float>(curved_rows[r].zem);
      plane[Column::wavelength][r] = static_cast<float>(curved_rows[r].lambda);
    }
    const CurvedPath path{kR, h_t, 0.0, 0.0, 1.0};
    std::vector<double> p_plane(curved_rows.size(), 0.0), p_sphere(curved_rows.size(), 0.0);
    obdeect::propagation::atmospheric_acceptance_curved<double>(plane, table, depths, path,
                                                                p_plane);
    obdeect::propagation::SphereGeometry geometry{kObs, kR, 0.0, 0.0, telescope_z};
    obdeect::propagation::atmospheric_acceptance_sphere<double>(sphere, table, depths, geometry,
                                                                p_sphere);
    for (std::size_t r = 0; r < curved_rows.size(); ++r) {
      CHECK(p_plane[r] > 0.0);
      CHECK_NEAR(p_sphere[r], p_plane[r], 1e-6 * p_plane[r]);
    }
  }
}

// Sphere row 77 of the validation run: near-horizontal (u = -2e-6), path 3714 m at 2031 m, 401.27
// nm, direction with cz = 0 and |d|^2 - 1 = -5.7e-8 (float32 cx, cy). Such a norm error moves the
// chord formulas (which assume a unit d) by about 0.3 m in altitude at the Earth radius, a depth
// error of 3e-5 relative. The reference integrates the table's extinction coefficient by Simpson's
// rule along the chord from E = Q - L d (d normalised in double), with the altitudes from the
// vector norms (no chord.hpp code).
void test_sphere_row_chord_normalised(const TransmissionTable &table) {
  const Row row{-8.1801405, 9.4517355, 0.0,        0.654408336, -0.756141305,
                0.0,        2031.0896, 3714.27661, 401.269989};
  const PhotonChunk chunk = sphere_chunk({row});
  const double x = chunk[Column::x][0], y = chunk[Column::y][0], z = chunk[Column::z][0];
  const double cx = chunk[Column::cx][0], cy = chunk[Column::cy][0], cz = chunk[Column::cz][0];
  const double length = chunk[Column::path][0], lambda = chunk[Column::wavelength][0];
  // Reference: Earth-centred coordinates (the Earth centre at (0, 0, -(R + H2)) in the sky frame).
  const std::array<double, 3> q{x, y, kR + kObs + kTelescopeZ + z};
  const double norm = std::sqrt(cx * cx + cy * cy + cz * cz);
  const std::array<double, 3> d{cx / norm, cy / norm, cz / norm};
  const std::array<double, 3> e{q[0] - length * d[0], q[1] - length * d[1], q[2] - length * d[2]};
  const auto altitude = [&](double s) {
    const double a = e[0] + s * d[0], b = e[1] + s * d[1], c = e[2] + s * d[2];
    return std::sqrt(a * a + b * b + c * c) - kR;
  };
  CHECK_NEAR(altitude(length), std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]) - kR, 1e-6);
  constexpr int panels = 20000;
  const double step = length / panels;
  double simpson = 0.0;
  for (int i = 0; i <= panels; ++i) {
    const double s = i * step;
    const double weight = (i == 0 || i == panels) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
    simpson += weight * table.extinction_coefficient<double>(lambda, altitude(s));
  }
  simpson *= step / 3.0;
  CHECK(simpson > 0.0);
  // The chord integral of obdeect (flat telescope depths: no in-telescope path).
  const TelescopeDepths depths = table.telescope_depths(kObs + kTelescopeZ, 0.0);
  const obdeect::propagation::SphereGeometry geometry{kObs, kR, 0.0, 0.0, kTelescopeZ};
  std::vector<double> p(1, 0.0);
  obdeect::propagation::atmospheric_acceptance_sphere<double>(chunk, table, depths, geometry, p);
  CHECK(p[0] > 0.0);
  CHECK_NEAR(-std::log(p[0]), simpson, 1e-9 * simpson);
}

// Rows of a wide mixture (both directions, |u| on both sides of the switch, zero paths, a row below
// the observation level, a wavelength above the table) give the same p for any row order and any
// split into ranges.
void test_order_and_range_invariance(const TransmissionTable &table) {
  std::vector<Row> rows;
  const double cos_values[] = {0.02, 0.09, 0.1, 0.4, 0.8, 1.0};
  for (std::size_t i = 0; i < 240; ++i) {
    const double cos_n = cos_values[i % 6];
    const bool upward = (i / 6) % 2 == 1;
    const double x = 1000.0 * std::sin(0.37 * static_cast<double>(i));
    const double y = 800.0 * std::cos(0.91 * static_cast<double>(i));
    const double z = 3.0 * static_cast<double>(i % 7) - 9.0;
    const double path = i % 11 == 0 ? 0.0 : 500.0 + 1500.0 * static_cast<double>(i % 37);
    Row row = direction_row(cos_n, upward, x, y, z, path, 300.0 + static_cast<double>(i % 40));
    row.zem = flat_emission(kObs + 200.0 + z, path, cos_n, upward);
    rows.push_back(row);
  }
  rows.push_back(
      direction_row(0.5, false, 0.0, 0.0, 0.0, 1000.0, 1200.0)); // wavelength above table
  rows.back().zem = 5000.0;
  rows.push_back(direction_row(0.5, false, 0.0, 0.0, 0.0, 1000.0, 355.0)); // below observation
  rows.back().zem = kObs - 50.0;
  const obdeect::propagation::SphereGeometry geometry{kObs, kR, 100.0, -50.0, 200.0};
  const TelescopeDepths depths = table.telescope_depths(kObs + 200.0, 20.0);
  const PhotonChunk chunk = sphere_chunk(rows);
  std::vector<double> whole(rows.size(), -1.0);
  obdeect::propagation::atmospheric_acceptance_sphere<double>(chunk, table, depths, geometry,
                                                              whole);
  std::vector<double> pieces(rows.size(), -1.0);
  const std::size_t cuts[] = {0, 1, 17, 100, 101, 200, rows.size()};
  for (std::size_t k = 6; k-- > 0;) {
    obdeect::propagation::atmospheric_acceptance_sphere<double>(chunk, table, depths, geometry,
                                                                pieces, cuts[k], cuts[k + 1]);
  }
  for (std::size_t r = 0; r < rows.size(); ++r) {
    CHECK(whole[r] == pieces[r]);
  }
  CHECK(whole.back() == 0.0);           // below the observation level
  CHECK(whole[rows.size() - 2] == 0.0); // wavelength above the table
  double positive = 0.0;
  for (const double value : whole) {
    positive += value > 0.0 ? 1.0 : 0.0;
  }
  CHECK(positive > 100.0);
}

// Input checks: the cz and path columns, finite zem, path >= 0, a unit direction.
void test_sphere_input_checks(const TransmissionTable &table) {
  const TelescopeDepths depths = table.telescope_depths(kObs + 200.0, 20.0);
  const obdeect::propagation::SphereGeometry geometry{kObs, 0.0, 0.0, 0.0, 200.0};
  std::vector<double> p(1, 0.0);
  const PhotonChunk good =
      sphere_chunk({direction_row(0.5, false, 0.0, 0.0, 0.0, 1000.0, kLambda)});
  {
    PhotonChunk no_cz = good;
    no_cz.columns[static_cast<std::size_t>(Column::cz)].clear();
    expect_invalid(
        [&] {
          obdeect::propagation::atmospheric_acceptance_sphere<double>(no_cz, table, depths,
                                                                      geometry, p);
        },
        "missing cz");
  }
  {
    PhotonChunk no_z = good;
    no_z.columns[static_cast<std::size_t>(Column::z)].clear();
    expect_invalid(
        [&] {
          obdeect::propagation::atmospheric_acceptance_sphere<double>(no_z, table, depths, geometry,
                                                                      p);
        },
        "missing z");
  }
  {
    PhotonChunk negative = sphere_chunk({direction_row(0.5, false, 0.0, 0.0, 0.0, -5.0, kLambda)});
    expect_invalid(
        [&] {
          obdeect::propagation::atmospheric_acceptance_sphere<double>(negative, table, depths,
                                                                      geometry, p);
        },
        "negative path");
  }
  {
    Row bad = direction_row(0.5, false, 0.0, 0.0, 0.0, 1000.0, kLambda);
    bad.cz = -0.9; // not a unit vector with cx
    expect_invalid(
        [&] {
          obdeect::propagation::atmospheric_acceptance_sphere<double>(sphere_chunk({bad}), table,
                                                                      depths, geometry, p);
        },
        "inconsistent cz");
  }
  expect_invalid(
      [&] {
        obdeect::propagation::atmospheric_acceptance_sphere<double>(
            good, table, depths, obdeect::propagation::SphereGeometry{kObs, -1.0, 0.0, 0.0, 0.0},
            p);
      },
      "negative Earth radius");
  expect_invalid(
      [&] {
        std::vector<double> short_p;
        obdeect::propagation::atmospheric_acceptance_sphere<double>(good, table, depths, geometry,
                                                                    short_p);
      },
      "acceptance size");
}

} // namespace

int main() {
  return run_checks([] {
    const TransmissionTable linear = site_table(WavelengthInterpolation::linear);
    const TransmissionTable floor_table = site_table(WavelengthInterpolation::floor);
    test_chord_integral_closed_forms();
    test_extinction_coefficient(linear);
    test_extinction_coefficient(floor_table);
    test_horizontal_flat_chord(linear);
    test_flat_chord_matches_plane_formula(linear);
    test_flat_chord_matches_plane_formula(floor_table);
    test_spherical_chord_dip(linear);
    test_continuity_at_switch(linear);
    test_continuity_at_switch(floor_table);
    test_sphere_rows_match_plane_rows(linear);
    test_sphere_rows_match_plane_rows(floor_table);
    test_sphere_row_chord_normalised(linear);
    test_order_and_range_invariance(linear);
    test_sphere_input_checks(linear);
  });
}
