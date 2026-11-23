// Upward photons (cz > 0) in the propagation stage: straight geometry of sky_to_ground, every loss
// category of DropCounts, bend_upward against a midpoint integration of Snell's law, reciprocity
// with bend_quadrature, the round trip through ground_to_sky, extinction of upward photons (flat
// and spherical), and independence of the output from chunking, row order and threads. See
// docs/propagation/refraction.md ("Upward photons") and docs/propagation/extinction.md.
#include "check.hpp"

#include "obdeect/propagation/extinction.hpp"
#include "obdeect/propagation/ground_to_sky.hpp"
#include "obdeect/propagation/pipeline.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/refraction_horizon.hpp"
#include "obdeect/propagation/sky.hpp"
#include "obdeect/propagation/timing.hpp"
#include "obdeect/propagation/transmission.hpp"
#include "obdeect/propagation/wavelength.hpp"

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using obdeect::propagation::AtmosphereProfile;
using obdeect::propagation::bend_quadrature;
using obdeect::propagation::bend_upward;
using obdeect::propagation::DropCounts;
using obdeect::propagation::Envelope;
using obdeect::propagation::GroundPlane;
using obdeect::propagation::RefractionTables;
using obdeect::propagation::SkyPlane;
using obdeect::propagation::TransmissionTable;
using obdeect::propagation::upward_turns_over;
using obdeect::stream::Column;
using obdeect::stream::PhotonChunk;

constexpr double kObs = 2000.0;                                 // observation level [m]
constexpr double kC = obdeect::stream::speed_of_light_m_per_ns; // [m/ns]
constexpr double kR = obdeect::stream::earth_radius_m;          // [m]

bool throws_invalid(auto &&action) {
  try {
    action();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

std::string data_file(const std::string &name) {
  return std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/" + name;
}

AtmosphereProfile ecmwf_profile() {
  return AtmosphereProfile::from_file(data_file("atmprof_ecmwf_north_winter_fixed.dat"));
}

// The ECMWF profile with the index n - 1 = nu set to a constant at every level (density and
// thickness unchanged, so the profile is valid): no bending.
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

PhotonChunk sky_chunk(const std::vector<std::array<double, 8>> &rows, bool with_cz) {
  // Each row: x y z cx cy cz t w (cz ignored without with_cz).
  PhotonChunk chunk;
  chunk.frame = obdeect::stream::Frame::sky;
  std::uint32_t mask =
      obdeect::stream::column_bit(Column::x) | obdeect::stream::column_bit(Column::y) |
      obdeect::stream::column_bit(Column::z) | obdeect::stream::column_bit(Column::cx) |
      obdeect::stream::column_bit(Column::cy) | obdeect::stream::column_bit(Column::t) |
      obdeect::stream::column_bit(Column::w);
  if (with_cz) {
    mask |= obdeect::stream::column_bit(Column::cz);
  }
  chunk.resize(rows.size(), mask);
  for (std::size_t r = 0; r < rows.size(); ++r) {
    const auto &row = rows[r];
    chunk[Column::x][r] = static_cast<float>(row[0]);
    chunk[Column::y][r] = static_cast<float>(row[1]);
    chunk[Column::z][r] = static_cast<float>(row[2]);
    chunk[Column::cx][r] = static_cast<float>(row[3]);
    chunk[Column::cy][r] = static_cast<float>(row[4]);
    if (with_cz) {
      chunk[Column::cz][r] = static_cast<float>(row[5]);
    }
    chunk[Column::t][r] = static_cast<float>(row[6]);
    chunk[Column::w][r] = static_cast<float>(row[7]);
  }
  return chunk;
}

// Direction cosine cz of the given sign from cx, cy (float32 storage, as in the chunks).
double cz_of(double cx, double cy, double sign) {
  const double fx = static_cast<float>(cx), fy = static_cast<float>(cy);
  const double cz = std::sqrt(1.0 - fx * fx - fy * fy);
  return sign * static_cast<float>(cz);
}

// Constant-index planes and tables (sky geometry helpers).
SkyPlane flat_plane(double telescope_z_m) {
  SkyPlane plane;
  plane.telescope_x_m = 5.0;
  plane.telescope_y_m = -3.0;
  plane.telescope_z_m = telescope_z_m;
  plane.observation_m = kObs;
  return plane;
}

// Straight geometry (no refraction, no profile): an upward photon 40 m above the observation level
// reaches a plane 100 m above it along its straight line; the time is the vacuum path.
void test_straight_geometry() {
  const double cx = 0.3, cy = -0.2;
  const double cz = cz_of(cx, cy, +1.0);
  auto chunk = sky_chunk({{10.0, -5.0, 40.0, cx, cy, cz, 7.0, 1.0}}, true);
  const SkyPlane plane = flat_plane(100.0);
  const DropCounts lost = obdeect::propagation::sky_to_ground(chunk, plane);
  CHECK(lost.total() == 0);
  CHECK(chunk.frame == obdeect::stream::Frame::ground);
  const double path = (100.0 - 40.0) / std::fabs(cz);
  CHECK_NEAR(chunk[Column::x][0], 10.0 + path * cx - 5.0, 1e-3);
  CHECK_NEAR(chunk[Column::y][0], -5.0 + path * cy + 3.0, 1e-3);
  CHECK_NEAR(chunk[Column::t][0], 7.0 + path / kC, 1e-3);
  CHECK_NEAR(chunk[Column::zem][0], kObs + 40.0, 1e-3);
  CHECK(chunk[Column::cx][0] == static_cast<float>(cx));
  CHECK(chunk[Column::cy][0] == static_cast<float>(cy));
  // The cz column is kept with the sign of the direction.
  CHECK(chunk.has(Column::cz) && chunk[Column::cz][0] > 0.0F);
  CHECK(chunk[Column::w][0] == 1.0F);
}

// Every loss category, with one row each (plane 1000 m above the observation level, tables built
// for the plane). A reaching upward row and a reaching downward row are kept.
void test_loss_categories() {
  const AtmosphereProfile profile = ecmwf_profile();
  const RefractionTables tables(profile, kObs + 1000.0);
  SkyPlane plane = flat_plane(1000.0);
  plane.tables = &tables;
  const double s_grazing = 0.99999; // turns over between 2 km and 3 km
  const std::vector<std::array<double, 8>> rows{
      {0.0, 0.0, 40.0, 1.0, 0.0, 0.0, 0.0, 1.0},                     // horizontal
      {0.0, 0.0, -5.0, 0.3, 0.1, cz_of(0.3, 0.1, +1.0), 0.0, 1.0},   // upward, below ground
      {0.0, 0.0, 1500.0, 0.3, 0.1, cz_of(0.3, 0.1, +1.0), 0.0, 1.0}, // upward above plane
      {0.0, 0.0, 500.0, 0.3, 0.1, cz_of(0.3, 0.1, -1.0), 0.0, 1.0},  // downward below plane
      {0.0, 0.0, -5.0, 0.3, 0.1, cz_of(0.3, 0.1, -1.0), 0.0, 1.0},   // downward below ground
      {0.0, 0.0, 0.0, s_grazing, 0.0, cz_of(s_grazing, 0.0, +1.0), 0.0, 1.0}, // turns over
      {0.0, 0.0, 0.0, 0.3, 0.1, cz_of(0.3, 0.1, +1.0), 0.0, 1.0},             // reaches (upward)
      {0.0, 0.0, 1200.0, 0.3, 0.1, cz_of(0.3, 0.1, -1.0), 0.0, 1.0},          // reaches (downward)
  };
  auto chunk = sky_chunk(rows, true);
  const DropCounts lost = obdeect::propagation::sky_to_ground(chunk, plane);
  CHECK(lost.horizontal == 1);
  CHECK(lost.below_ground == 2);
  CHECK(lost.moving_away == 2);
  CHECK(lost.turns_over == 1);
  CHECK(lost.not_bendable == 0);
  CHECK(lost.misses_level == 0);
  CHECK(lost.total() == 6);
  // Dropped rows have w = 0; the two reaching rows keep their weight.
  std::size_t kept = 0;
  for (std::size_t r = 0; r < chunk.rows(); ++r) {
    kept += chunk[Column::w][r] > 0.0F ? 1U : 0U;
  }
  CHECK(kept == 2);
  CHECK(chunk[Column::w][6] > 0.0F && chunk[Column::w][7] > 0.0F);
  // Dropped rows keep cz = -1 (as before for downward photons).
  CHECK(chunk[Column::cz][0] == -1.0F && chunk[Column::cz][5] == -1.0F);
  // Downward photons below the plane and upward photons above it are not confused.
  CHECK(chunk[Column::cz][7] < 0.0F && chunk[Column::cz][6] > 0.0F);

  // Spherical geometry: a downward photon from 17 km at zenith 89 degrees passes above the
  // telescope level (impact parameter at least R + h_t): misses_level.
  SkyPlane spherical = flat_plane(0.0);
  spherical.earth_radius_m = kR;
  spherical.telescope_altitude_m = kObs;
  const double cx = std::sin(89.0 * 3.14159265358979323846 / 180.0);
  auto far = sky_chunk({{0.0, 0.0, 17000.0, cx, 0.0, cz_of(cx, 0.0, -1.0), 0.0, 1.0}}, true);
  const DropCounts far_lost = obdeect::propagation::sky_to_ground(far, spherical);
  CHECK(far_lost.misses_level == 1);
  CHECK(far_lost.total() == 1);
}

// bend_upward with a constant index: no bending, the time is n L / (c cos theta_e).
void test_constant_index() {
  const double nu = 2.9e-4;
  const AtmosphereProfile profile = constant_index_profile(nu);
  const double h_e = kObs, h_t = kObs + 500.0;
  const double cx0 = 0.6, cy0 = 0.0;
  const double c_e = 0.8;
  double cx = cx0, cy = cy0, x = 3.0, y = -4.0, t = 1.0;
  CHECK(bend_upward(profile, h_e, h_t, cx, cy, x, y, t, 1.0, 1.0));
  CHECK_NEAR(x, 3.0, 1e-9);
  CHECK_NEAR(y, -4.0, 1e-9);
  CHECK_NEAR(cx, cx0, 1e-12);
  CHECK_NEAR(cy, cy0, 1e-12);
  CHECK_NEAR(t, 1.0 + (1.0 + nu) * (h_t - h_e) / (kC * c_e), 1e-6);
}

// Midpoint integration of Snell's law n sin(theta) = n_e sin(theta_e) from h_e up to h_t (metres),
// with n = 1 + nu from the profile and k = k_group = 1. Independent of the quadrature in
// refraction_horizon.hpp. Returns the displacement along (cx, cy)/s_e, the travel time [ns] and the
// direction scale s_t / s_e.
struct Snell {
  double displacement_m;
  double time_ns;
  double scale;
};

Snell snell_midpoint(const AtmosphereProfile &profile, double h_e, double h_t, double cx, double cy,
                     int steps) {
  const double s_e = std::sqrt(cx * cx + cy * cy);
  const double c_e = std::sqrt(1.0 - s_e * s_e);
  const double nu_e = profile.index_minus_one_cm(100.0 * h_e);
  const double dh = (h_t - h_e) / steps;
  double displacement = 0.0, excess = 0.0;
  for (int i = 0; i < steps; ++i) {
    const double h = h_e + (i + 0.5) * dh;
    const double nu = profile.index_minus_one_cm(100.0 * h);
    const double s = s_e * (1.0 + nu_e) / (1.0 + nu);
    const double c = std::sqrt(1.0 - s * s);
    displacement += (s / c - s_e / c_e) * dh;
    excess += ((1.0 + nu) / c - 1.0 / c_e) * dh;
  }
  const double nu_t = profile.index_minus_one_cm(100.0 * h_t);
  Snell out{};
  out.displacement_m = displacement;
  out.time_ns = ((h_t - h_e) / c_e + excess) / kC;
  out.scale = s_e * (1.0 + nu_e) / (1.0 + nu_t) / s_e;
  return out;
}

// bend_upward against the midpoint integration of Snell's law (1e-6 m, 1e-6 ns).
void test_bend_upward_matches_snell() {
  const AtmosphereProfile profile = ecmwf_profile();
  const double h_e = kObs, h_t = kObs + 1500.0;
  const double directions[][2] = {{0.3, 0.1}, {0.9, 0.2}, {0.99, 0.05}};
  for (const auto &direction : directions) {
    const double cx0 = direction[0], cy0 = direction[1];
    const Snell reference = snell_midpoint(profile, h_e, h_t, cx0, cy0, 400000);
    double cx = cx0, cy = cy0, x = 0.0, y = 0.0, t = 0.0;
    CHECK(bend_upward(profile, h_e, h_t, cx, cy, x, y, t, 1.0, 1.0));
    const double s_e = std::sqrt(cx0 * cx0 + cy0 * cy0);
    CHECK_NEAR(x, reference.displacement_m * cx0 / s_e, 1e-6);
    CHECK_NEAR(y, reference.displacement_m * cy0 / s_e, 1e-6);
    CHECK_NEAR(t, reference.time_ns, 1e-6);
    CHECK_NEAR(cx, cx0 * reference.scale, 1e-9);
    CHECK_NEAR(cy, cy0 * reference.scale, 1e-9);
    // The bent ray leaves towards the horizon: sin theta grows with height (n decreases).
    CHECK(cx * cx + cy * cy > cx0 * cx0 + cy0 * cy0);
  }
}

// Reciprocity: the exit ray of bend_upward, reversed and bent downwards by bend_quadrature from
// the tables at the emission level, returns to the start (position, time and direction).
void test_reciprocity() {
  const AtmosphereProfile profile = ecmwf_profile();
  const double h_e = kObs, h_t = kObs + 1500.0, length = h_t - h_e;
  for (const double cx0 : {0.3, 0.9}) {
    const double cy0 = 0.2;
    const double c_e = std::sqrt(1.0 - cx0 * cx0 - cy0 * cy0);
    // Straight-line position at the plane of a photon emitted at the origin.
    double fx = length * cx0 / c_e, fy = length * cy0 / c_e, fcx = cx0, fcy = cy0, ft = 0.0;
    CHECK(bend_upward(profile, h_e, h_t, fcx, fcy, fx, fy, ft, 1.0, 1.0));
    const double c_f = std::sqrt(1.0 - fcx * fcx - fcy * fcy);
    // Reverse: from the plane point downwards; its straight line reaches the emission level at
    // the back-projection (x_f - L cx_f / c_f). The time at emission is -t_f.
    double rx = fx - length * fcx / c_f, ry = fy - length * fcy / c_f;
    double rcx = -fcx, rcy = -fcy, rt = -ft;
    const RefractionTables tables(profile, h_e);
    CHECK(bend_quadrature(tables, h_t, rcx, rcy, rx, ry, rt));
    CHECK_NEAR(rx, 0.0, 1e-6);
    CHECK_NEAR(ry, 0.0, 1e-6);
    CHECK_NEAR(rt, 0.0, 1e-6);
    CHECK_NEAR(rcx, -cx0, 1e-9);
    CHECK_NEAR(rcy, -cy0, 1e-9);
  }
}

// A near-horizontal photon turns over before a plane 1 km above it: bend_upward returns false and
// leaves its arguments unchanged; upward_turns_over says so. A steep photon does not turn over.
void test_turning_over() {
  const AtmosphereProfile profile = ecmwf_profile();
  const double h_e = kObs, h_t = kObs + 1000.0;
  double cx = 0.99999, cy = 0.0, x = 1.0, y = 2.0, t = 3.0;
  CHECK(upward_turns_over(profile, h_e, h_t, cx, cy, 1.0));
  CHECK(!bend_upward(profile, h_e, h_t, cx, cy, x, y, t, 1.0, 1.0));
  CHECK(cx == 0.99999 && cy == 0.0 && x == 1.0 && y == 2.0 && t == 3.0);
  double steep_x = 0.3, steep_y = 0.0;
  CHECK(!upward_turns_over(profile, h_e, h_t, steep_x, steep_y, 1.0));
  CHECK(bend_upward(profile, h_e, h_t, steep_x, steep_y, x, y, t, 1.0, 1.0));
  // Invalid inputs are not bendable (false), and do not count as turning over.
  double bad_x = 1.0, bad_y = 0.0;
  CHECK(!upward_turns_over(profile, h_e, h_t, bad_x, bad_y, 1.0));
  CHECK(!bend_upward(profile, h_t, h_e, bad_x, bad_y, x, y, t, 1.0, 1.0)); // emission above plane
}

// Round trip sky -> G -> sky for upward photons (profile, refraction applied: no bending).
void test_round_trip() {
  const AtmosphereProfile profile = ecmwf_profile();
  SkyPlane plane = flat_plane(100.0);
  plane.profile = &profile;
  const std::vector<std::array<double, 8>> rows{
      {12.0, -7.0, 0.0, 0.2, 0.3, cz_of(0.2, 0.3, +1.0), 5.0, 1.0},
      {-30.0, 40.0, 40.0, 0.5, -0.4, cz_of(0.5, -0.4, +1.0), 9.0, 1.0},
      {100.0, 0.0, 100.0, 0.1, 0.0, cz_of(0.1, 0.0, +1.0), 1.5, 1.0},
  };
  auto chunk = sky_chunk(rows, true);
  const PhotonChunk original = chunk;
  const DropCounts lost = obdeect::propagation::sky_to_ground(chunk, plane);
  CHECK(lost.total() == 0);
  const GroundPlane ground{plane.telescope_x_m, plane.telescope_y_m, plane.telescope_z_m, kObs,
                           &profile};
  std::uint64_t dropped = 0;
  const auto pieces = obdeect::propagation::ground_to_sky(chunk, ground, dropped);
  CHECK(dropped == 0);
  CHECK(pieces.size() == 1);
  if (pieces.size() != 1) {
    return;
  }
  const PhotonChunk &back = pieces[0];
  for (std::size_t r = 0; r < rows.size(); ++r) {
    CHECK_NEAR(back[Column::x][r], original[Column::x][r], 1e-3);
    CHECK_NEAR(back[Column::y][r], original[Column::y][r], 1e-3);
    CHECK_NEAR(back[Column::z][r], original[Column::z][r], 1e-3);
    CHECK_NEAR(back[Column::t][r], original[Column::t][r], 1e-3);
    CHECK(back[Column::cz][r] > 0.0F);
  }
  // An upward G photon emitted above the telescope plane cannot come from below: rejected.
  auto above = chunk;
  above[Column::zem][0] = static_cast<float>(kObs + 150.0);
  std::uint64_t none = 0;
  CHECK(throws_invalid([&] { obdeect::propagation::ground_to_sky(above, ground, none); }));
}

// Chunking, row order and the per-row kernel: a mixed-direction chunk gives the same rows for one
// chunk, chunks of 7 and 1 rows, and reversed order. Extinction with 1 and 4 threads agrees too.
void test_mixed_directions_independent() {
  const AtmosphereProfile profile = ecmwf_profile();
  const RefractionTables tables(profile, kObs + 100.0);
  SkyPlane plane = flat_plane(100.0);
  plane.tables = &tables;
  std::vector<std::array<double, 8>> rows;
  for (std::size_t i = 0; i < 300; ++i) {
    const auto k = static_cast<double>(i);
    const double cx = i % 9 == 0 ? 0.95 : 0.5 * std::sin(0.01 * k);
    const double cy = 0.3 * std::cos(0.013 * k);
    const double sign = i % 3 == 0 ? -1.0 : 1.0;
    double z = 3000.0 + 20.0 * static_cast<double>(i % 11);
    if (i % 3 != 0) {
      // Upward rows: below the plane (some below ground), or above it (moving away).
      z = i % 5 == 0 ? 150.0 : static_cast<double>(i % 7) * 15.0 - 5.0;
    }
    rows.push_back({1000.0 * std::sin(0.1 * k), 1000.0 * std::cos(0.07 * k), z, cx, cy,
                    cz_of(cx, cy, sign), 0.5 * k, 1.0 + static_cast<double>(i % 4)});
  }
  rows.push_back({0.0, 0.0, 10.0, 1.0, 0.0, 0.0, 0.0, 1.0}); // horizontal
  auto whole = sky_chunk(rows, true);
  const DropCounts counts = obdeect::propagation::sky_to_ground(whole, plane);
  CHECK(counts.total() > 0);
  CHECK(counts.horizontal == 1);
  CHECK(counts.moving_away > 0 && counts.below_ground > 0);

  const auto same_row = [](const PhotonChunk &a, std::size_t i, const PhotonChunk &b,
                           std::size_t j) {
    for (const Column c : {Column::x, Column::y, Column::cx, Column::cy, Column::t, Column::zem,
                           Column::w, Column::cz}) {
      if (a[c][i] != b[c][j]) {
        return false;
      }
    }
    return true;
  };
  // Chunks of 7 rows (each converted on its own).
  for (std::size_t first = 0; first < rows.size(); first += 7) {
    const std::size_t end = std::min(rows.size(), first + 7);
    std::vector<std::array<double, 8>> piece(rows.begin() + static_cast<std::ptrdiff_t>(first),
                                             rows.begin() + static_cast<std::ptrdiff_t>(end));
    auto chunk = sky_chunk(piece, true);
    obdeect::propagation::sky_to_ground(chunk, plane);
    for (std::size_t r = 0; r < piece.size(); ++r) {
      CHECK(same_row(chunk, r, whole, first + r));
    }
  }
  // Reversed order.
  std::vector<std::array<double, 8>> reversed(rows.rbegin(), rows.rend());
  auto back = sky_chunk(reversed, true);
  obdeect::propagation::sky_to_ground(back, plane);
  for (std::size_t r = 0; r < rows.size(); ++r) {
    CHECK(same_row(back, r, whole, rows.size() - 1 - r));
  }

  // Extinction of the converted chunk with 1 and 4 threads (weighted, identical w and p).
  const auto table =
      TransmissionTable::from_file(data_file("atm_trans_2000_1_10_0_0_2000.dat"), kObs,
                                   obdeect::propagation::WavelengthInterpolation::linear);
  const auto depths = table.telescope_depths(kObs + 100.0, 0.0);
  const Envelope envelope({200.0, 800.0}, {0.9, 0.3});
  auto weights = [&](unsigned threads) {
    auto chunk = whole;
    obdeect::propagation::ExtinctionSettings settings;
    settings.table = &table;
    settings.depths = &depths;
    settings.envelope = &envelope;
    settings.wavelengths = {250.0, 700.0};
    settings.seed = 7;
    settings.event = 1;
    settings.telescope = 1;
    settings.weighted = true;
    settings.threads = threads;
    const auto out = obdeect::propagation::apply_extinction_chunk(std::move(chunk), settings, 0);
    std::vector<float> w;
    for (const auto &piece : out) {
      w.insert(w.end(), piece[Column::w].begin(), piece[Column::w].end());
    }
    return w;
  };
  const auto one = weights(1);
  const auto four = weights(4);
  CHECK(one.size() == rows.size());
  CHECK(one == four);
  // Dropped rows have zero weight in the input, so zero after extinction too.
  for (std::size_t r = 0; r < rows.size() && r < one.size(); ++r) {
    if (whole[Column::w][r] == 0.0F) {
      CHECK(one[r] == 0.0F);
    }
  }
}

// Flat extinction of an upward photon: transmission exp(-|tau(h_t) - tau(zem)| / |cz|), with the
// optical depths from the table at the two heights (independent of the transmission code).
void test_flat_extinction_upward() {
  const auto table =
      TransmissionTable::from_file(data_file("atm_trans_2000_1_10_0_0_2000.dat"), kObs,
                                   obdeect::propagation::WavelengthInterpolation::floor);
  const double h_t = kObs + 100.0;
  const double zem = kObs + 30.0;
  const double cx = std::sqrt(0.75); // cz = 0.5, airmass 2
  const double cz = cz_of(cx, 0.0, +1.0);
  const double lambda = 400.0;
  const double expected = std::exp(-std::fabs(table.optical_depth<double>(lambda, h_t) -
                                              table.optical_depth<double>(lambda, zem)) /
                                   std::fabs(cz));
  const auto depths = table.telescope_depths(h_t, 0.0);
  const double od_tel = table.telescope_depths_at(depths, lambda).od_tel;
  CHECK_NEAR(table.transmission2<double>(lambda, zem, 1.0 / std::fabs(cz), od_tel, 0.0), expected,
             1e-12);

  // The same through the chunk kernel (cz column present, upward row). The kernel reads cx in
  // float32, so the airmass differs from the double value by about 1e-9.
  PhotonChunk chunk;
  chunk.frame = obdeect::stream::Frame::ground;
  const std::uint32_t mask =
      obdeect::stream::column_bit(Column::cx) | obdeect::stream::column_bit(Column::cy) |
      obdeect::stream::column_bit(Column::cz) | obdeect::stream::column_bit(Column::zem) |
      obdeect::stream::column_bit(Column::wavelength) | obdeect::stream::column_bit(Column::x) |
      obdeect::stream::column_bit(Column::y) | obdeect::stream::column_bit(Column::t) |
      obdeect::stream::column_bit(Column::w);
  chunk.resize(1, mask);
  chunk[Column::cx][0] = static_cast<float>(cx);
  chunk[Column::cz][0] = static_cast<float>(cz);
  chunk[Column::zem][0] = static_cast<float>(zem);
  chunk[Column::wavelength][0] = static_cast<float>(lambda);
  std::vector<double> p(1, 0.0);
  obdeect::propagation::atmospheric_acceptance<double>(chunk, table, depths, p);
  CHECK_NEAR(p[0], expected, 1e-7);
}

// Spherical extinction of an upward photon: the slant depth over [zem, h_t] against a fine
// integral of |d tau| / cos theta with tau from the table (cos theta(h) = sqrt(1 - p^2/(R+h)^2),
// p = (R + h_t) sin theta_t), and the plane-parallel limit R -> infinity.
void test_curved_extinction_upward() {
  const auto table =
      TransmissionTable::from_file(data_file("atm_trans_2000_1_10_0_0_2000.dat"), kObs,
                                   obdeect::propagation::WavelengthInterpolation::floor);
  const double lambda = 400.0;
  const double h_t = 2100.0, zem = 2000.0, u = 0.5;
  const double p = (kR + h_t) * std::sqrt(1.0 - u * u);
  constexpr int steps = 40000;
  double fine = 0.0;
  double h_prev = zem;
  double tau_prev = table.optical_depth<double>(lambda, zem);
  for (int i = 1; i <= steps; ++i) {
    const double h = zem + (h_t - zem) * i / steps;
    const double tau = table.optical_depth<double>(lambda, h);
    const double mid = 0.5 * (h + h_prev);
    const double cos_mid = std::sqrt(1.0 - p * p / ((kR + mid) * (kR + mid)));
    fine += std::fabs(tau - tau_prev) / cos_mid;
    h_prev = h;
    tau_prev = tau;
  }
  const double slant = table.slant_optical_depth_upward<double>(lambda, h_t, zem, u, kR);
  CHECK(slant > 0.0);
  CHECK_NEAR(slant, fine, 1e-4 * fine);
  CHECK_NEAR(table.transmission_curved_upward<double>(lambda, zem, u, h_t, 0.0, kR),
             std::exp(-slant), 1e-12);

  // Downward photons keep zero slant depth below the telescope (transmission exp(-od_focus)).
  CHECK(table.slant_optical_depth<double>(lambda, h_t, zem, u, kR) == 0.0);
  CHECK_NEAR(table.transmission_curved<double>(lambda, zem, u, h_t, 0.0, kR), 1.0, 1e-15);

  // R -> infinity: the plane-parallel slant depth |tau(h_t) - tau(zem)| / cos theta_t.
  const double flat = std::fabs(table.optical_depth<double>(lambda, h_t) -
                                table.optical_depth<double>(lambda, zem)) /
                      u;
  CHECK_NEAR(table.slant_optical_depth_upward<double>(lambda, h_t, zem, u, 1e12), flat,
             1e-7 * flat);
}

} // namespace

int main() {
  return run_checks([] {
    test_straight_geometry();
    test_loss_categories();
    test_constant_index();
    test_bend_upward_matches_snell();
    test_reciprocity();
    test_turning_over();
    test_round_trip();
    test_mixed_directions_independent();
    test_flat_extinction_upward();
    test_curved_extinction_upward();
  });
}
