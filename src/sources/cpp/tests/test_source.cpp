#include "check.hpp"
#include "source_fixtures.hpp"

#include "obdeect/sources/emit.hpp"
#include "obdeect/sources/source.hpp"

#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace obdeect::sources;
using obdeect::stream::Column;

namespace {

using obdeect_test::accepts;

constexpr double inf = std::numeric_limits<double>::infinity();

// All photons of a test source, one array per column. A PhotonChunk holds at most 2^16
// rows, so the test keeps the columns itself.
struct Photons {
  std::array<std::vector<float>, obdeect::stream::n_columns> columns{};

  std::size_t rows() const { return columns[static_cast<std::size_t>(Column::x)].size(); }
  const std::vector<float> &operator[](Column c) const {
    return columns[static_cast<std::size_t>(c)];
  }
};

// Emit candidates [0, n) of a source in chunks of `chunk` candidates. Rows dropped by emit()
// (line photons that miss their target) are not stored.
Photons emit_all(const Source &source, std::uint64_t seed, std::uint64_t n,
                 std::size_t chunk = 1 << 16, std::uint32_t event = 1,
                 std::uint32_t telescope = 2) {
  Photons all;
  obdeect::stream::PhotonChunk part;
  for (std::uint64_t first = 0; first < n; first += chunk) {
    const auto candidates = static_cast<std::size_t>(std::min<std::uint64_t>(chunk, n - first));
    const std::size_t rows = emit(source, seed, event, telescope, first, candidates, part);
    for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c) {
      if (!(source_columns & (std::uint32_t{1} << c)))
        continue;
      const auto &src = part.columns[c];
      all.columns[c].insert(all.columns[c].end(), src.begin(),
                            src.begin() + static_cast<std::ptrdiff_t>(rows));
    }
  }
  return all;
}

bool same_chunk(const Photons &a, const Photons &b) {
  if (a.rows() != b.rows())
    return false;
  for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c)
    if (a.columns[c] != b.columns[c])
      return false;
  return true;
}

Source star_source(double distance, double photons = 1000.0) {
  Source s;
  s.photons = photons;
  s.geometry = PlaneWave{normalised({0.01, 0.0, 1.0}), distance, 5.0, 0.0, 30.0};
  return s;
}

Source lamp_source(double photons = 1e5) {
  Source s;
  s.photons = photons;
  s.geometry = Point{{0.0, 0.0, 10.0}, {0.0, 0.0, -1.0}, Isotropic{}};
  s.pulse = UniformRange{-1.0, 1.0};
  return s;
}

// Position and direction of a photon at the plane z = z0 (for rays from the emission point).
Vec3 at_plane(const Photons &c, std::size_t i, double z0) {
  const double x = c[Column::x][i], y = c[Column::y][i], z = c[Column::z][i];
  const double cx = c[Column::cx][i], cy = c[Column::cy][i], cz = c[Column::cz][i];
  const double s = (z0 - z) / cz;
  return {x + s * cx, y + s * cy, z0};
}

void stars() {
  // Plane wave over a 5 m disc (obdeect ray tracing test_kernels.cpp sources): all rays start at
  // the start plane z = 30 m, a quarter of the photons land within 2.5 m of the centre.
  const auto c = emit_all(star_source(inf, 100000.0), 3, 100000);
  double inside = 0.0;
  for (std::size_t i = 0; i < c.rows(); ++i) {
    CHECK_NEAR(c[Column::z][i], 30.0, 1e-4);
    const Vec3 target = at_plane(c, i, 0.0);
    CHECK(std::hypot(target.x, target.y) <= 5.0 + 1e-3);
    inside += std::hypot(target.x, target.y) < 2.5;
  }
  CHECK_NEAR(inside / 100000.0, 0.25, 0.005);
  // Times are the arrival of the wavefront through the disc centre: t = dot(start, dir) / c.
  const Vec3 towards = normalised({0.01, 0.0, 1.0});
  for (std::size_t i = 0; i < 100; ++i) {
    const Vec3 start{c[Column::x][i], c[Column::y][i], c[Column::z][i]};
    const Vec3 dir{c[Column::cx][i], c[Column::cy][i], c[Column::cz][i]};
    CHECK_NEAR(dot(dir, towards), -1.0, 1e-6);
    CHECK_NEAR(c[Column::t][i], dot(start, dir) / speed_of_light_m_per_ns, 1e-3);
  }
  // Point source at 20 m, 1.5 m above the entrance plane, below the start plane:
  // photons leave the source itself (obdeect ray tracing PointSource, finite distance).
  Source lamp;
  lamp.photons = 50.0;
  lamp.geometry = PlaneWave{normalised({0.6, 0.0, 0.8}), 20.0, 1.0, 1.5, 30.0};
  // The source sits at (0.6, 0, 0.8) * 20 m above the disc centre (0, 0, 1.5 m): (12, 0, 17.5).
  const auto p = emit_all(lamp, 3, 50);
  for (std::size_t i = 0; i < p.rows(); ++i) {
    CHECK_NEAR(p[Column::x][i], 12.0, 1e-4);
    CHECK_NEAR(p[Column::y][i], 0.0, 1e-4);
    CHECK_NEAR(p[Column::z][i], 17.5, 1e-4);
    // Start at the source, which is 20 m before the disc-centre wavefront: t = -d / c.
    CHECK_NEAR(p[Column::t][i], -20.0 / speed_of_light_m_per_ns, 1e-3);
  }
}

// Fraction of the light of a point that reaches a disc, by a midpoint rule over the disc area with
// the same integrand as disc_fraction (1 / (2 pi I) int I(u) h / d^3 dA).
double disc_fraction_midpoint(const Point &p, const DiscTarget &disc, int steps) {
  const Basis basis = orthonormal_basis(disc.normal);
  const double h = std::abs(dot(disc.centre - p.position, disc.normal));
  double sum = 0.0;
  for (int i = 0; i < steps; ++i) {
    const double rho = disc.radius_m * (i + 0.5) / steps;
    for (int k = 0; k < steps; ++k) {
      const double phi = two_pi * (k + 0.5) / steps;
      const Vec3 v = disc.centre + basis.e1 * (rho * std::cos(phi)) +
                     basis.e2 * (rho * std::sin(phi)) - p.position;
      const double d = norm(v);
      const double cosine = std::clamp(dot(v, p.axis) / d, -1.0, 1.0);
      sum += std::visit([&](const auto &a) { return a.intensity(1.0 - cosine); }, p.angular) * rho *
             h / (d * d * d);
    }
  }
  const double dr = disc.radius_m / steps, dphi = two_pi / steps;
  return sum * dr * dphi / (two_pi * integral(p.angular));
}

void disc_fraction_closed_forms() {
  // On the axis: the cdf of the profile at 1 - cos(atan(R / h)). Isotropic: solid angle / 4 pi.
  const double h = 10.0, r = 2.0;
  const double solid_angle = 2.0 * pi * (1.0 - h / std::sqrt(h * h + r * r));
  const DiscTarget disc{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, r};
  const Point iso{{0.0, 0.0, h}, {0.0, 0.0, -1.0}, Isotropic{}};
  CHECK_NEAR(disc_fraction(iso, disc), solid_angle / (4.0 * pi), 1e-15);
  // Lambertian: sin^2(alpha) = R^2 / (R^2 + h^2).
  const Point lamb{{0.0, 0.0, h}, {0.0, 0.0, -1.0}, Lambertian{}};
  CHECK_NEAR(disc_fraction(lamb, disc), r * r / (r * r + h * h), 1e-14);
  // Filled cone of 5 deg inside the disc's cap (alpha = 11.3 deg): all of its light.
  const Point cone{{0.0, 0.0, h}, {0.0, 0.0, -1.0}, FilledCone{1.0 - std::cos(5.0 * pi / 180.0)}};
  CHECK_NEAR(disc_fraction(cone, disc), 1.0, 1e-14);
  // Off the axis of the profile: the disc is not the cap about the axis; a 5 deg cone there does
  // not reach a disc 60 deg away (exact zero).
  const Point tilted{{0.0, 0.0, h},
                     normalised({std::sin(1.0), 0.0, -std::cos(1.0)}),
                     FilledCone{1.0 - std::cos(5.0 * pi / 180.0)}};
  CHECK(disc_fraction(tilted, disc) == 0.0);
}

void disc_fraction_quadrature() {
  // Off-axis disc (tilted normal, centre off the point's axis): the Gauss-Legendre quadrature
  // (64 x 64) against a midpoint rule of 2000 x 2000 cells, for smooth profiles, to 1e-9 relative.
  const DiscTarget disc{{1.0, -0.5, 3.0}, normalised({0.2, 0.1, 1.0}), 4.0};
  const Point iso{{0.0, 0.0, 0.0}, normalised({0.05, 0.0, 1.0}), Isotropic{}};
  // The midpoint rule is O(h^2): Richardson extrapolation of the 2000 and 4000 cell rules.
  const auto richardson = [&](const Point &p) {
    return (4.0 * disc_fraction_midpoint(p, disc, 4000) - disc_fraction_midpoint(p, disc, 2000)) /
           3.0;
  };
  const double gl = disc_fraction(iso, disc);
  const double mid = richardson(iso);
  CHECK(gl > 0.0 && gl < 0.5);
  CHECK(std::abs(gl - mid) <= 1e-9 * mid);
  // A Gaussian profile (table, sigma 0.1 rad) on the same disc.
  const Point gauss{
      {0.0, 0.0, 0.0}, normalised({0.05, 0.0, 1.0}), TabulatedAngular{gaussian_angular_table(0.1)}};
  const double g_gl = disc_fraction(gauss, disc);
  // The axis (5 sigma of the table, 28 deg) lies inside the disc: all of the light reaches it.
  CHECK(g_gl > 0.99 && g_gl <= 1.0 + 1e-9);
  // The table is piecewise linear (kinks every 0.05 sigma), so the rules agree to 1e-5 only.
  const double g_mid = (4.0 * disc_fraction_midpoint(gauss, disc, 4000) -
                        disc_fraction_midpoint(gauss, disc, 2000)) /
                       3.0;
  CHECK(std::abs(g_gl - g_mid) <= 1e-5 * g_gl);
}

void emitters_isotropic_disc() {
  // Isotropic emitter 10 m above a 2 m disc (obdeect ray tracing test_kernels.cpp emitters): the
  // photons that reach the disc are its solid angle over 4 pi: 2 pi (1 - h / sqrt(h^2 + R^2)) / 4
  // pi.
  Source s = lamp_source(4e5);
  s.target = DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0};
  const double h = 10.0, r = 2.0;
  const double fraction = (1.0 - h / std::sqrt(h * h + r * r)) / 2.0;
  CHECK_NEAR(expected_count(s), 4e5 * fraction, 1e-9);
  CHECK_NEAR(expected_weight(s), expected_count(s), 1e-9);
  const auto c = emit_all(s, 5, 4000);
  CHECK(c.rows() == 4000);
  for (std::size_t i = 0; i < c.rows(); ++i) {
    // Every candidate is a hit: weight 1, on the disc, moving down onto it.
    CHECK(c[Column::w][i] == 1.0f);
    const Vec3 hit = at_plane(c, i, 0.0);
    CHECK(std::hypot(hit.x, hit.y) <= 2.0 + 1e-3);
    CHECK(c[Column::cz][i] < 0.0f);
  }
  CHECK(same_chunk(c, emit_all(s, 5, 4000, 333)));
}

void disc_target_angular_profile() {
  // A 5 deg filled cone lights only the inner part of the disc (2 m radius at 10 m: the disc
  // subtends 11.3 deg): every candidate is inside the cone, with weight 1.
  Source s = lamp_source(1e5);
  s.geometry =
      Point{{0.0, 0.0, 10.0}, {0.0, 0.0, -1.0}, FilledCone{1.0 - std::cos(5.0 * pi / 180.0)}};
  s.target = DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0};
  const auto c = emit_all(s, 5, 10000);
  CHECK(c.rows() == 10000);
  CHECK_NEAR(expected_weight(s), expected_count(s), 1e-9);
  for (std::size_t i = 0; i < c.rows(); ++i) {
    const double angle = std::acos(-c[Column::cz][i]) * 180.0 / pi;
    CHECK(angle <= 5.0 + 1e-3);
    CHECK(c[Column::w][i] == 1.0f);
  }
}

// A narrow beam seen off the disc axis: the acceptance is integrated over the beam, not over a
// fixed grid of the disc area that could miss its footprint (review of PR 83: the whole 1e-8 cone
// reaches the radius-2 disc, but a disc-area grid gave 0 and the emission was empty).
void disc_target_narrow_off_axis_beam() {
  Source s = lamp_source(1000);
  s.geometry = Point{{0.0, 0.0, 10.0}, normalised({0.1, 0.0, -1.0}), FilledCone{1e-8}};
  s.target = DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0};
  CHECK_NEAR(disc_fraction(std::get<Point>(s.geometry), std::get<DiscTarget>(s.target)), 1.0, 1e-9);
  CHECK_NEAR(expected_count(s), 1000.0, 1e-6);
  const auto c = emit_all(s, 3, 1000);
  CHECK(c.rows() == 1000);
  for (std::size_t i = 0; i < c.rows(); ++i)
    CHECK(c[Column::w][i] == 1.0f);
}

// A broad profile and a small, distant disc: the directions are drawn in the cone enclosing the
// disc (rejection on the profile), so a valid run needs few attempts per photon (review of PR 83:
// drawing from the full Lambertian profile needed ~1e8 attempts per row and aborted). The disc is
// off the axis (quadrature branch): F = cos^2(theta) R^2 / d^2 for a small disc facing up.
void disc_target_small_distant_disc() {
  const double h = 10000.0, x0 = 50.0, r = 1.0;
  Source s = lamp_source(1.0);
  s.geometry = Point{{0.0, 0.0, h}, {0.0, 0.0, -1.0}, Lambertian{}};
  s.target = DiscTarget{{x0, 0.0, 0.0}, {0.0, 0.0, 1.0}, r};
  const double d2 = h * h + x0 * x0, cos_theta = h / std::sqrt(d2);
  const double expected = cos_theta * cos_theta * r * r / d2;
  const double fraction =
      disc_fraction(std::get<Point>(s.geometry), std::get<DiscTarget>(s.target));
  CHECK(std::abs(fraction / expected - 1.0) < 1e-4);
  s.photons = 1000.0 / fraction; // about 1e11 photons for 1000 rows
  const auto c = emit_all(s, 1, static_cast<std::uint64_t>(std::llround(expected_count(s))));
  CHECK(c.rows() >= 990);
  // Every photon reaches the disc.
  std::size_t outside = 0;
  for (std::size_t i = 0; i < c.rows(); ++i) {
    const Vec3 q = at_plane(c, i, 0.0);
    outside += (q.x - x0) * (q.x - x0) + q.y * q.y > r * r * (1.0 + 1e-6) ? 1 : 0;
  }
  CHECK(outside == 0);
}

void chunks_and_determinism() {
  Source s = lamp_source(1e6);
  s.pulse = ComplexPulse{1.0, 0.5, 0.5, 2.0};
  s.spectrum = InverseSquare{300.0, 600.0};
  s.target = DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0};
  const auto whole = emit_all(s, 9, 3000, 1 << 16);
  CHECK(same_chunk(whole, emit_all(s, 9, 3000, 1000)));
  CHECK(same_chunk(whole, emit_all(s, 9, 3000, 7)));
  CHECK(!same_chunk(whole, emit_all(s, 10, 3000, 1 << 16)));   // Another seed differs.
  CHECK(!same_chunk(whole, emit_all(s, 9, 3000, 1 << 16, 2))); // Another event differs.
  // A chunk starting at photon 2^20 is the same photons whatever the chunk length.
  obdeect::stream::PhotonChunk late, check;
  const std::size_t rows_late = emit(s, 9, 1, 2, 1U << 20, 100, late);
  const std::size_t rows_check = emit(s, 9, 1, 2, 1U << 20, 50, check);
  CHECK(rows_late == 100 && rows_check == 50);
  for (std::size_t i = 0; i < 50; ++i)
    CHECK(late[Column::x][i] == check[Column::x][i] && late[Column::t][i] == check[Column::t][i]);
  bool rejected = false;
  try {
    obdeect::stream::PhotonChunk big;
    emit(s, 9, 1, 2, 0, (1U << 16) + 1, big);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
}

void directions_and_angles() {
  // Directions are unit vectors for every profile; positions are the emission point's.
  const std::vector<Angular> profiles{Isotropic{}, FixedAngle{1.0 - std::cos(0.2)}, FilledCone{0.3},
                                      Lambertian{}, TabulatedAngular{gaussian_angular_table(0.2)}};
  for (const auto &profile : profiles) {
    Source s = lamp_source(1e4);
    s.geometry = Point{{1.0, -2.0, 3.0}, normalised({0.2, 0.1, -1.0}), profile};
    const auto c = emit_all(s, 1, 2000, 512);
    for (std::size_t i = 0; i < c.rows(); ++i) {
      const double length = std::hypot(c[Column::cx][i], c[Column::cy][i], c[Column::cz][i]);
      CHECK_NEAR(length, 1.0, 1e-6);
      CHECK_NEAR(c[Column::x][i], 1.0, 0.0);
      CHECK_NEAR(c[Column::y][i], -2.0, 0.0);
      CHECK_NEAR(c[Column::z][i], 3.0, 0.0);
    }
  }
  // Fixed-angle (option "Cone:") emission is at exactly the given angle.
  Source cone = lamp_source(1e3);
  const Vec3 axis = normalised({0.2, 0.1, -1.0});
  cone.geometry = Point{{0.0, 0.0, 10.0}, axis, FixedAngle{1.0 - std::cos(0.2)}};
  const auto c = emit_all(cone, 4, 500);
  for (std::size_t i = 0; i < c.rows(); ++i) {
    const double cosine =
        c[Column::cx][i] * axis.x + c[Column::cy][i] * axis.y + c[Column::cz][i] * axis.z;
    CHECK_NEAR(cosine, std::cos(0.2), 1e-6);
  }
  // A plane wave keeps its start-plane height (30 m) in the emitted positions.
  const auto s = emit_all(star_source(inf, 10.0), 1, 10);
  CHECK(s.rows() == 10);
  CHECK_NEAR(s[Column::z][0], 30.0, 1e-4);
}

void counts() {
  // Counts: floor of the mean plus one with the fractional probability.
  const Source s = star_source(inf, 2.5);
  double total = 0.0;
  const int trials = 40000;
  int outside = 0;
  for (int seed = 0; seed < trials; ++seed) {
    const auto k = photon_count(s, static_cast<std::uint64_t>(seed), 0, 0);
    outside += (k == 2 || k == 3) ? 0 : 1;
    total += static_cast<double>(k);
  }
  CHECK(outside == 0);
  CHECK_NEAR(total / trials, 2.5, 4.0 * std::sqrt(0.25 / trials));
  CHECK(photon_count(star_source(inf, 7.0), 3, 0, 0) == 7);
  // A cone target: photons * (1 - cos) / 2 of the 4 pi emission.
  Source sphere = lamp_source(1e6);
  sphere.target = SphereTarget{{0.0, 0.0, 0.0}, 2.0};
  const double cos_max = std::sqrt(1.0 - 4.0 / 100.0);
  CHECK_NEAR(expected_count(sphere), 1e6 * 0.5 * (1.0 - cos_max), 1e-6);
}

void sphere_target() {
  // Every photon of the cone hits the sphere (perpendicular distance to the centre < R).
  Source s = lamp_source(1e5);
  s.target = SphereTarget{{0.0, 0.0, 0.0}, 2.0};
  const auto c = emit_all(s, 7, 5000);
  for (std::size_t i = 0; i < c.rows(); ++i) {
    CHECK(c[Column::w][i] == 1.0f);
    const Vec3 origin{c[Column::x][i], c[Column::y][i], c[Column::z][i]};
    const Vec3 d{c[Column::cx][i], c[Column::cy][i], c[Column::cz][i]};
    const Vec3 to_centre = Vec3{0.0, 0.0, 0.0} - origin;
    const double along = dot(to_centre, d);
    const double perpendicular = norm(to_centre - d * along);
    CHECK(perpendicular <= 2.0 + 1e-4);
  }
}

void validation() {
  const auto rejects_with = [](const Source &s, const char *fragment) {
    return throws_with([&] { validate(s); }, fragment);
  };
  Source s = lamp_source(10.0);
  CHECK(accepts([&] { validate(s); }));
  Source bad = s;
  bad.photons = -1.0;
  CHECK(rejects_with(bad, "photons must be finite"));
  bad = s;
  bad.geometry = Point{{0, 0, 1}, {0, 0, 2}, Isotropic{}};
  CHECK(rejects_with(bad, "axis must be a unit vector"));
  bad = s;
  bad.geometry = Point{{0, 0, 10}, {0, 0, -1}, FixedAngle{0.1}};
  bad.target = DiscTarget{{0, 0, 0}, {0, 0, 1}, 2.0};
  CHECK(rejects_with(bad, "needs a density")); // No density for a disc target.
  bad = s;
  bad.geometry = Point{{0, 0, 10}, {0, 0, -1}, Lambertian{}};
  bad.target = SphereTarget{{0, 0, 0}, 2.0};
  CHECK(rejects_with(bad, "needs an isotropic point")); // A sphere target needs an isotropic point.
  bad = s;
  bad.target = SphereTarget{{0, 0, 9}, 2.0};
  CHECK(rejects_with(bad, "sphere")); // The point lies inside the sphere.
  bad = star_source(inf);
  bad.target = SphereTarget{{0, 0, 9}, 2.0};
  CHECK(rejects_with(bad, "takes no target")); // A plane wave takes no target.
  bad = star_source(inf);
  bad.geometry = PlaneWave{{0.0, 0.0, -1.0}, inf, 5.0, 0.0, 30.0};
  CHECK(rejects_with(
      bad, "source side")); // Upward rays from a start plane above the disc: the wrong side.
  bad = star_source(inf);
  bad.geometry = PlaneWave{{1.0, 0.0, 0.0}, inf, 5.0, 0.0, 30.0};
  CHECK(rejects_with(bad, "start_distance_m")); // Horizontal with start_z_m: no start plane z =
                                                // const (use start_distance_m).
  bad = star_source(inf);
  bad.geometry = PlaneWave{{0.0, 0.0, 1.0}, inf, 5.0, 0.0, -30.0};
  CHECK(rejects_with(bad, "source side")); // Downward rays must start above the disc.
  bad = star_source(inf);
  bad.geometry = PlaneWave{{0.0, 0.0, -1.0}, 20.0, 5.0, 0.0, 30.0};
  CHECK(rejects_with(bad, "source side")); // Finite distance: the same side rule.
  bad = s;
  bad.spectrum = InverseSquare{600.0, 300.0};
  CHECK(rejects_with(bad, "need 0 < low"));
  bad = s;
  bad.pulse = Exponential{0.0, -1.0};
  CHECK(rejects_with(bad, "pulse:"));
}

// The start rules of plane waves: exactly one of start_z_m and start_distance_m; start_z_m only
// for non-horizontal directions (the message points to start_distance_m); start_distance_m at
// least the disc's extent along towards_source.
void plane_wave_start_rules() {
  Source s;
  s.photons = 10.0;
  s.geometry = PlaneWave{{1.0, 0.0, 0.0}, inf, 5.0, 0.0, 30.0, std::nullopt};
  CHECK(throws_with([&] { validate(s); }, "start_distance_m"));
  try {
    validate(s);
  } catch (const std::invalid_argument &e) {
    CHECK(std::string(e.what()).find("start_distance_m") != std::string::npos);
  }
  s.geometry = PlaneWave{{1.0, 0.0, 0.0}, inf, 5.0, 0.0, std::nullopt, 30.0};
  CHECK(accepts([&] { validate(s); }));
  s.geometry = PlaneWave{{0.0, 0.0, 1.0}, inf, 5.0, 0.0, std::nullopt, std::nullopt};
  CHECK(accepts([&] { validate(s); })); // neither start: start_z_m = 0 (as before)
  s.geometry = PlaneWave{{1.0, 0.0, 0.0}, inf, 5.0, 0.0, std::nullopt, std::nullopt};
  CHECK(throws_with([&] { validate(s); },
                    "start_distance_m")); // horizontal with neither: needs start_distance_m
  s.geometry = PlaneWave{{0.0, 0.0, 1.0}, inf, 5.0, 0.0, 30.0, 30.0};
  CHECK(throws_with([&] { validate(s); }, "start")); // both starts
  s.geometry = PlaneWave{{1.0, 0.0, 0.0}, inf, 5.0, 0.0, std::nullopt, 4.0};
  CHECK(throws_with([&] { validate(s); },
                    "start_distance_m must be at least")); // closer than the disc's extent (5 m)
  s.geometry = PlaneWave{{1.0, 0.0, 0.0}, inf, 5.0, 0.0, std::nullopt, 5.0};
  CHECK(accepts([&] { validate(s); }));
  s.geometry = PlaneWave{{1.0, 0.0, 0.0}, inf, 5.0, 0.0, std::nullopt, -1.0};
  CHECK(throws_with([&] { validate(s); }, "start_distance_m must be finite"));
  // Direction 60 deg from the vertical: the extent is 5 sqrt(1 - 0.25) = 4.33 m.
  s.geometry = PlaneWave{normalised({0.866, 0.0, 0.5}), inf, 5.0, 0.0, std::nullopt, 4.5};
  CHECK(accepts([&] { validate(s); }));
  s.geometry = PlaneWave{normalised({0.866, 0.0, 0.5}), inf, 5.0, 0.0, std::nullopt, 4.0};
  CHECK(throws_with([&] { validate(s); }, "start_distance_m must be at least"));
}

void upward_plane_waves() {
  // Beam from below: towards_source = (0.01, 0, -1), disc at z = 0, start plane at z = -30 m.
  // Photons travel up along (-0.01, 0, 1) normalised, start on the plane, times relative to the
  // wavefront through the disc centre.
  Source up;
  up.photons = 20000.0;
  up.geometry = PlaneWave{normalised({0.01, 0.0, -1.0}), inf, 5.0, 0.0, -30.0};
  CHECK(accepts([&] { validate(up); }));
  const auto c = emit_all(up, 7, 20000);
  CHECK(c.rows() == 20000);
  const Vec3 towards = normalised({0.01, 0.0, -1.0});
  double inside = 0.0;
  for (std::size_t i = 0; i < c.rows(); ++i) {
    CHECK_NEAR(c[Column::z][i], -30.0, 1e-4);
    CHECK(c[Column::cz][i] > 0.0);
    CHECK_NEAR(c[Column::cz][i], normalised({-0.01, 0.0, 1.0}).z, 1e-6);
    const Vec3 target = at_plane(c, i, 0.0);
    CHECK(std::hypot(target.x, target.y) <= 5.0 + 1e-3);
    inside += std::hypot(target.x, target.y) < 2.5;
    const Vec3 start{c[Column::x][i], c[Column::y][i], c[Column::z][i]};
    const Vec3 dir{c[Column::cx][i], c[Column::cy][i], c[Column::cz][i]};
    CHECK_NEAR(dot(dir, towards), -1.0, 1e-6);
    CHECK_NEAR(c[Column::t][i], dot(start, dir) / speed_of_light_m_per_ns, 1e-3);
  }
  CHECK_NEAR(inside / 20000.0, 0.25, 0.01);

  // Point source below the disc at a finite distance: source at (0, 0, 1.5) - 20 * (0.6, 0, 0.8)
  // = (-12, 0, -14.5). A start plane above the source (z = -5) is crossed by the photons: they
  // start there. A start plane below the source (z = -30): they start at the source itself.
  Source lamp;
  lamp.photons = 50.0;
  lamp.geometry = PlaneWave{normalised({-0.6, 0.0, -0.8}), 20.0, 1.0, 1.5, -5.0};
  const auto p = emit_all(lamp, 3, 50);
  CHECK(p.rows() == 50);
  for (std::size_t i = 0; i < p.rows(); ++i) {
    CHECK_NEAR(p[Column::z][i], -5.0, 1e-4);
    CHECK(p[Column::cz][i] > 0.0);
    CHECK(p[Column::x][i] > -12.0); // on the way from the source (x = -12) to the disc
  }
  lamp.geometry = PlaneWave{normalised({-0.6, 0.0, -0.8}), 20.0, 1.0, 1.5, -30.0};
  const auto q = emit_all(lamp, 3, 50);
  for (std::size_t i = 0; i < q.rows(); ++i) {
    CHECK_NEAR(q[Column::x][i], -12.0, 1e-4);
    CHECK_NEAR(q[Column::y][i], 0.0, 1e-4);
    CHECK_NEAR(q[Column::z][i], -14.5, 1e-4);
    // Start at the source, 20 m before the disc-centre wavefront: t = -d / c.
    CHECK_NEAR(q[Column::t][i], -20.0 / speed_of_light_m_per_ns, 1e-3);
  }
}

void validation_bounds() {
  // Photon counts: non-finite and beyond 2^53 (not exact as a count) are rejected.
  CHECK(throws_with([] { validate(lamp_source(std::numeric_limits<double>::quiet_NaN())); },
                    "photons must be finite"));
  CHECK(throws_with([] { validate(lamp_source(inf)); }, "photons must be finite"));
  CHECK(throws_with([] { validate(lamp_source(1e17)); }, "at most 2^53"));
  CHECK(accepts([] { validate(lamp_source(9007199254740992.0)); }));
  // A disc target at h = 1e-170 (h^2 underflows to zero): the fraction is the hemisphere, not NaN.
  Source near = lamp_source(1.0);
  near.geometry = Point{{0.0, 0.0, 1e-170}, {0.0, 0.0, -1.0}, Isotropic{}};
  near.target = DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0};
  // The disc then subtends the hemisphere: fraction 1/2 (no overflow).
  CHECK_NEAR(expected_count(near), 0.5, 1e-12);
  CHECK(photon_count(near, 1, 0, 0) <= 1);
  // A large but valid count is still accepted.
  CHECK(accepts([] { (void)expected_count(lamp_source(1e15)); }));

  // Complex pulse: every width must be finite (+inf used to pass the >= 0 test).
  Source pulse = lamp_source(10.0);
  pulse.pulse = ComplexPulse{inf, 0.5, 0.5, 0.0};
  CHECK(throws_with([&] { validate(pulse); }, "pulse:"));
  pulse.pulse = ComplexPulse{1.0, inf, 0.5, 0.0};
  CHECK(throws_with([&] { validate(pulse); }, "pulse:"));
  pulse.pulse = ComplexPulse{1.0, 0.5, inf, 0.0};
  CHECK(throws_with([&] { validate(pulse); }, "pulse:"));
  pulse.pulse = ComplexPulse{1.0, 0.5, 0.5, inf};
  CHECK(throws_with([&] { validate(pulse); }, "pulse:"));
  pulse.pulse = ComplexPulse{1.0, 0.5, 0.5, 0.0};
  CHECK(accepts([&] { validate(pulse); }));

  // Magnitudes: lengths, times and wavelengths stay within 1e30 (squares and float32
  // columns finite); positive wavelengths are at least 1e-30 (inverses finite).
  Source big = lamp_source(10.0);
  big.geometry = Point{{1e31, 0.0, 0.0}, {0.0, 0.0, -1.0}, Isotropic{}};
  CHECK(throws_with([&] { validate(big); }, "out of range"));
  big = lamp_source(10.0);
  big.pulse = Fixed{1e31};
  CHECK(throws_with([&] { validate(big); }, "time out of range"));
  big = lamp_source(10.0);
  big.spectrum = Fixed{1e-31};
  CHECK(throws_with([&] { validate(big); }, "wavelength out of range"));
  big = lamp_source(10.0);
  big.spectrum = InverseSquare{1e-320, 600.0}; // 1 / low would overflow.
  CHECK(throws_with([&] { validate(big); }, "need 0 < low < high"));
  big = star_source(1e200);
  CHECK(throws_with([&] { validate(big); }, "distance must be infinite or positive"));
  big = star_source(inf); // A plane wave is not a length: accepted.
  CHECK(accepts([&] { validate(big); }));
}

} // namespace

int main() {
  return run_checks([] {
    stars();
    upward_plane_waves();
    plane_wave_start_rules();
    emitters_isotropic_disc();
    disc_fraction_closed_forms();
    disc_fraction_quadrature();
    disc_target_angular_profile();
    disc_target_narrow_off_axis_beam();
    disc_target_small_distant_disc();
    chunks_and_determinism();
    directions_and_angles();
    counts();
    sphere_target();
    validation();
    validation_bounds();
  });
}
