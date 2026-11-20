#include "check.hpp"

#include "obdeect/sources/emit.hpp"
#include "obdeect/sources/source.hpp"

#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace obdeect::sources;
using obdeect::stream::Column;

namespace {

// All photons of a test source, one array per column (rows dropped by emit are not stored).
struct Photons {
  std::array<std::vector<float>, obdeect::stream::n_columns> columns{};

  std::size_t rows() const { return columns[static_cast<std::size_t>(Column::x)].size(); }
  const std::vector<float> &operator[](Column c) const {
    return columns[static_cast<std::size_t>(c)];
  }
};

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

template <class F> bool throws_invalid(F f) {
  try {
    f();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

Vec3 position_of(const Photons &c, std::size_t i) {
  return {c[Column::x][i], c[Column::y][i], c[Column::z][i]};
}

Vec3 direction_of(const Photons &c, std::size_t i) {
  return {c[Column::cx][i], c[Column::cy][i], c[Column::cz][i]};
}

// Distance along the line from its start to photon i of c (the projection on the direction).
double along_line(const Line &l, const Photons &c, std::size_t i) {
  return dot(position_of(c, i) - l.start, l.direction);
}

// A line with a fixed or isotropic angular profile and a sphere (or no) target.
Source line_source(const Line &line, double photons, const Target &target = std::monostate{}) {
  Source s;
  s.photons = photons;
  s.geometry = line;
  s.target = target;
  return s;
}

// Line that passes 30 m from a sphere of radius 15 m at the origin, from z = 300 m down to 0.
Line sphere_line(const Angular &angular) {
  Line l;
  l.start = {30.0, 0.0, 300.0};
  l.direction = {0.0, 0.0, -1.0};
  l.length_m = 300.0;
  l.angular = angular;
  return l;
}

const SphereTarget sphere_15{{0.0, 0.0, 0.0}, 15.0};

// Distance from the sphere centre to the ray of photon i (infinite line, through the photon).
double ray_distance_to_centre(const Photons &c, std::size_t i, const SphereTarget &sphere) {
  const Vec3 to_centre = sphere.centre - position_of(c, i);
  const Vec3 d = direction_of(c, i);
  return norm(to_centre - d * dot(to_centre, d));
}

void line_fixed_angle_no_target() {
  // Cherenkov-like line with no target: every candidate is emitted at the cone angle, on the
  // segment, at time = pulse + s / (beta c); the spectrum follows 1/lambda^2.
  const double theta = 2.0 * pi / 180.0;
  Line l;
  l.start = {10.0, -5.0, 200.0};
  l.direction = normalised({0.2, -0.1, -1.0});
  l.length_m = 100.0;
  l.beta = 0.5;
  l.angular = FixedAngle{1.0 - std::cos(theta)};
  Source s = line_source(l, 2000.0);
  s.pulse = Fixed{5.0};
  s.spectrum = InverseSquare{300.0, 600.0};
  const std::size_t n = 2000;
  const auto c = emit_all(s, 21, n);
  CHECK(c.rows() == n); // nothing is dropped without a target
  double sum_s = 0.0, s1 = 0.0, s2 = 0.0;
  for (std::size_t i = 0; i < c.rows(); ++i) {
    CHECK_NEAR(dot(direction_of(c, i), l.direction), std::cos(theta), 1e-6);
    const double along = along_line(l, c, i);
    CHECK(along >= -1e-3 && along <= l.length_m + 1e-3);
    CHECK(norm(position_of(c, i) - l.start - l.direction * along) < 1e-3);
    CHECK_NEAR(c[Column::t][i], 5.0 + along / (l.beta * speed_of_light_m_per_ns), 1e-3);
    CHECK(c[Column::w][i] == 1.0f);
    sum_s += along;
    const double wl = c[Column::wavelength][i];
    CHECK(wl >= 300.0 && wl <= 600.0);
    s1 += wl;
    s2 += wl * wl;
  }
  // Positions are uniform along the segment: mean L / 2 within 4 sigma.
  const double mean_s = sum_s / n, sd_s = std::sqrt(l.length_m * l.length_m / 12.0);
  CHECK(std::abs(mean_s - 0.5 * l.length_m) <= 4.0 * sd_s / std::sqrt(static_cast<double>(n)));
  // Wavelength moments of the 1/lambda^2 law within 4 sigma.
  const InverseSquare law{300.0, 600.0};
  const double mean = s1 / n, var = s2 / n - mean * mean;
  CHECK(std::abs(mean - law.mean()) <= 4.0 * std::sqrt(law.variance() / n));
  CHECK(std::abs(var - law.variance()) <= 4.0 * law.variance() * std::sqrt(2.0 / n));
}

void line_sphere_fixed_angle() {
  // Fixed-angle cone from a line that passes the sphere at the side: only candidates whose cone
  // hits the sphere are kept, each within R of the centre, with weight = arc fraction (half / pi)
  // over the reference (largest arc fraction of the scan, so weights are at most ~1).
  const Line l = sphere_line(FixedAngle{1.0 - std::cos(5.0 * pi / 180.0)});
  const double reference = line_reference_fraction(l, sphere_15);
  CHECK(reference > 0.0 && reference <= 1.0);
  const std::size_t n = 20000;
  // Candidates = expected_count = photons * reference.
  Source s = line_source(l, static_cast<double>(n) / reference, sphere_15);
  const auto c = emit_all(s, 22, n);
  CHECK(c.rows() <= n);
  CHECK(c.rows() > 0 && c.rows() < n); // some candidates miss the sphere
  double sum = 0.0, sum2 = 0.0;
  for (std::size_t i = 0; i < c.rows(); ++i) {
    CHECK(ray_distance_to_centre(c, i, sphere_15) <= sphere_15.radius_m + 1e-3);
    const double w = c[Column::w][i];
    CHECK(w > 0.0 && w <= 1.01); // the scan maximum may miss the true maximum slightly
    sum += w;
    sum2 += w * w;
  }
  const double expected = expected_weight(s);
  CHECK(expected > 0.0);
  CHECK(std::abs(sum - expected) <= 4.0 * std::sqrt(sum2));
}

void line_sphere_isotropic() {
  // Isotropic line with a sphere target: every candidate is kept, its cone covers the sphere,
  // and the weight is its sphere fraction over the largest one along the line.
  const Line l = sphere_line(Isotropic{});
  const double fraction_max = line_sphere_fraction_max(l, sphere_15);
  const std::size_t n = 20000;
  // Candidates = expected_count = photons * fraction_max.
  Source s = line_source(l, static_cast<double>(n) / fraction_max, sphere_15);
  const auto c = emit_all(s, 23, n);
  CHECK(c.rows() == n);
  double sum = 0.0, sum2 = 0.0;
  for (std::size_t i = 0; i < c.rows(); ++i) {
    CHECK(ray_distance_to_centre(c, i, sphere_15) <= sphere_15.radius_m + 1e-3);
    const double w = c[Column::w][i];
    CHECK(w > 0.0 && w <= 1.0);
    sum += w;
    sum2 += w * w;
  }
  const double expected = expected_weight(s);
  CHECK(std::abs(sum - expected) <= 4.0 * std::sqrt(sum2));
}

void line_profile() {
  // Triangular density on [0, 10] m with its peak at 2 m: mean (0 + 2 + 10) / 3 = 4 m and
  // P(s < 2 m) = 1 / 5.
  Line l;
  l.start = {0.0, 0.0, 100.0};
  l.direction = {0.0, 0.0, -1.0};
  l.length_m = 10.0;
  l.profile = LinearDensity({0.0, 2.0, 10.0}, {0.0, 1.0, 0.0});
  const std::size_t n = 20000;
  const auto c = emit_all(line_source(l, double(n)), 24, n);
  CHECK(c.rows() == n);
  double sum = 0.0;
  std::size_t below = 0;
  for (std::size_t i = 0; i < c.rows(); ++i) {
    const double along = along_line(l, c, i);
    CHECK(along >= -1e-3 && along <= 10.0 + 1e-3);
    sum += along;
    below += along < 2.0;
  }
  const double var = (0.0 + 4.0 + 100.0 - 20.0) / 18.0; // triangular variance
  CHECK(std::abs(sum / n - 4.0) <= 4.0 * std::sqrt(var / n));
  CHECK(std::abs(static_cast<double>(below) / n - 0.2) <= 4.0 * std::sqrt(0.16 / n));
}

void rows_independent_of_chunk() {
  // Dropped rows included: the same rows for chunks of 1000 and 65536 candidates.
  const std::size_t n = 5000;
  const Source fixed =
      line_source(sphere_line(FixedAngle{1.0 - std::cos(5.0 * pi / 180.0)}), double(n), sphere_15);
  const auto whole = emit_all(fixed, 25, n, 1 << 16);
  CHECK(whole.rows() > 0 && whole.rows() < n);
  CHECK(same_chunk(whole, emit_all(fixed, 25, n, 1000)));

  const Source iso = line_source(
      sphere_line(Isotropic{}),
      double(n) / line_sphere_fraction_max(sphere_line(Isotropic{}), sphere_15), sphere_15);
  const auto iso_whole = emit_all(iso, 25, n, 1 << 16);
  CHECK(iso_whole.rows() == n);
  CHECK(same_chunk(iso_whole, emit_all(iso, 25, n, 1000)));
}

void sky_background() {
  // Photons start on the sphere and travel inwards, within the view cone about the axis.
  SkyBackground b;
  b.centre = {0.0, 0.0, 100.0};
  b.radius_m = 500.0;
  b.axis = normalised({0.0, 0.2, -1.0});
  b.viewcone_deg = 10.0;
  Source s;
  s.photons = 5000.0;
  s.geometry = b;
  const auto c = emit_all(s, 26, 5000);
  CHECK(c.rows() == 5000);
  const double cos_cone = std::cos(b.viewcone_deg * pi / 180.0);
  for (std::size_t i = 0; i < c.rows(); ++i) {
    const Vec3 rel = position_of(c, i) - b.centre;
    CHECK_NEAR(norm(rel), b.radius_m, 1e-3);
    const Vec3 d = direction_of(c, i);
    CHECK_NEAR(norm(d), 1.0, 1e-6);
    CHECK(dot(d, b.axis) >= cos_cone - 1e-6);
    CHECK(dot(d, rel) / b.radius_m <= 1e-5); // travelling inwards
  }
}

void arc_matches_brute_force() {
  // Arc of a fixed cone from a point towards a sphere against brute force over the azimuth.
  const Vec3 start{20.0, -5.0, 1000.0};
  const Vec3 direction = normalised({-0.2, 0.1, -1.0});
  const double cos_theta = std::cos(19.0 * pi / 180.0);
  const SphereTarget sphere{{0.0, 0.0, 0.0}, 15.0};
  const double s = 400.0;
  const Vec3 point = start + direction * s;
  const Arc arc = detail::arc_towards(point, direction, cos_theta, sphere);
  const Basis b = orthonormal_basis(direction);
  const double sin_t = std::sin(std::acos(cos_theta));
  const std::size_t n = 200000;
  std::size_t hits = 0;
  const Vec3 to = normalised(sphere.centre - point);
  const double cos_alpha = std::sqrt(1.0 - std::pow(15.0 / norm(sphere.centre - point), 2));
  for (std::size_t k = 0; k < n; ++k) {
    const double phi = two_pi * (static_cast<double>(k) + 0.5) / n;
    const Vec3 e = direction * cos_theta + (b.e1 * std::cos(phi) + b.e2 * std::sin(phi)) * sin_t;
    hits += dot(e, to) > cos_alpha;
  }
  CHECK(arc.hit);
  CHECK_NEAR(arc.half / pi, static_cast<double>(hits) / n, 2e-5);
}

void arc_on_axis() {
  // Point on the axis above the sphere: all azimuths hit if the cone is inside asin(R / d).
  const Vec3 start{0.0, 0.0, 3000.0};
  const Vec3 down{0.0, 0.0, -1.0};
  const SphereTarget sphere{{0.0, 0.0, 0.0}, 15.0};
  const auto inside = detail::arc_towards(start, down, std::cos(0.2 * pi / 180.0), sphere);
  CHECK(inside.hit);
  CHECK_NEAR(inside.half, pi, 1e-12);
  CHECK(!detail::arc_towards(start, down, std::cos(1.0 * pi / 180.0), sphere).hit);
  // A point inside the sphere is hit in every direction.
  const auto within =
      detail::arc_towards({0.0, 0.0, 5.0}, down, std::cos(60.0 * pi / 180.0), sphere);
  CHECK(within.hit);
  CHECK_NEAR(within.half, pi, 1e-12);
}

void line_validation() {
  const auto base = [] {
    Line l;
    l.start = {0.0, 0.0, 100.0};
    l.length_m = 10.0;
    return line_source(l, 10.0);
  };
  CHECK(!throws_invalid([&] { validate(base()); }));
  Source bad = base();
  std::get<Line>(bad.geometry).length_m = 1e100;
  CHECK(throws_invalid([&] { validate(bad); }));
  bad = base();
  std::get<Line>(bad.geometry).start.x = 1e100;
  CHECK(throws_invalid([&] { validate(bad); }));
  bad = base();
  std::get<Line>(bad.geometry).direction = {0.0, 0.0, -2.0};
  CHECK(throws_invalid([&] { validate(bad); }));
  bad = base();
  std::get<Line>(bad.geometry).beta = 1.5;
  CHECK(throws_invalid([&] { validate(bad); }));
  bad = base();
  std::get<Line>(bad.geometry).profile = LinearDensity({0.0, 20.0}, {1.0, 1.0});
  CHECK(throws_invalid([&] { validate(bad); })); // profile beyond the length
  bad = base();
  bad.target = DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0};
  CHECK(throws_invalid([&] { validate(bad); })); // a line takes no disc target
  bad = base();
  std::get<Line>(bad.geometry).angular = Lambertian{};
  bad.target = sphere_15;
  CHECK(throws_invalid([&] { validate(bad); })); // sphere target needs fixed or isotropic
  // Photon counts above 2^53 are rejected, as for the other geometries.
  bad = base();
  bad.photons = 1e16;
  CHECK(throws_invalid([&] { validate(bad); }));
  CHECK(throws_invalid([&] { (void)emit_all(bad, 1, 10); }));
  bad.photons = 9007199254740992.0;
  CHECK(!throws_invalid([&] { validate(bad); }));
}

} // namespace

int main() {
  return run_checks([] {
    line_fixed_angle_no_target();
    line_sphere_fixed_angle();
    line_sphere_isotropic();
    line_profile();
    rows_independent_of_chunk();
    sky_background();
    arc_matches_brute_force();
    arc_on_axis();
    line_validation();
  });
}
