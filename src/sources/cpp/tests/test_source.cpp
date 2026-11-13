#include "check.hpp"

#include "dlamb/source.hpp"
#include "dlamb/streams.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace dlamb;
using gstrickd::Column;

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();

// All photons of a test source, one array per column. A PhotonChunk holds at most 2^16
// rows, so the test keeps the columns itself.
struct Photons {
  gstrickd::Frame frame{gstrickd::Frame::ground};
  std::array<std::vector<float>, gstrickd::n_columns> columns{};

  std::size_t rows() const { return columns[static_cast<std::size_t>(Column::x)].size(); }
  const std::vector<float> &operator[](Column c) const {
    return columns[static_cast<std::size_t>(c)];
  }
};

// Emit photons [0, n) of a source in chunks of `chunk` rows.
Photons emit_all(const Source &source, std::uint64_t seed, std::uint64_t n,
                 std::size_t chunk = 1 << 16, std::uint32_t event = 1,
                 std::uint32_t telescope = 2) {
  Photons all;
  all.frame = source.frame;
  for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
    if (source_columns & (std::uint32_t{1} << c))
      all.columns[c].resize(n);
  }
  gstrickd::PhotonChunk part;
  for (std::uint64_t first = 0; first < n; first += chunk) {
    const std::size_t rows = static_cast<std::size_t>(std::min<std::uint64_t>(chunk, n - first));
    emit(source, seed, event, telescope, first, rows, part);
    for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
      const auto bit = std::uint32_t{1} << c;
      if (!(source_columns & bit))
        continue;
      auto &dst = all.columns[c];
      const auto &src = part.columns[c];
      for (std::size_t i = 0; i < rows; ++i)
        dst[first + i] = src[i];
    }
  }
  return all;
}

bool same_chunk(const Photons &a, const Photons &b) {
  if (a.rows() != b.rows())
    return false;
  for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c)
    if (a.columns[c] != b.columns[c])
      return false;
  return true;
}

Source star_source(double distance, double photons = 1000.0) {
  Source s;
  s.frame = gstrickd::Frame::telescope;
  s.photons = photons;
  s.geometry = Star{normalised({0.01, 0.0, 1.0}), distance, 5.0, 0.0, 30.0};
  return s;
}

Source lamp_source(double photons = 1e5) {
  Source s;
  s.frame = gstrickd::Frame::camera;
  s.photons = photons;
  s.geometry = Emitter{{0.0, 0.0, 10.0}, {0.0, 0.0, -1.0}, Isotropic{}};
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
  // Plane wave over a 5 m disc (obdeect test_kernels.cpp sources): all rays start at the
  // start plane z = 30 m, a quarter of the photons land within 2.5 m of the centre.
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
  // photons leave the source itself (obdeect PointSource, finite distance).
  Source lamp;
  lamp.frame = gstrickd::Frame::telescope;
  lamp.photons = 50.0;
  lamp.geometry = Star{normalised({0.6, 0.0, 0.8}), 20.0, 1.0, 1.5, 30.0};
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

void emitters_isotropic_disc() {
  // Isotropic emitter 10 m above a 2 m disc (obdeect test_kernels.cpp emitters): the
  // candidates per emitted photon are R^2 / (4 h^2) = 1/100.
  Source s = lamp_source(4e5);
  s.target = DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0};
  CHECK_NEAR(expected_count(s), 4e5 / 100.0, 1e-9);
  const auto c = emit_all(s, 5, 4000);
  double sum = 0.0, sum2 = 0.0;
  const std::size_t rows = c.rows();
  for (std::size_t i = 0; i < rows; ++i) {
    const double w = c[Column::w][i];
    CHECK(w > 0.0 && w <= 1.0);
    const Vec3 hit = at_plane(c, i, 0.0);
    CHECK(std::hypot(hit.x, hit.y) <= 2.0 + 1e-3);
    sum += w;
    sum2 += w * w;
  }
  // Mean weight per candidate = (fraction of emitted photons on the disc) / (candidates per
  // photon) = (solid angle / 4 pi) * 100. Checked at four sigma.
  const double mean = sum / rows, sd = std::sqrt((sum2 / rows - mean * mean) / rows);
  const double solid_angle = 2.0 * pi * (1.0 - 10.0 / std::sqrt(104.0));
  const double expected = solid_angle / (4.0 * pi) * 100.0;
  CHECK(std::abs(mean - expected) <= 4.0 * sd);
  CHECK(same_chunk(c, emit_all(s, 5, 4000, 333)));
}

void disc_target_angular_profile() {
  // A 5 deg filled cone lights only the inner part of the disc (obdeect test: weight is
  // positive inside the cone and zero outside it).
  Source s = lamp_source(1e5);
  s.geometry =
      Emitter{{0.0, 0.0, 10.0}, {0.0, 0.0, -1.0}, FilledCone{1.0 - std::cos(5.0 * pi / 180.0)}};
  s.target = DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0};
  const auto c = emit_all(s, 5, 10000);
  for (std::size_t i = 0; i < c.rows(); ++i) {
    const double angle = std::acos(-c[Column::cz][i]) * 180.0 / pi;
    // Away from the cone edge (5 deg, margin 1e-3 deg): weight > 0 inside, 0 outside.
    if (angle < 5.0 - 1e-3)
      CHECK(c[Column::w][i] > 0.0);
    else if (angle > 5.0 + 1e-3)
      CHECK(c[Column::w][i] == 0.0);
  }
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
  // A chunk starting at photon 2^20 is the same photons whatever came before.
  gstrickd::PhotonChunk late, check;
  emit(s, 9, 1, 2, 1U << 20, 100, late);
  emit(s, 9, 1, 2, 1U << 20, 50, check);
  CHECK(late.first_index == (1U << 20) && check.first_index == (1U << 20));
  for (std::size_t i = 0; i < 50; ++i)
    CHECK(late[Column::x][i] == check[Column::x][i] && late[Column::t][i] == check[Column::t][i]);
  bool rejected = false;
  try {
    gstrickd::PhotonChunk big;
    emit(s, 9, 1, 2, 0, (1U << 16) + 1, big);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
}

void directions_and_frames() {
  // Directions are unit vectors for every profile; positions and frames are the source's.
  const std::vector<Angular> profiles{Isotropic{}, FixedAngle{1.0 - std::cos(0.2)}, FilledCone{0.3},
                                      Lambertian{}, TabulatedAngular{gaussian_angular_table(0.2)}};
  for (const auto &profile : profiles) {
    Source s = lamp_source(1e4);
    s.frame = gstrickd::Frame::camera;
    s.geometry = Emitter{{1.0, -2.0, 3.0}, normalised({0.2, 0.1, -1.0}), profile};
    const auto c = emit_all(s, 1, 2000, 512);
    CHECK(c.frame == gstrickd::Frame::camera || c.rows() == 0);
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
  cone.geometry = Emitter{{0.0, 0.0, 10.0}, axis, FixedAngle{1.0 - std::cos(0.2)}};
  const auto c = emit_all(cone, 4, 500);
  for (std::size_t i = 0; i < c.rows(); ++i) {
    const double cosine =
        c[Column::cx][i] * axis.x + c[Column::cy][i] * axis.y + c[Column::cz][i] * axis.z;
    CHECK_NEAR(cosine, std::cos(0.2), 1e-6);
  }
  // A star in frame T keeps its start-plane height in the telescope frame.
  const auto s = emit_all(star_source(inf, 10.0), 1, 10);
  CHECK(s.frame == gstrickd::Frame::telescope);
  CHECK_NEAR(s[Column::z][0], 30.0, 1e-4);
}

void counts() {
  // Counts: floor of the mean plus one with the fractional probability.
  const Source s = star_source(inf, 2.5);
  double total = 0.0;
  const int trials = 40000;
  for (int seed = 0; seed < trials; ++seed) {
    const auto k = photon_count(s, static_cast<std::uint64_t>(seed), 0, 0);
    CHECK(k == 2 || k == 3);
    total += static_cast<double>(k);
  }
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
  const auto rejects = [](const Source &s) {
    try {
      validate(s);
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  Source s = lamp_source(10.0);
  CHECK(!rejects(s));
  Source bad = s;
  bad.photons = -1.0;
  CHECK(rejects(bad));
  bad = s;
  bad.frame = gstrickd::Frame::ground;
  CHECK(rejects(bad));
  bad = s;
  bad.geometry = Emitter{{0, 0, 1}, {0, 0, 2}, Isotropic{}};
  CHECK(rejects(bad));
  bad = s;
  bad.geometry = Emitter{{0, 0, 10}, {0, 0, -1}, FixedAngle{0.1}};
  bad.target = DiscTarget{{0, 0, 0}, {0, 0, 1}, 2.0};
  CHECK(rejects(bad)); // No density for a disc target.
  bad = s;
  bad.geometry = Emitter{{0, 0, 10}, {0, 0, -1}, Lambertian{}};
  bad.target = SphereTarget{{0, 0, 0}, 2.0};
  CHECK(rejects(bad)); // A sphere target needs an isotropic emitter.
  bad = s;
  bad.target = SphereTarget{{0, 0, 9}, 2.0};
  CHECK(rejects(bad)); // The emitter lies inside the sphere.
  bad = star_source(inf);
  bad.target = SphereTarget{{0, 0, 9}, 2.0};
  CHECK(rejects(bad)); // A star takes no target.
  bad = star_source(inf);
  bad.geometry = Star{{0.0, 0.0, -1.0}, inf, 5.0, 0.0, 30.0};
  CHECK(rejects(bad));
  bad = s;
  bad.spectrum = InverseSquare{600.0, 300.0};
  CHECK(rejects(bad));
  bad = s;
  bad.pulse = Exponential{0.0, -1.0};
  CHECK(rejects(bad));
}

template <class F> bool throws_invalid(F f) {
  try {
    f();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

void validation_bounds() {
  // Photon counts: non-finite and beyond 2^53 (not exact as a count) are rejected.
  CHECK(throws_invalid([] { validate(lamp_source(std::numeric_limits<double>::quiet_NaN())); }));
  CHECK(throws_invalid([] { validate(lamp_source(inf)); }));
  CHECK(throws_invalid([] { validate(lamp_source(1e17)); }));
  CHECK(!throws_invalid([] { validate(lamp_source(9007199254740992.0)); }));
  // A disc target at h = 1e-170 overflows the candidate count (h^2 underflows to zero):
  // expected_count and photon_count throw instead of returning inf.
  Source near = lamp_source(1.0);
  near.geometry = Emitter{{0.0, 0.0, 1e-170}, {0.0, 0.0, -1.0}, Isotropic{}};
  near.target = DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0};
  CHECK(throws_invalid([&] { (void)expected_count(near); }));
  CHECK(throws_invalid([&] { (void)photon_count(near, 1, 0, 0); }));
  // A large but valid count is still accepted.
  CHECK(!throws_invalid([] { (void)expected_count(lamp_source(1e15)); }));

  // Frames: only sky, T and C are source frames (G and unknown codes are rejected).
  Source frame = lamp_source(10.0);
  frame.frame = static_cast<gstrickd::Frame>(7);
  CHECK(throws_invalid([&] { validate(frame); }));
  frame.frame = gstrickd::Frame::sky;
  CHECK(!throws_invalid([&] { validate(frame); }));
  frame.frame = gstrickd::Frame::telescope;
  CHECK(!throws_invalid([&] { validate(frame); }));

  // Complex pulse: every width must be finite (+inf used to pass the >= 0 test).
  Source pulse = lamp_source(10.0);
  pulse.pulse = ComplexPulse{inf, 0.5, 0.5, 0.0};
  CHECK(throws_invalid([&] { validate(pulse); }));
  pulse.pulse = ComplexPulse{1.0, inf, 0.5, 0.0};
  CHECK(throws_invalid([&] { validate(pulse); }));
  pulse.pulse = ComplexPulse{1.0, 0.5, inf, 0.0};
  CHECK(throws_invalid([&] { validate(pulse); }));
  pulse.pulse = ComplexPulse{1.0, 0.5, 0.5, inf};
  CHECK(throws_invalid([&] { validate(pulse); }));
  pulse.pulse = ComplexPulse{1.0, 0.5, 0.5, 0.0};
  CHECK(!throws_invalid([&] { validate(pulse); }));

  // Magnitudes: lengths, times and wavelengths stay within 1e30 (squares and float32
  // columns finite); positive wavelengths are at least 1e-30 (inverses finite).
  Source big = lamp_source(10.0);
  big.geometry = Emitter{{1e31, 0.0, 0.0}, {0.0, 0.0, -1.0}, Isotropic{}};
  CHECK(throws_invalid([&] { validate(big); }));
  big = lamp_source(10.0);
  big.pulse = Fixed{1e31};
  CHECK(throws_invalid([&] { validate(big); }));
  big = lamp_source(10.0);
  big.spectrum = Fixed{1e-31};
  CHECK(throws_invalid([&] { validate(big); }));
  big = lamp_source(10.0);
  big.spectrum = InverseSquare{1e-320, 600.0}; // 1 / low would overflow.
  CHECK(throws_invalid([&] { validate(big); }));
  big = star_source(1e200);
  CHECK(throws_invalid([&] { validate(big); }));
  big = star_source(inf); // A plane wave is not a length: accepted.
  CHECK(!throws_invalid([&] { validate(big); }));
}

} // namespace

int main() {
  return run_checks([] {
    stars();
    emitters_isotropic_disc();
    disc_target_angular_profile();
    chunks_and_determinism();
    directions_and_frames();
    counts();
    sphere_target();
    validation();
    validation_bounds();
  });
}
