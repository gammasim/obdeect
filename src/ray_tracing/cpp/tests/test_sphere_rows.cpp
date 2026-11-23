// Sphere rows (G rows of a sphere run) traced by trace_bunches: each row starts at its own point.
#include "check.hpp"

#include "obdeect/ray_tracing/bunches.hpp"
#include "obdeect/ray_tracing/trace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <vector>

using namespace obdeect::ray_tracing;

namespace {

constexpr double kLightSpeed = 0.299792458; // m/ns
constexpr double kSphereRadius = 25.0; // Encloses the test optics (dish at z = 0, focus at 10).
constexpr double kPlaneHeight = 20.0;  // Input plane of plane rows (above the dish, inside R).

// One large spherical facet of focal length f; focal plane at z = f, facet vertex at the origin.
OpticalModel single_facet(double f, double half_width) {
  OpticalModel m;
  m.primary.facets.push_back({1,
                              {0, 0, 0},
                              {0, 0, 1},
                              {1, 0, 0},
                              {0, 1, 0},
                              {1, 0, 0},
                              {0, 1, 0},
                              half_width,
                              2 * f,
                              Aperture::circle});
  m.primary.grid = FacetGrid(m.primary.facets);
  m.focal_surface = {2, f, std::nullopt};
  return m;
}

// Bunch columns of n photons (one photon per bunch), owned by the caller.
struct Columns {
  std::vector<float> x, y, cx, cy, time, zem, size, wavelength, z, cz;
  std::vector<std::uint64_t> id;

  [[nodiscard]] BunchBlock block(bool with_z) const {
    BunchBlock b{x.data(),    y.data(),   cx.data(),   cy.data(),
                 time.data(), zem.data(), size.data(), wavelength.data(),
                 id.data(),   x.size(),   cz.data(),   with_z ? z.data() : nullptr};
    return b;
  }
};

// Result columns of a trace_bunches call, as vectors.
struct Output {
  std::vector<double> position, direction, time, weight, wavelength;
  std::vector<std::uint64_t> index;
  std::vector<std::uint32_t> sub;
  std::size_t rows() const { return time.size(); }
  // Keeps the first n rows (trace_bunches returns the rows written; the allocation may be larger).
  void resize_to(std::size_t n) {
    position.resize(3 * n);
    direction.resize(3 * n);
    time.resize(n);
    weight.resize(n);
    wavelength.resize(n);
    index.resize(n);
    sub.resize(n);
  }
};

Output run(const OpticalModel &model, const Columns &in, const BunchFrame &frame, bool with_z,
           unsigned threads, const RandomKey &seed = {20261123, 1, 1}) {
  Output out;
  ThreadPool pool(threads);
  BunchOptions options;
  options.roulette = false;
  const auto allocate = [&](std::size_t rows) {
    out.position.resize(3 * rows);
    out.direction.resize(3 * rows);
    out.time.resize(rows);
    out.weight.resize(rows);
    out.wavelength.resize(rows);
    out.index.resize(rows);
    out.sub.resize(rows);
    return BunchResult{out.position.data(), out.direction.data(),  out.time.data(),
                       out.weight.data(),   out.wavelength.data(), out.index.data(),
                       out.sub.data()};
  };
  const std::size_t written =
      trace_bunches(model, in.block(with_z), frame, options, seed, pool, allocate);
  out.resize_to(written);
  return out;
}

// Downward photons of a beam: the plane start (on the input plane at height kPlaneHeight, time
// there) and the sphere row (the entry point of the same line into the sphere of radius
// kSphereRadius centred on the origin, with its time). Both come from the same ground point
// (x, y, 0) and direction d, as the stream tracer computes them.
struct Pair {
  Columns plane, sphere;
};

void push_common(Columns &c, double x, double y, double z, double dx, double dy, double dz,
                 double time, std::uint64_t id) {
  c.x.push_back(static_cast<float>(x));
  c.y.push_back(static_cast<float>(y));
  c.z.push_back(static_cast<float>(z));
  c.cx.push_back(static_cast<float>(dx));
  c.cy.push_back(static_cast<float>(dy));
  c.cz.push_back(static_cast<float>(dz));
  c.time.push_back(static_cast<float>(time));
  c.zem.push_back(0.0F);
  c.size.push_back(1.0F);
  c.wavelength.push_back(400.0F);
  c.id.push_back(id);
}

Pair beam_pair(std::size_t n, double cx, double cy, double radius_m, std::uint64_t seed) {
  Pair pair;
  std::uint64_t state = seed;
  const auto uniform = [&state] {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state >> 11) * (1.0 / 9007199254740992.0);
  };
  for (std::size_t i = 0; i < n; ++i) {
    const double r = radius_m * std::sqrt(uniform()), phi = 2.0 * 3.141592653589793 * uniform();
    // Ground point of the row, as float32 (the stream stores float32 rows).
    const double ox = static_cast<double>(static_cast<float>(r * std::cos(phi)));
    const double oy = static_cast<double>(static_cast<float>(r * std::sin(phi)));
    const double dx = static_cast<double>(static_cast<float>(cx));
    const double dy = static_cast<double>(static_cast<float>(cy));
    const double dz = -std::sqrt(1.0 - dx * dx - dy * dy);
    // Plane start: P = o + s_P d on the input plane n.p = kPlaneHeight (s_P < 0 for dz < 0).
    const double s_plane = kPlaneHeight / dz;
    push_common(pair.plane, ox + s_plane * dx, oy + s_plane * dy, kPlaneHeight, dx, dy, dz,
                s_plane / kLightSpeed, i);
    // Sphere row: the entry point of the line into the sphere (the smaller root of |o + s d| = R).
    const double b = ox * dx + oy * dy;
    const double c = ox * ox + oy * oy - kSphereRadius * kSphereRadius;
    const double s = -b - std::sqrt(b * b - c);
    push_common(pair.sphere, ox + s * dx, oy + s * dy, s * dz, dx, dy, dz, s / kLightSpeed, i);
  }
  return pair;
}

// Max |a - b| over matching rows of two outputs, and the difference in bunch rows.
void check_same_detection(const Output &a, const Output &b, double tolerance_m) {
  CHECK(a.rows() == b.rows());
  CHECK(a.rows() > 0);
  double position = 0.0, direction = 0.0, time = 0.0, weight = 0.0;
  for (std::size_t i = 0; i < std::min(a.rows(), b.rows()); ++i) {
    CHECK(a.index[i] == b.index[i]);
    for (int k = 0; k < 3; ++k) {
      position = std::max(position, std::abs(a.position[3 * i + k] - b.position[3 * i + k]));
      direction = std::max(direction, std::abs(a.direction[3 * i + k] - b.direction[3 * i + k]));
    }
    time = std::max(time, std::abs(a.time[i] - b.time[i]));
    weight = std::max(weight, std::abs(a.weight[i] - b.weight[i]));
  }
  CHECK(position < tolerance_m);
  CHECK(direction < 1e-6);
  CHECK(time < 1e-3);
  CHECK(weight == 0.0);
}

BunchFrame identity_frame(double plane) {
  BunchFrame frame;
  frame.plane_z = plane;
  return frame;
}

// Plane rows and their sphere entry points give the same detected photons (downward).
void sphere_rows_match_plane_rows() {
  const auto model = single_facet(10.0, 0.6);
  const Pair pair = beam_pair(20000, 0.01, -0.005, 0.5, 7);
  const Output plane = run(model, pair.plane, identity_frame(kPlaneHeight), false, 1);
  const Output sphere = run(model, pair.sphere, identity_frame(0.0), true, 1);
  CHECK(sphere.rows() > 5000);
  check_same_detection(plane, sphere, 1e-4);
}

// The same rows give the same output for any thread count, and for the reversed row order.
void sphere_rows_invariance() {
  const auto model = single_facet(10.0, 0.6);
  const Pair pair = beam_pair(6000, -0.02, 0.01, 0.5, 11);
  const Output one = run(model, pair.sphere, identity_frame(0.0), true, 1);
  const Output four = run(model, pair.sphere, identity_frame(0.0), true, 4);
  CHECK(one.rows() == four.rows());
  CHECK(std::memcmp(one.position.data(), four.position.data(),
                    one.position.size() * sizeof(double)) == 0);
  CHECK(std::memcmp(one.time.data(), four.time.data(), one.time.size() * sizeof(double)) == 0);

  // Reversed order: the detected photons are the same set, in reverse input order.
  Columns reversed;
  const auto n = pair.sphere.x.size();
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t i = n - 1 - k;
    reversed.x.push_back(pair.sphere.x[i]);
    reversed.y.push_back(pair.sphere.y[i]);
    reversed.z.push_back(pair.sphere.z[i]);
    reversed.cx.push_back(pair.sphere.cx[i]);
    reversed.cy.push_back(pair.sphere.cy[i]);
    reversed.cz.push_back(pair.sphere.cz[i]);
    reversed.time.push_back(pair.sphere.time[i]);
    reversed.zem.push_back(pair.sphere.zem[i]);
    reversed.size.push_back(pair.sphere.size[i]);
    reversed.wavelength.push_back(pair.sphere.wavelength[i]);
    reversed.id.push_back(pair.sphere.id[i]);
  }
  const Output back = run(model, reversed, identity_frame(0.0), true, 1);
  CHECK(back.rows() == one.rows());
  // Bunch ids (not positions in the input) decide the random numbers: the same photons.
  std::vector<std::uint64_t> forward_ids, backward_ids;
  for (std::size_t i = 0; i < one.rows(); ++i)
    forward_ids.push_back(one.index[i]);
  for (std::size_t i = 0; i < back.rows(); ++i)
    backward_ids.push_back(n - 1 - back.index[i]);
  std::sort(forward_ids.begin(), forward_ids.end());
  std::sort(backward_ids.begin(), backward_ids.end());
  CHECK(forward_ids == backward_ids);
}

// A photon that starts between the primary and the focal surface, moving down, takes the
// sequential path: the dish reflects it up to the focus.
void start_inside_moving_towards_primary() {
  const auto model = single_facet(10.0, 0.6);
  const Vec3 start{0.3, 0.0, 5.0};
  const auto t = trace_photon(model, {start, {0, 0, -1}}, 400.0, 0.0, 1.0, 0, 1);
  CHECK(t.status == Status::detected);
  CHECK(t.primary_surface == 1);
  CHECK(std::abs(t.primary_point.z - 0.00225) < 1e-3); // Sag of the 20 m sphere at r = 0.3 m.
  CHECK_NEAR(t.position.z, 10.0, 1e-9);
  CHECK(std::hypot(t.position.x, t.position.y) < 0.01); // Focused on the axis.
  // The same start as a sphere row (cz = -1, z = 5) gives the same result.
  Columns row;
  row.x = {0.3F};
  row.y = {0.0F};
  row.z = {5.0F};
  row.cx = {0.0F};
  row.cy = {0.0F};
  row.cz = {-1.0F};
  row.time = {0.0F};
  row.zem = {0.0F};
  row.size = {1.0F};
  row.wavelength = {400.0F};
  row.id = {0};
  const Output out = run(model, row, identity_frame(0.0), true, 1);
  CHECK(out.rows() == 1);
  CHECK_NEAR(out.position[2], t.position.z, 1e-9);
  CHECK_NEAR(out.position[0], t.position.x, 1e-6);
}

// A photon that starts between the primary and the focal surface, moving up, takes the first-hit
// path and reaches the focal surface from below: direct_focal, not detected.
void start_inside_moving_towards_camera() {
  const auto model = single_facet(10.0, 0.6);
  const auto t = trace_photon(model, {{0.3, 0.0, 5.0}, {0, 0, 1}}, 400.0, 0.0, 1.0, 0, 1, true,
                              false, true, true, 0, true);
  CHECK(t.status == Status::direct_focal);
  CHECK(t.terminal_surface == 2);
  CHECK_NEAR(t.position.z, 10.0, 1e-9);
  // As a sphere row with cz = +1 (first-hit branch by the traced direction).
  Columns row;
  row.x = {0.3F};
  row.y = {0.0F};
  row.z = {5.0F};
  row.cx = {0.0F};
  row.cy = {0.0F};
  row.cz = {1.0F};
  row.time = {0.0F};
  row.zem = {0.0F};
  row.size = {1.0F};
  row.wavelength = {400.0F};
  row.id = {0};
  // direct_focal is a loss: no detected row; the bunch trace yields nothing.
  CHECK(run(model, row, identity_frame(0.0), true, 1).rows() == 0);
}

// Horizontal photons (cz = 0 in the ground frame) on a tilted telescope: the traced direction has
// z < 0, so they take the sequential path and hit the dish front from the side.
void horizontal_rows_on_tilted_telescope() {
  const auto model = single_facet(10.0, 0.6);
  // Ground-to-telescope rotation about y by 30 deg: the traced z of (1, 0, 0) is -0.5.
  const double phi = 30.0 * 3.141592653589793 / 180.0;
  BunchFrame frame;
  const double c = std::cos(phi), s = std::sin(phi);
  frame.rotation[0] = c;
  frame.rotation[2] = s;
  frame.rotation[4] = 1.0;
  frame.rotation[6] = -s;
  frame.rotation[8] = c;
  const std::array<double, 3> d_ground{1.0, 0.0, 0.0};
  const std::array<double, 3> d_tel{c * d_ground[0] + s * d_ground[2], d_ground[1],
                                    -s * d_ground[0] + c * d_ground[2]};
  CHECK(d_tel[2] < -0.4);
  // Target on the dish front near the axis; start 3 m before it along the traced direction.
  const Vec3 target{0.2, 0.05, 0.0};
  const Vec3 start_tel{target.x - 3.0 * d_tel[0], target.y - 3.0 * d_tel[1],
                       target.z - 3.0 * d_tel[2]};
  // Inverse rotation (rows of R^T): the ground point of the start.
  const double gx = c * start_tel.x - s * start_tel.z;
  const double gy = start_tel.y;
  const double gz = s * start_tel.x + c * start_tel.z;
  Columns row;
  row.x = {static_cast<float>(gx)};
  row.y = {static_cast<float>(gy)};
  row.z = {static_cast<float>(gz)};
  row.cx = {1.0F};
  row.cy = {0.0F};
  row.cz = {0.0F};
  row.time = {0.0F};
  row.zem = {0.0F};
  row.size = {1.0F};
  row.wavelength = {400.0F};
  row.id = {0};
  const Output out = run(model, row, frame, true, 1);
  CHECK(out.rows() == 1);
  // The sequential trace of the same ray (double precision start from the float row).
  const Vec3 start{c * static_cast<double>(row.x[0]) + s * static_cast<double>(row.z[0]),
                   static_cast<double>(row.y[0]),
                   -s * static_cast<double>(row.x[0]) + c * static_cast<double>(row.z[0])};
  const auto seq = trace_photon(model, {start, normalized(Vec3{d_tel[0], d_tel[1], d_tel[2]})},
                                400.0, 0.0, 1.0, 0, 1, true, false, true, true, 0, false);
  CHECK(seq.status == Status::detected);
  CHECK(seq.primary_surface == 1);
  CHECK(std::abs(out.position[2] - seq.position.z) < 1e-6);
  CHECK(std::abs(out.position[0] - seq.position.x) < 1e-6);
}

// Rows with cz = 0 need z; without z they are still rejected (plane rows).
// The output path is the input path plus the traced path: a photon from (0, 0, 30) moving down
// reaches the dish at z = 0 (30 m) and the focus at z = 10 (10 m), so its path is 40 m more.
void sphere_row_path_is_carried() {
  const auto model = single_facet(10.0, 0.6);
  const float x = 0.0F, y = 0.0F, z = 30.0F, cx = 0.0F, cy = 0.0F, cz = -1.0F, t = 0.0F;
  const float zem = 0.0F, size = 1.0F, nm = 400.0F;
  const double path_in = 5.0;
  const std::uint64_t id = 0;
  BunchBlock in{&x, &y, &cx, &cy, &t, &zem, &size, &nm, &id, 1, &cz, &z, &path_in};
  ThreadPool pool(1);
  double position[3], direction[3], time_ns = 0.0, weight = 0.0, wavelength = 0.0, path = 0.0;
  std::uint64_t index = 0;
  std::uint32_t sub = 0;
  BunchResult out{position, direction, &time_ns, &weight, &wavelength, &index, &sub, &path};
  const std::size_t written = trace_bunches(model, in, identity_frame(0.0), BunchOptions{},
                                            {1, 1, 1}, pool, [&](std::size_t) { return out; });
  CHECK(written == 1);
  CHECK_NEAR(path, path_in + 40.0, 1e-9);
  CHECK_NEAR(position[2], 10.0, 1e-9);
}

void horizontal_rows_need_a_sphere_start() {
  const auto model = single_facet(10.0, 0.6);
  Columns row;
  row.x = {0.1F};
  row.y = {0.0F};
  row.z = {0.0F};
  row.cx = {1.0F};
  row.cy = {0.0F};
  row.cz = {0.0F};
  row.time = {0.0F};
  row.zem = {0.0F};
  row.size = {1.0F};
  row.wavelength = {400.0F};
  row.id = {0};
  ThreadPool pool(1);
  const auto allocate = [](std::size_t) { return BunchResult{}; };
  BunchBlock plane = row.block(false);
  bool rejected = false;
  try {
    trace_bunches(model, plane, identity_frame(kPlaneHeight), BunchOptions{}, {1, 1, 1}, pool,
                  allocate);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
  // A start height without the direction column is invalid.
  BunchBlock no_cz = row.block(true);
  no_cz.cz = nullptr;
  rejected = false;
  try {
    trace_bunches(model, no_cz, identity_frame(0.0), BunchOptions{}, {1, 1, 1}, pool, allocate);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
}

} // namespace

void all_tests() {
  sphere_rows_match_plane_rows();
  sphere_rows_invariance();
  start_inside_moving_towards_primary();
  start_inside_moving_towards_camera();
  horizontal_rows_on_tilted_telescope();
  horizontal_rows_need_a_sphere_start();
  sphere_row_path_is_carried();
}

int main() { return run_checks(all_tests); }
