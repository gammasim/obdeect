#include "check.hpp"
#include "sky_stream_check.hpp"

#include "obdeect/sources/background.hpp"
#include "obdeect/sources/programs/night-sky.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <variant>
#include <vector>

using namespace obdeect::sources;
using namespace obdeect::sources::programs;
using obdeect::stream::Column;

namespace {

constexpr double deg = pi / 180.0;

// E[z] of the sphere point: projected radius r = sqrt(U) (area-uniform disc), z = sqrt(1 - r^2).
// Analytic value: integral of sqrt(1 - u) over [0, 1] = 2/3.
double mean_zi() {
  const int n = 2000000;
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    const double u = (i + 0.5) / n;
    const double r = std::sqrt(u);
    s += std::sqrt(1.0 - r * r);
  }
  return s / n;
}

// Kolmogorov-Smirnov statistic of sorted values against the uniform law on [0, 1].
double ks_uniform(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const double n = static_cast<double>(v.size());
  double d = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    const double x = std::clamp(v[i], 0.0, 1.0);
    d = std::max(d,
                 std::max((static_cast<double>(i) + 1.0) / n - x, x - static_cast<double>(i) / n));
  }
  return d;
}

struct Moments {
  double mean{};
  double var{};
};

Moments moments_of(const std::vector<double> &v) {
  double s = 0.0, s2 = 0.0;
  for (const double x : v) {
    s += x;
    s2 += x * x;
  }
  const double n = static_cast<double>(v.size());
  const double mean = s / n;
  return {mean, (s2 - n * mean * mean) / (n - 1.0)};
}

// Rows of one event of a background plan, written in full through fill().
obdeect::stream::PhotonChunk rows_of(const NsbOptions &o, std::uint32_t event) {
  const auto plan = plan_night_sky(o);
  const auto &em = plan.cases[0][0];
  obdeect::stream::PhotonChunk chunk;
  em.fill(event, 0, static_cast<std::size_t>(em.count(event)), chunk);
  return chunk;
}

// Photons are emitted into the viewing cone and start on the telescope sphere.
void cone_and_sphere() {
  NsbOptions o;
  o.photons = 6e4;
  o.bunch_size = 2.0;
  o.spectrum = Fixed{400.0};
  const auto plan = plan_night_sky(o);
  CHECK(plan.cases[0][0].count(1) == 30000U); // mean 3e4, no fractional part
  const auto c = rows_of(o, 1);
  const double cos_vc = std::cos(10.0 * deg);
  std::vector<double> cos_theta;
  for (std::size_t i = 0; i < c.rows(); ++i) {
    const Vec3 d{c[Column::cx][i], c[Column::cy][i], c[Column::cz][i]};
    CHECK(-d.z >= cos_vc - 1e-6);
    CHECK_NEAR(norm(d), 1.0, 1e-6);
    CHECK_NEAR(norm(Vec3{c[Column::x][i], c[Column::y][i], c[Column::z][i]}), 15.0, 1e-3);
    CHECK_NEAR(c[Column::w][i], 2.0, 1e-6);
    CHECK_NEAR(c[Column::wavelength][i], 400.0, 1e-4);
    cos_theta.push_back(-d.z);
  }
  // cos(theta) is uniform in [cos(vc), 1] for the cone (solid-angle sampling).
  const auto m = moments_of(cos_theta);
  const double n = static_cast<double>(cos_theta.size());
  const double width = 1.0 - cos_vc;
  CHECK_NEAR(m.mean, 0.5 * (1.0 + cos_vc), 4.0 * std::sqrt(width * width / 12.0 / n));
  CHECK_NEAR(m.var, width * width / 12.0, 4.0 * (width * width / 12.0) * std::sqrt(0.8 / n));
}

// Entry points are uniform over the projected disc of the sphere: for each direction d the
// perpendicular radius ri = |u - (u.d) d| has ri^2 uniform in [0, 1] (isotropic flux through
// the telescope sphere), and E[zi] = 2/3.
void projected_area_uniform() {
  BackgroundSource s;
  s.centre = {0.0, 0.0, 0.0};
  s.radius_m = 15.0;
  s.pointing = {0.0, 0.0, -1.0};
  s.viewcone_rad = pi / 2.0;
  const std::uint64_t n = 40000;
  std::vector<double> r2;
  std::vector<double> zl;
  for (std::uint64_t i = 0; i < n; ++i) {
    const Photon ph = sample_background(s, 11, 2, 1, i);
    const Vec3 u = ph.position * (1.0 / s.radius_m);
    const double ud = dot(u, ph.direction);
    r2.push_back(1.0 - ud * ud);
    zl.push_back(-ud);
  }
  const auto m = moments_of(r2);
  const double nd = static_cast<double>(n);
  CHECK_NEAR(m.mean, 0.5, 4.0 * std::sqrt(1.0 / 12.0 / nd));
  CHECK_NEAR(m.var, 1.0 / 12.0, 4.0 * (1.0 / 12.0) * std::sqrt(0.8 / nd));
  // KS critical value at alpha = 0.001: 1.95 / sqrt(n).
  CHECK(ks_uniform(r2) < 1.95 / std::sqrt(nd));
  CHECK_NEAR(moments_of(zl).mean, 2.0 / 3.0, 4.0 * std::sqrt(1.0 / 18.0 / nd));
  CHECK_NEAR(mean_zi(), 2.0 / 3.0, 1e-6);
}

// Emission time: -zi R / c (zi of the sphere point, pulse 0, delay 0). With a pulse width
// and a delay the mean is shifted by the delay only.
void timing_and_delay() {
  NsbOptions o;
  o.photons = 6e4;
  o.bunch_size = 2.0;
  o.spectrum = Fixed{400.0};
  o.pulse_width_ns = 4.0;
  o.pulse_delay_ns = 10.0;
  const auto c = rows_of(o, 1);
  std::vector<double> t;
  for (const float x : c[Column::t])
    t.push_back(x);
  const auto m = moments_of(t);
  const double n = static_cast<double>(t.size());
  const double expected = 10.0 - mean_zi() * 15.0 / speed_of_light_m_per_ns;
  CHECK_NEAR(m.mean, expected, 4.0 * std::sqrt(m.var / n));
}

// A tabulated spectrum is sampled from its density: uniform 300..600 nm.
void table_spectrum() {
  NsbOptions o;
  o.photons = 6e4;
  o.bunch_size = 2.0;
  o.spectrum = Tabulated{LinearDensity({300.0, 600.0}, {1.0, 1.0})};
  const auto c = rows_of(o, 1);
  std::vector<double> wl;
  for (const float x : c[Column::wavelength])
    wl.push_back(x);
  const auto m = moments_of(wl);
  const double n = static_cast<double>(wl.size());
  CHECK_NEAR(m.mean, 450.0, 4.0 * std::sqrt(7500.0 / n));
  CHECK_NEAR(m.var, 7500.0, 4.0 * 7500.0 * std::sqrt(0.8 / n));
}

// Tilted viewing directions: every photon starts on the hemisphere that faces its own
// direction, so it travels inward (dot(pos - centre, d) < 0). Its emission distance
// zi = -dot(u, d) has the mean of the hemisphere point distribution, E[zi].
void tilted_directions_enter_sphere() {
  const Vec3 pointings[] = {
      {0.5, 0.0, -0.866}, {0.3, -0.4, -0.866}, {0.0, 0.7, -0.714}, {-0.6, 0.6, -0.529}};
  const double cone_deg = 60.0;
  for (const Vec3 &pointing : pointings) {
    BackgroundSource s;
    s.centre = {1.0, -2.0, 3.0};
    s.radius_m = 15.0;
    s.pointing = normalised(pointing);
    s.viewcone_rad = cone_deg * deg;
    const std::uint64_t n = 20000;
    int inside = 0, outside_cone = 0;
    std::vector<double> zl;
    for (std::uint64_t i = 0; i < n; ++i) {
      const Photon ph = sample_background(s, 7, 1, 1, i);
      const Vec3 u = (ph.position - s.centre) * (1.0 / s.radius_m);
      if (dot(u, ph.direction) < 0.0)
        ++inside;
      if (dot(ph.direction, s.pointing) < std::cos(cone_deg * deg) - 1e-9)
        ++outside_cone;
      zl.push_back(-dot(u, ph.direction));
    }
    CHECK(inside == static_cast<int>(n));
    CHECK(outside_cone == 0);
    const auto m = moments_of(zl);
    CHECK_NEAR(m.mean, mean_zi(), 4.0 * std::sqrt(m.var / static_cast<double>(n)));
  }
}

// Stream: run, one telescope block on the sphere, photon total and events.
void stream_layout() {
  NsbOptions o;
  o.photons = 4e4;
  o.events = 2;
  o.spectrum = Fixed{400.0};
  const auto s = read_stream(write_to_string([&](std::ostream &out) { write_night_sky(out, o); }));
  CHECK(s.ended && s.events == 2);
  CHECK(s.types.front() == 1 && s.types.back() == 8);
  CHECK(std::get<std::int64_t>(s.run.at("n_telescopes")) == 1);
  CHECK_NEAR(std::get<std::vector<double>>(s.run.at("telescope_r_m")).at(0), 15.0, 1e-12);
  CHECK(s.photons.at(1) == 2 * 20000U);
}

// Emission altitude: zem = observation level + z (format.md, zem), set by the stream.
void emission_altitude() {
  NsbOptions o;
  o.photons = 4e4;
  o.altitude_m = 1000.0;
  o.spectrum = Fixed{400.0};
  const auto s = read_stream(write_to_string([&](std::ostream &out) { write_night_sky(out, o); }));
  CHECK(!s.chunks.empty());
  for (const auto &c : s.chunks) {
    CHECK((c.mask() & obdeect::stream::column_bit(Column::zem)) != 0U);
    if ((c.mask() & obdeect::stream::column_bit(Column::zem)) == 0U)
      continue;
    for (std::size_t i = 0; i < c.rows(); ++i)
      CHECK_NEAR(c[Column::zem][i] - c[Column::z][i], 1000.0, 1e-2);
  }
}

// Negative photon counts are rejected, not clamped to zero.
void rejects_negative_photons() {
  NsbOptions o;
  o.photons = -5.0;
  bool threw = false;
  try {
    (void)plan_night_sky(o);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
}

// Source parameters are checked when the plan is built, also with zero photons.
template <class Change> bool plan_throws(Change change) {
  NsbOptions o;
  o.photons = 0.0;
  change(o);
  try {
    (void)plan_night_sky(o);
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

void rejects_bad_source_parameters() {
  CHECK(plan_throws([](NsbOptions &o) { o.telescope_radius_m = -1.0; }));
  CHECK(plan_throws([](NsbOptions &o) { o.telescope_radius_m = 0.0; }));
  CHECK(plan_throws([](NsbOptions &o) { o.viewcone_deg = 200.0; }));
  CHECK(plan_throws([](NsbOptions &o) { o.viewcone_deg = 0.0; }));
  CHECK(plan_throws([](NsbOptions &o) { o.pulse_width_ns = -1.0; }));
  CHECK(plan_throws([](NsbOptions &o) { o.pulse_delay_ns = std::nan(""); }));
  CHECK(plan_throws([](NsbOptions &o) { o.bunch_size = 0.0; }));
  CHECK(plan_throws([](NsbOptions &o) { o.altitude_m = std::nan(""); }));
  CHECK(!plan_throws([](NsbOptions &) {}));
  // Viewcone 180 deg (the full sphere) is valid.
  CHECK(!plan_throws([](NsbOptions &o) { o.viewcone_deg = 180.0; }));
}

} // namespace

int main() {
  return run_checks([] {
    cone_and_sphere();
    projected_area_uniform();
    timing_and_delay();
    table_spectrum();
    stream_layout();
    tilted_directions_enter_sphere();
    emission_altitude();
    rejects_negative_photons();
    rejects_bad_source_parameters();
  });
}
