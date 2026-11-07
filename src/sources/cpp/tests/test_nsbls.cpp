#include "check.hpp"
#include "sky_stream_check.hpp"

#include "dlamb/programs/nsbls.hpp"

#include <cmath>
#include <cstdint>
#include <variant>
#include <vector>

using namespace dlamb;
using namespace dlamb::programs;
using gstrickd::Column;

namespace {

constexpr double deg = pi / 180.0;

// E[z] of the sphere point of the original: projected radius r = 1 - sqrt(U), z = sqrt(1 - r^2).
double mean_zi() {
  const int n = 2000000;
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    const double u = (i + 0.5) / n;
    const double r = 1.0 - std::sqrt(u);
    s += std::sqrt(1.0 - r * r);
  }
  return s / n;
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
gstrickd::PhotonChunk rows_of(const NsbOptions &o, std::uint32_t event) {
  const auto plan = plan_nsbls(o);
  const auto &em = plan.cases[0][0];
  gstrickd::PhotonChunk chunk;
  em.fill(event, 0, static_cast<std::size_t>(em.count(event)), chunk);
  return chunk;
}

// Photons are emitted into the viewing cone and start on the telescope sphere.
void cone_and_sphere() {
  NsbOptions o;
  o.photons = 6e4;
  o.bunch_size = 2.0;
  o.spectrum = Fixed{400.0};
  const auto plan = plan_nsbls(o);
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

// Stream: run, one telescope block on the sphere, photon total and events.
void stream_layout() {
  NsbOptions o;
  o.photons = 4e4;
  o.events = 2;
  o.spectrum = Fixed{400.0};
  const auto s = read_stream(write_to_string([&](std::ostream &out) { write_nsbls(out, o); }));
  CHECK(s.ended && s.events == 2);
  CHECK(s.types.front() == 1 && s.types.back() == 8);
  CHECK(std::get<std::int64_t>(s.run.at("n_telescopes")) == 1);
  CHECK_NEAR(std::get<std::vector<double>>(s.run.at("telescope_r_m")).at(0), 15.0, 1e-12);
  CHECK(s.photons.at(1) == 2 * 20000U);
}

} // namespace

int main() {
  return run_checks([] {
    cone_and_sphere();
    timing_and_delay();
    table_spectrum();
    stream_layout();
  });
}
