#include "check.hpp"
#include "sky_stream_check.hpp"

#include "dlamb/programs/laser-beam.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <variant>

using namespace dlamb;
using namespace dlamb::programs;
using gstrickd::Column;

namespace {

constexpr double deg = pi / 180.0;

// Telescope axis of ls-beam.cc: DirectionVect((90 - theta) deg, phi deg).
Vec3 telescope_axis(const LsBeamOptions &o) {
  return direction_from_altitude((90.0 - o.telescope_theta_deg) * deg, o.telescope_phi_deg * deg);
}

// Cosine of the angle between a point (seen from the telescope at the origin) and the axis.
double cos_to_axis(Vec3 p, const LsBeamOptions &o) { return dot(normalised(p), telescope_axis(o)); }

// The beam enters and leaves the cone of full opening angle fov: both ends lie on it.
void geometry_default() {
  const LsBeamOptions o;
  const BeamGeometry g = laser_beam_geometry(o);
  const double half = 0.5 * o.fov_deg * deg;
  CHECK(g.path_m > 0.0);
  CHECK_NEAR(cos_to_axis(g.start_m, o), std::cos(half), 1e-9);
  CHECK_NEAR(cos_to_axis(g.start_m + g.direction * g.path_m, o), std::cos(half), 1e-9);
  // The beam direction is the laser pointing: zenith 60 deg, azimuth 180 deg.
  CHECK_NEAR(g.direction.z, std::sin(30.0 * deg), 1e-12);
  CHECK_NEAR(g.direction.x, -std::cos(30.0 * deg), 1e-12);
}

// Default run: steps inside the cone, step count and time, primary photon sum, colours.
void default_steps() {
  const LsBeamOptions o;
  const auto plan = plan_laser_beam(o);
  const BeamGeometry g = laser_beam_geometry(o);
  const auto &steps = plan.cases[0];
  const double expected_steps = std::ceil(g.path_m / o.step_m);
  CHECK(std::abs(static_cast<double>(steps.size()) - expected_steps) <= 1.0);
  const double half_cos = std::cos(0.5 * o.fov_deg * deg);
  double sum_expected = 0.0;
  for (std::size_t k = 0; k < steps.size(); ++k) {
    const auto &e = std::get<Emitter>(steps[k].source.geometry);
    CHECK(cos_to_axis(e.position, o) >= half_cos - 1e-9);
    CHECK_NEAR(steps[k].time_offset_ns, static_cast<double>(k) * o.step_m / speed_of_light_m_per_ns,
               1e-9);
    sum_expected += o.primary_photons_per_m * o.step_m;
  }
  // Photons = bunches x bunch size: the expected sum of weights is the primary total.
  const double bunch = std::max(1.0, o.bunch_size);
  std::uint64_t events_sum = 0;
  for (const auto &st : steps)
    events_sum += st.count(1);
  const double sum_w = static_cast<double>(events_sum) * bunch;
  CHECK_NEAR(sum_w, sum_expected, 4.0 * bunch * std::sqrt(static_cast<double>(steps.size())) / 2.0);

  const auto s = read_stream(write_to_string([&](std::ostream &out) { write_laser_beam(out, o); }));
  CHECK(s.ended && s.events == 1);
  CHECK(std::get<std::string>(s.run.at("scattering_angular")) == "Rayleigh");
  CHECK_NEAR(std::get<double>(s.run.at("beam_path_m")), g.path_m, 1e-9);
  const Vec3 axis = g.direction;
  for (const auto &c : s.chunks) {
    for (std::size_t i = 0; i < c.rows(); ++i) {
      CHECK_NEAR(c[Column::cz][i], axis.z, 1e-6);
      CHECK_NEAR(c[Column::wavelength][i], 355.0, 1e-4);
      CHECK_NEAR(c[Column::w][i], bunch, 1e-6);
    }
  }
}

// Non-default: fov 20 deg, 0.5 m step, bunch size 1, Gauss pulse of sigma 0.5 ns. The
// photon times minus the path offset have the pulse's moments.
void non_default_pulse() {
  LsBeamOptions o;
  o.fov_deg = 20.0;
  o.step_m = 0.5;
  o.bunch_size = 1.0;
  o.primary_photons_per_m = 1e3;
  o.pulse = GaussianShape{0.0, 0.5};
  const auto plan = plan_laser_beam(o);
  const BeamGeometry g = laser_beam_geometry(o);
  const double half_cos = std::cos(0.5 * o.fov_deg * deg);
  double sum = 0.0, sum2 = 0.0;
  std::uint64_t n = 0;
  gstrickd::PhotonChunk chunk;
  for (const auto &st : plan.cases[0]) {
    CHECK(cos_to_axis(std::get<Emitter>(st.source.geometry).position, o) >= half_cos - 1e-9);
    const std::uint64_t rows = st.count(1);
    if (rows == 0)
      continue;
    st.fill(1, 0, static_cast<std::size_t>(rows), chunk);
    for (const float t : chunk[Column::t]) {
      const double residual = t - st.time_offset_ns;
      sum += residual;
      sum2 += residual * residual;
      ++n;
    }
  }
  CHECK(n > 10000);
  const double mean = sum / static_cast<double>(n);
  const double var = (sum2 - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1);
  CHECK_NEAR(mean, 0.0, 4.0 * 0.5 / std::sqrt(static_cast<double>(n)));
  // Gaussian: kurtosis 3, variance of the sample variance 2 sigma^4 / n.
  CHECK_NEAR(var, 0.25, 4.0 * 0.25 * std::sqrt(2.0 / static_cast<double>(n)));
  CHECK(g.path_m > 0.0);
}

// A beam that misses the cone is rejected.
void rejects_missing_beam() {
  LsBeamOptions o;
  o.laser_position_m = {0.0, 1000.0, 0.0};
  o.laser_theta_deg = 90.0; // horizontal beam along +x, 1 km off the telescope axis
  bool threw = false;
  try {
    (void)laser_beam_geometry(o);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
}

// The field of view must lie in (0, 90) deg: a negative cone is not accepted as its mirror.
void rejects_fov_outside_range() {
  const auto throws_for = [](double fov) {
    LsBeamOptions o;
    o.fov_deg = fov;
    try {
      (void)laser_beam_geometry(o);
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  CHECK(throws_for(-15.0));
  CHECK(throws_for(0.0));
  CHECK(throws_for(90.0));
  CHECK(throws_for(180.0));
  CHECK(throws_for(std::nan("")));
  CHECK(throws_for(std::numeric_limits<double>::infinity()));
  CHECK(!throws_for(15.0));
  CHECK(!throws_for(89.9));
}

} // namespace

int main() {
  return run_checks([] {
    geometry_default();
    default_steps();
    non_default_pulse();
    rejects_missing_beam();
    rejects_fov_outside_range();
  });
}
