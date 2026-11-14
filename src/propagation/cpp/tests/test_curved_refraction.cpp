// Spherical-shell ray bending and air delay: flat limit against bend() and
// straight_path_air_delay_ns, an independent Simpson reference on the real Earth radius, physics
// sanity, the interpolated tables against direct integration, and input rejection.
#include "check.hpp"

#include "obdeect/propagation/curved_refraction.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/timing.hpp"

#include "obdeect/stream/constants.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using obdeect::propagation::AtmosphereProfile;
using obdeect::propagation::curved_ray_integrals;
using obdeect::propagation::CurvedRayIntegrals;
using obdeect::propagation::CurvedRefractionTables;
using obdeect::propagation::RefractionTables;
using obdeect::propagation::speed_of_light_m_per_ns;
using obdeect::propagation::straight_path_air_delay_ns;

constexpr double kPi = obdeect::stream::pi;
constexpr double kEarthRadius = obdeect::stream::earth_radius_m;
constexpr double kFlatRadius = 1e12; // plane-parallel limit
constexpr double kTelescope = 2180.; // observation level of the fixture [m a.s.l.]
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

const std::string kFixture =
    std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/atmprof_ecmwf_north_winter_fixed.dat";

double radians(double degrees) { return degrees * kPi / 180.; }

template <class Function> bool throws_invalid(Function function) {
  try {
    function();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

// Test-local integrals: composite Simpson in u = sqrt(x), x = r - r_t, with n intervals. The
// integrand has a limit (not a sample) at u = 0 for theta_s = 90 degrees, so the first node is
// taken at a tiny u. Written with other algebra than the kernel: cos theta_vac^2 from
// (r - b)(r + b), cos theta_ref^2 from (n r - K)(n r + K).
CurvedRayIntegrals reference_integrals(const AtmosphereProfile &profile, double R, double h_t,
                                       double theta, double h_e, int n) {
  const double rt = R + h_t;
  const double s = std::sin(theta), c = std::cos(theta);
  const double b = rt * s;
  const double one_minus_s = c * c / (1. + s); // 1 - sin(theta), stable near 90 degrees
  const double ne = 1. + profile.index_minus_one_cm(100. * h_e);
  const double K = ne * b;
  const double U = std::sqrt(h_e - h_t);
  const double du = U / n;
  double sums[3] = {0., 0., 0.};
  for (int k = 0; k <= n; ++k) {
    const double u = (k == 0) ? 1e-12 * U : k * du;
    const double x = u * u;
    const double r = rt + x;
    const double nu = profile.index_minus_one_cm(100. * (h_t + x));
    const double nn = 1. + nu;
    const double r_minus_b = x + rt * one_minus_s;
    const double cos_vac = std::sqrt(r_minus_b * (r + b)) / r;
    const double sin_vac = b / r;
    const double nr_minus_K = nn * r - K;
    const double cos_ref = std::sqrt(nr_minus_K * (nn * r + K)) / (nn * r);
    const double sin_ref = K / (nn * r);
    const double f[3] = {rt * (sin_ref / cos_ref - sin_vac / cos_vac) / r,
                         nn / cos_ref - 1. / cos_vac, nu / cos_vac};
    const double weight = (k == 0 || k == n) ? 1. : ((k % 2) ? 4. : 2.);
    for (int q = 0; q < 3; ++q) {
      sums[q] += weight * 2. * u * f[q]; // dx = 2 u du
    }
  }
  const double scale = du / 3.;
  const double c_ns = speed_of_light_m_per_ns;
  return {sums[0] * scale, sums[1] * scale / c_ns, sums[2] * scale / c_ns};
}

// Absolute floor or relative tolerance, whichever is larger.
bool close(double value, double expected, double relative, double floor) {
  return std::abs(value - expected) <= std::max(relative * std::abs(expected), floor);
}

// Deterministic 64-bit LCG (Knuth MMIX constants), uniform in [0, 1).
class Lcg {
public:
  explicit Lcg(std::uint64_t seed) : state_(seed) {}
  double next() {
    state_ = state_ * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state_ >> 11) * (1.0 / 9007199254740992.0);
  }

private:
  std::uint64_t state_;
};

void test_flat_limit_against_bend(const AtmosphereProfile &profile) {
  const RefractionTables tables(profile, kTelescope);
  for (const double theta_deg : {0., 20., 45., 60., 75.}) {
    const double theta = radians(theta_deg);
    for (const double above_km : {1., 5., 20., 60.}) {
      const double h_e = kTelescope + 1000. * above_km;
      const CurvedRayIntegrals flat =
          curved_ray_integrals(profile, kFlatRadius, kTelescope, theta, h_e);

      // The ray bend() follows: emitted with the straight-line direction, as at the telescope
      // level.
      double cx = std::sin(theta), cy = 0., x = 0., y = 0., t = 0.;
      if (!obdeect::propagation::bend(tables, h_e, cx, cy, x, y, t)) {
        continue; // too close to the horizon for the series
      }
      const double vacuum_ns = (h_e - kTelescope) / (speed_of_light_m_per_ns * std::cos(theta));
      CHECK(close(flat.shift_m, x, 1e-4, 1e-3));
      CHECK(close(flat.excess_bent_ns, t - vacuum_ns, 1e-4, 1e-3));
      if (theta_deg == 0.) {
        CHECK(flat.shift_m == 0.);
        const double vertical = straight_path_air_delay_ns(profile, kTelescope, h_e, 0.);
        CHECK(close(flat.excess_bent_ns, vertical, 1e-6, 1e-9));
        CHECK(close(flat.excess_straight_ns, vertical, 1e-6, 1e-9));
      }
    }
  }
}

void test_flat_limit_straight_delay(const AtmosphereProfile &profile) {
  for (const double theta_deg : {0., 20., 45., 60., 75.}) {
    const double theta = radians(theta_deg);
    for (const double above_km : {1., 5., 20., 60.}) {
      const double h_e = kTelescope + 1000. * above_km;
      const CurvedRayIntegrals flat =
          curved_ray_integrals(profile, kFlatRadius, kTelescope, theta, h_e);
      const double expected = straight_path_air_delay_ns(profile, kTelescope, h_e, theta);
      CHECK(close(flat.excess_straight_ns, expected, 1e-6, 1e-12));
    }
  }
}

void test_reference_real_radius(const AtmosphereProfile &profile) {
  for (const double theta_deg : {0., 30., 60., 80., 85., 89., 90.}) {
    const double theta = radians(theta_deg);
    for (const double h_e_km : {3., 10., 30., 80.}) {
      const double h_e = 1000. * h_e_km;
      const CurvedRayIntegrals kernel =
          curved_ray_integrals(profile, kEarthRadius, kTelescope, theta, h_e);
      const CurvedRayIntegrals ref =
          reference_integrals(profile, kEarthRadius, kTelescope, theta, h_e, 200000);
      CHECK(close(kernel.shift_m, ref.shift_m, 1e-6, 1e-6));
      CHECK(close(kernel.excess_bent_ns, ref.excess_bent_ns, 1e-6, 1e-9));
      CHECK(close(kernel.excess_straight_ns, ref.excess_straight_ns, 1e-6, 1e-9));
    }
  }
}

void test_physics(const AtmosphereProfile &profile) {
  for (const double theta_deg : {0., 20., 30., 45., 60., 80., 89., 90.}) {
    const double theta = radians(theta_deg);
    for (const double h_e_km : {3., 10., 30., 80.}) {
      const CurvedRayIntegrals v =
          curved_ray_integrals(profile, kEarthRadius, kTelescope, theta, 1000. * h_e_km);
      CHECK(std::isfinite(v.shift_m) && std::isfinite(v.excess_bent_ns) &&
            std::isfinite(v.excess_straight_ns));
      CHECK(v.shift_m <= 0.);
      CHECK(v.excess_straight_ns >= 0.);
      // The bent-ray excess is non-negative for small zenith angles only: in the flat limit it is
      // proportional to nu (cos^2 - sin^2) / cos^3, which is negative beyond 45 degrees.
      if (theta_deg <= 30.) {
        CHECK(v.excess_bent_ns >= 0.);
      }
    }
  }
  // Zero shift at zenith angle zero, and the vertical delay equals the straight-line delay.
  const CurvedRayIntegrals vertical =
      curved_ray_integrals(profile, kEarthRadius, kTelescope, 0., 20000.);
  CHECK(vertical.shift_m == 0.);
  CHECK(close(vertical.excess_bent_ns, vertical.excess_straight_ns, 1e-9, 1e-12));

  // Constant above the profile top: emission heights above the top give the values at the top.
  const double top = profile.top_m();
  const CurvedRayIntegrals at_top =
      curved_ray_integrals(profile, kEarthRadius, kTelescope, radians(70.), top);
  const CurvedRayIntegrals above =
      curved_ray_integrals(profile, kEarthRadius, kTelescope, radians(70.), top + 5000.);
  CHECK(above.shift_m == at_top.shift_m);
  CHECK(above.excess_bent_ns == at_top.excess_bent_ns);
  CHECK(above.excess_straight_ns == at_top.excess_straight_ns);

  // Curvature at 80 degrees and an emission at 20 km a.s.l.: the shift differs from the flat value
  // by more than 1 %.
  const CurvedRayIntegrals curved80 =
      curved_ray_integrals(profile, kEarthRadius, kTelescope, radians(80.), 20000.);
  const CurvedRayIntegrals flat80 =
      curved_ray_integrals(profile, kFlatRadius, kTelescope, radians(80.), 20000.);
  std::fprintf(stderr, "info: shift at 80 deg, 20 km: curved %.6f m, flat %.6f m (ratio %.6f)\n",
               curved80.shift_m, flat80.shift_m, curved80.shift_m / flat80.shift_m);
  CHECK(std::abs(curved80.shift_m - flat80.shift_m) > 0.01 * std::abs(flat80.shift_m));

  // Finite at 90 degrees.
  const CurvedRayIntegrals horizon =
      curved_ray_integrals(profile, kEarthRadius, kTelescope, kPi / 2., 20000.);
  CHECK(std::isfinite(horizon.shift_m) && std::isfinite(horizon.excess_bent_ns) &&
        std::isfinite(horizon.excess_straight_ns));
}

void test_tables_against_direct(const AtmosphereProfile &profile) {
  const CurvedRefractionTables tables(profile, kEarthRadius, kTelescope);
  CHECK(tables.earth_radius_m() == kEarthRadius);
  CHECK(tables.telescope_altitude_m() == kTelescope);
  CHECK(tables.profile().levels() == profile.levels());
  CHECK(tables.theta_nodes() > 2 && tables.height_nodes() > 2);

  // Pseudo-random points over the full range of zenith angles and emission heights.
  Lcg random(20261114);
  std::size_t worst_shift_index = 0;
  double worst_shift_ratio = 0.;
  for (std::size_t k = 0; k < 2000; ++k) {
    const double theta = random.next() * kPi / 2.;
    const double h_e = kTelescope + random.next() * (profile.top_m() - kTelescope);
    CurvedRayIntegrals looked{};
    CHECK(tables.lookup(theta, h_e, looked));
    const CurvedRayIntegrals direct =
        curved_ray_integrals(profile, kEarthRadius, kTelescope, theta, h_e);
    const double shift_ratio =
        std::abs(looked.shift_m - direct.shift_m) / std::max(0.01, 1e-3 * std::abs(direct.shift_m));
    if (shift_ratio > worst_shift_ratio) {
      worst_shift_ratio = shift_ratio;
      worst_shift_index = k;
    }
    CHECK(shift_ratio <= 1.);
    CHECK(std::abs(looked.excess_bent_ns - direct.excess_bent_ns) <=
          0.01 + 1e-3 * std::abs(direct.excess_bent_ns));
    CHECK(std::abs(looked.excess_straight_ns - direct.excess_straight_ns) <=
          0.01 + 1e-3 * std::abs(direct.excess_straight_ns));
  }
  std::fprintf(stderr, "info: table vs direct, worst shift error ratio %.3f at sample %zu\n",
               worst_shift_ratio, worst_shift_index);

  // Exact at grid corners: zenith angle 90 degrees at the telescope level and at the top of the
  // grid.
  const double grid_top = std::min(profile.top_m(), 120000.);
  for (const double h_e : {kTelescope, grid_top}) {
    CurvedRayIntegrals looked{};
    CHECK(tables.lookup(kPi / 2., h_e, looked));
    const CurvedRayIntegrals direct =
        curved_ray_integrals(profile, kEarthRadius, kTelescope, kPi / 2., h_e);
    CHECK(close(looked.shift_m, direct.shift_m, 1e-12, 1e-12));
    CHECK(close(looked.excess_bent_ns, direct.excess_bent_ns, 1e-12, 1e-12));
    CHECK(close(looked.excess_straight_ns, direct.excess_straight_ns, 1e-12, 1e-12));
  }

  // Lookup rejections and clamping.
  CurvedRayIntegrals out{};
  CHECK(!tables.lookup(0.3, kTelescope - 1., out));
  CHECK(!tables.lookup(kNan, 20000., out));
  CHECK(!tables.lookup(0.3, kNan, out));
  CHECK(!tables.lookup(kInf, 20000., out));
  CHECK(!tables.lookup(-0.01, 20000., out));
  CHECK(!tables.lookup(kPi / 2. + 0.01, 20000., out));
  CHECK(tables.lookup(-1e-13, 20000., out));
  CHECK(tables.lookup(kPi / 2. + 1e-13, 20000., out));
  CurvedRayIntegrals at_top{}, above_top{};
  CHECK(tables.lookup(1.0, grid_top, at_top));
  CHECK(tables.lookup(1.0, grid_top + 5000., above_top));
  CHECK(above_top.shift_m == at_top.shift_m);
  CHECK(above_top.excess_bent_ns == at_top.excess_bent_ns);
  CHECK(above_top.excess_straight_ns == at_top.excess_straight_ns);
}

void test_rejections(const AtmosphereProfile &profile) {
  const double theta = 0.5, h_e = 20000.;
  CHECK(throws_invalid(
      [&] { static_cast<void>(curved_ray_integrals(profile, 0., kTelescope, theta, h_e)); }));
  CHECK(throws_invalid([&] {
    static_cast<void>(curved_ray_integrals(profile, -kEarthRadius, kTelescope, theta, h_e));
  }));
  CHECK(throws_invalid(
      [&] { static_cast<void>(curved_ray_integrals(profile, kNan, kTelescope, theta, h_e)); }));
  CHECK(throws_invalid(
      [&] { static_cast<void>(curved_ray_integrals(profile, kInf, kTelescope, theta, h_e)); }));
  CHECK(throws_invalid([&] {
    static_cast<void>(curved_ray_integrals(profile, kEarthRadius, kTelescope, -0.1, h_e));
  }));
  CHECK(throws_invalid([&] {
    static_cast<void>(
        curved_ray_integrals(profile, kEarthRadius, kTelescope, kPi / 2. + 0.01, h_e));
  }));
  CHECK(throws_invalid([&] {
    static_cast<void>(curved_ray_integrals(profile, kEarthRadius, kTelescope, kNan, h_e));
  }));
  CHECK(throws_invalid([&] {
    static_cast<void>(
        curved_ray_integrals(profile, kEarthRadius, kTelescope, theta, kTelescope - 1.));
  }));
  CHECK(throws_invalid([&] {
    static_cast<void>(curved_ray_integrals(profile, kEarthRadius, kTelescope, theta, kNan));
  }));
  CHECK(throws_invalid(
      [&] { static_cast<void>(curved_ray_integrals(profile, kEarthRadius, kNan, theta, h_e)); }));
  CHECK(throws_invalid(
      [&] { static_cast<void>(curved_ray_integrals(profile, kEarthRadius, -1., theta, h_e)); }));
  CHECK(throws_invalid([&] { CurvedRefractionTables bad(profile, 0., kTelescope); }));
  CHECK(throws_invalid([&] { CurvedRefractionTables bad(profile, kEarthRadius, kNan); }));
}

} // namespace

int main() {
  return run_checks([] {
    const AtmosphereProfile profile = AtmosphereProfile::from_file(kFixture);
    test_flat_limit_against_bend(profile);
    test_flat_limit_straight_delay(profile);
    test_reference_real_radius(profile);
    test_physics(profile);
    test_tables_against_direct(profile);
    test_rejections(profile);
  });
}
