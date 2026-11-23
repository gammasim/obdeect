// Wavelength dependence of refraction: Ciddor dry-air refractivity and its group index (analytic
// against a numerical derivative), the plane-parallel bend() with scale factors k and k_group
// (against numerical integration of the exact integrands), and the spherical lookup_scaled()
// (against direct quadrature with the n-1 column multiplied by k).
#include "check.hpp"

#include "obdeect/propagation/curved_refraction.hpp"
#include "obdeect/propagation/dispersion.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/timing.hpp"
#include "obdeect/stream/constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using obdeect::propagation::AtmosphereProfile;
using obdeect::propagation::bend;
using obdeect::propagation::curved_ray_integrals;
using obdeect::propagation::CurvedRayIntegrals;
using obdeect::propagation::CurvedRefractionTables;
using obdeect::propagation::group_refractivity_e8;
using obdeect::propagation::group_refractivity_ratio;
using obdeect::propagation::RefractionTables;
using obdeect::propagation::refractivity_e8;
using obdeect::propagation::refractivity_ratio;
using obdeect::propagation::speed_of_light_m_per_ns;
using obdeect::propagation::straight_path_air_delay_ns;

constexpr double kPi = obdeect::stream::pi;
constexpr double kEarthRadius = obdeect::stream::earth_radius_m;
constexpr double kObservation = 2180.; // observation level of the fixture [m a.s.l.]
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

// Fixture text with the n-1 column (4th column) multiplied by k; comments and other columns kept.
std::string scaled_profile_text(const std::string &path, double k) {
  std::ifstream file(path);
  if (!file) {
    throw std::runtime_error("cannot open " + path);
  }
  std::ostringstream out;
  std::string line;
  while (std::getline(file, line)) {
    const std::size_t first = line.find_first_not_of(" \t");
    if (first == std::string::npos || line[first] == '#') {
      out << line << '\n';
      continue;
    }
    std::istringstream fields(line);
    std::vector<std::string> tokens;
    std::string token;
    while (fields >> token) {
      tokens.push_back(token);
    }
    if (tokens.size() < 4) {
      throw std::runtime_error("short data line in " + path);
    }
    char buffer[64];
    std::snprintf(buffer, sizeof buffer, "%.9e", k * std::stod(tokens[3]));
    tokens[3] = buffer;
    for (std::size_t i = 0; i < tokens.size(); ++i) {
      out << tokens[i] << (i + 1 < tokens.size() ? ' ' : '\n');
    }
  }
  return out.str();
}

// Ciddor dry-air refractivity: value and derivative (numerical, central in ln lambda).
void test_ciddor(const AtmosphereProfile &) {
  // Standard air at 589.3 nm, (n-1) 1e8 = 27716 (n = 1.000277).
  const double n589 = refractivity_e8(589.3);
  std::fprintf(stderr, "info: (n-1) 1e8 at 589.3 nm = %.3f\n", n589);
  CHECK(n589 > 27700. && n589 < 27730.);

  // Ratios: exactly 1 at the reference; phase and group ratios decrease with wavelength.
  CHECK(refractivity_ratio(500., 500.) == 1.);
  CHECK(group_refractivity_ratio(500., 500.) > 1.);
  CHECK(refractivity_ratio(300., 500.) > refractivity_ratio(400., 500.));
  CHECK(refractivity_ratio(400., 500.) > refractivity_ratio(700., 500.));
  std::fprintf(stderr, "info: ref 500 nm: k(300)=%.6f k(400)=%.6f k(700)=%.6f\n",
               refractivity_ratio(300., 500.), refractivity_ratio(400., 500.),
               refractivity_ratio(700., 500.));
  std::fprintf(stderr, "info: ref 500 nm: kg(300)=%.6f kg(400)=%.6f kg(700)=%.6f\n",
               group_refractivity_ratio(300., 500.), group_refractivity_ratio(400., 500.),
               group_refractivity_ratio(700., 500.));

  // Group index n_g - 1 = (n-1) - dN/dln(lambda), against a central difference in ln(lambda).
  const double step = 1e-4;
  for (const double lambda : {250., 350., 500., 800., 1500.}) {
    const double up = refractivity_e8(lambda * std::exp(step));
    const double down = refractivity_e8(lambda * std::exp(-step));
    const double n = refractivity_e8(lambda);
    const double numerical = n - (up - down) / (2. * step);
    const double analytic = group_refractivity_e8(lambda);
    std::fprintf(stderr, "info: lambda %.0f nm: group analytic %.6f numerical %.6f rel diff %.2e\n",
                 lambda, analytic, numerical, std::abs(analytic - numerical) / analytic);
    CHECK(std::abs(analytic - numerical) <= 1e-7 * analytic);
  }

  // Validation: the wavelength range and non-finite values are rejected.
  for (const double bad : {199.9, 2000.1, 0., -500., kNan, kInf}) {
    CHECK(throws_invalid([&] { static_cast<void>(refractivity_e8(bad)); }));
    CHECK(throws_invalid([&] { static_cast<void>(group_refractivity_e8(bad)); }));
    CHECK(throws_invalid([&] { static_cast<void>(refractivity_ratio(bad, 500.)); }));
    CHECK(throws_invalid([&] { static_cast<void>(group_refractivity_ratio(500., bad)); }));
  }
  CHECK(!throws_invalid([&] { static_cast<void>(refractivity_e8(200.)); }));
  CHECK(!throws_invalid([&] { static_cast<void>(refractivity_e8(2000.)); }));
}

// Plane-parallel reference: composite Simpson in height of the exact integrands with the phase
// index n = 1 + k nu for the geometry (Snell's law) and the group index n_g = 1 + k_group nu for
// the time. Returns the horizontal displacement [m] and the travel time from the observation level
// [ns] for a photon emitted at height h_e [m] with zenith angle theta_e.
struct PlaneReference {
  double x_m;
  double t_ns;
};

PlaneReference plane_reference(const AtmosphereProfile &profile, double h_e, double theta_e,
                               double k, double k_group, int intervals) {
  const double s_e = std::sin(theta_e);
  const double nu_e = profile.refractive_index_minus_one(h_e);
  const double dh = (h_e - kObservation) / intervals;
  double sum_x = 0., sum_t = 0.;
  for (int i = 0; i <= intervals; ++i) {
    const double h = (i == intervals) ? h_e : kObservation + i * dh;
    const double nu = profile.refractive_index_minus_one(h);
    const double s = s_e * (1. + k * nu_e) / (1. + k * nu);
    const double c = std::sqrt(1. - s * s);
    const double weight = (i == 0 || i == intervals) ? 1. : ((i % 2) ? 4. : 2.);
    sum_x += weight * (s / c - std::tan(theta_e));
    sum_t += weight * (1. + k_group * nu) / c;
  }
  const double scale = (h_e - kObservation) / intervals / 3.;
  return {sum_x * scale, sum_t * scale / speed_of_light_m_per_ns};
}

void test_plane_scaled(const AtmosphereProfile &profile) {
  const RefractionTables tables(profile, kObservation);

  // Scale pairs (k, k_group): Ciddor ratios at 300, 400, 700 nm for 500 nm reference, and
  // synthetic pairs. Equal pairs test the series alone, unequal pairs the group term too.
  struct Pair {
    double k, k_group;
  };
  const std::vector<Pair> pairs = {
      {refractivity_ratio(300., 500.), group_refractivity_ratio(300., 500.)},
      {refractivity_ratio(400., 500.), group_refractivity_ratio(400., 500.)},
      {refractivity_ratio(700., 500.), group_refractivity_ratio(700., 500.)},
      {0.97, 0.97},
      {1.03, 1.03},
      {0.97, 0.98},
      {1.03, 1.02},
  };

  double worst_x_series = 0., worst_t_series = 0., worst_x_group = 0., worst_t_group = 0.;
  double worst_identity = 0.;
  int used = 0, total = 0;
  for (const Pair &pair : pairs) {
    for (const double theta_deg : {0., 20., 40., 60., 70., 75.}) {
      for (const double h_e : {10000., 20000., 30000.}) {
        const double theta = radians(theta_deg);
        ++total;
        double cx = std::sin(theta), cy = 0., x = 0., y = 0., t = 0.;
        if (!bend(tables, h_e, cx, cy, x, y, t, pair.k, pair.k_group)) {
          continue;
        }
        ++used;
        const PlaneReference ref =
            plane_reference(profile, h_e, theta, pair.k, pair.k_group, 20000);
        const double ex = std::abs(x - ref.x_m), et = std::abs(t - ref.t_ns);
        if (pair.k == pair.k_group) {
          worst_x_series = std::max(worst_x_series, ex);
          worst_t_series = std::max(worst_t_series, et);
        } else {
          worst_x_group = std::max(worst_x_group, ex);
          worst_t_group = std::max(worst_t_group, et);
        }

        // The group term is exactly (k_group - k) times the straight-path air delay.
        double cx0 = std::sin(theta), cy0 = 0., x0 = 0., y0 = 0., t0 = 0.;
        CHECK(bend(tables, h_e, cx0, cy0, x0, y0, t0, pair.k, pair.k));
        const double group_term = t - t0;
        const double expected =
            (pair.k_group - pair.k) * straight_path_air_delay_ns(profile, kObservation, h_e, theta);
        worst_identity = std::max(worst_identity, std::abs(group_term - expected));
      }
    }
  }
  std::fprintf(stderr,
               "info: plane scaled bend (%d of %d photons): series x %.3e m t %.3e ns; "
               "group x %.3e m t %.3e ns; group identity %.1e ns\n",
               used, total, worst_x_series, worst_t_series, worst_x_group, worst_t_group,
               worst_identity);
  CHECK(used * 10 >= total * 8);
  CHECK(worst_x_series <= 1e-6);
  CHECK(worst_t_series <= 1e-6);
  CHECK(worst_identity <= 1e-9);
  CHECK(worst_x_group <= 1e-6);
  CHECK(worst_t_group <= 8e-3);

  // k = k_group = 1 is bit-identical to the 7-argument bend().
  for (const double h_e : {5000., 20000.}) {
    for (const double theta_deg : {10., 50.}) {
      const double theta = radians(theta_deg);
      double ax = std::sin(theta), ay = 0.3 * std::sin(theta), axx = 0.1, ayy = -0.2, at = 5.;
      double bx = ax, by = ay, bxx = axx, byy = ayy, bt = at;
      const bool a = bend(tables, h_e, ax, ay, axx, ayy, at);
      const bool b = bend(tables, h_e, bx, by, bxx, byy, bt, 1., 1.);
      CHECK(a == b);
      CHECK(ax == bx && ay == by && axx == bxx && ayy == byy && at == bt);
    }
  }

  // Invalid scales and inputs: returns false and leaves the photon unchanged.
  for (const double bad : {0., -1., kNan, kInf}) {
    double cx = 0.5, cy = 0., x = 1., y = 2., t = 3.;
    CHECK(!bend(tables, 15000., cx, cy, x, y, t, bad, 1.));
    CHECK(!bend(tables, 15000., cx, cy, x, y, t, 1., bad));
    CHECK(cx == 0.5 && x == 1. && t == 3.);
  }
}

// Spherical reference: composite Simpson in u = sqrt(x), x = r - r_t, of the Bouguer-law integrands
// with the phase index 1 + k nu and the group index 1 + k_group nu in the time integrand.
CurvedRayIntegrals reference_scaled(const AtmosphereProfile &profile, double h_t, double theta,
                                    double h_e, double k, double k_group, int n) {
  const double rt = kEarthRadius + h_t;
  const double s = std::sin(theta), c = std::cos(theta);
  const double b = rt * s;
  const double one_minus_s = c * c / (1. + s); // 1 - sin(theta), stable near 90 degrees
  const double ne = 1. + k * profile.index_minus_one_cm(100. * h_e);
  const double K = ne * b;
  const double U = std::sqrt(h_e - h_t);
  const double du = U / n;
  double sums[3] = {0., 0., 0.};
  for (int i = 0; i <= n; ++i) {
    const double u = (i == 0) ? 1e-12 * U : i * du;
    const double x = u * u;
    const double r = rt + x;
    const double nu0 = profile.index_minus_one_cm(100. * (h_t + x));
    const double nu = k * nu0;
    const double nn = 1. + nu;
    const double cos_vac = std::sqrt((x + rt * one_minus_s) * (r + b)) / r;
    const double sin_vac = b / r;
    const double cos_ref = std::sqrt((nn * r - K) * (nn * r + K)) / (nn * r);
    const double sin_ref = K / (nn * r);
    const double f[3] = {rt * (sin_ref / cos_ref - sin_vac / cos_vac) / r,
                         (1. + k_group * nu0) / cos_ref - 1. / cos_vac, k_group * nu0 / cos_vac};
    const double weight = (i == 0 || i == n) ? 1. : ((i % 2) ? 4. : 2.);
    for (int q = 0; q < 3; ++q) {
      sums[q] += weight * 2. * u * f[q];
    }
  }
  const double scale = du / 3.;
  return {sums[0] * scale, sums[1] * scale / speed_of_light_m_per_ns,
          sums[2] * scale / speed_of_light_m_per_ns};
}

void test_spherical_scaled(const AtmosphereProfile &profile) {
  const CurvedRefractionTables tables(profile, kEarthRadius, kObservation, false, 2);

  // The reference against the kernel on the same profile (k = 1): validates the quadrature.
  for (const double theta_deg : {0., 45., 85.}) {
    const double theta = radians(theta_deg);
    const CurvedRayIntegrals kernel =
        curved_ray_integrals(profile, kEarthRadius, kObservation, theta, 30000.);
    const CurvedRayIntegrals ref =
        reference_scaled(profile, kObservation, theta, 30000., 1., 1., 20000);
    CHECK(std::abs(kernel.shift_m - ref.shift_m) <= 1e-6 * std::abs(kernel.shift_m) + 1e-9);
    CHECK(std::abs(kernel.excess_bent_ns - ref.excess_bent_ns) <=
          1e-6 * std::abs(kernel.excess_bent_ns) + 1e-9);
  }

  // k = k_group = 1: the lookup is the table itself, exactly.
  {
    CurvedRayIntegrals base{}, scaled{};
    CHECK(tables.lookup(0.7, 20000., base));
    CHECK(tables.lookup_scaled(0.7, 20000., 1., 1., scaled));
    CHECK(base.shift_m == scaled.shift_m && base.excess_bent_ns == scaled.excess_bent_ns &&
          base.excess_straight_ns == scaled.excess_straight_ns);
  }

  // Equal scales: the profile multiplied by k in direct quadrature (kernel on the scaled profile).
  double worst_shift_abs = 0., worst_shift_rel = 0., worst_bent_rel = 0., worst_straight_rel = 0.;
  int count = 0;
  for (const double k : {0.97, 1.03}) {
    const AtmosphereProfile scaled(scaled_profile_text(kFixture, k));
    for (int deg = 0; deg <= 85; deg += 5) {
      const double theta = radians(deg);
      for (const double h_e : {10000., 30000.}) {
        const CurvedRayIntegrals direct =
            curved_ray_integrals(scaled, kEarthRadius, kObservation, theta, h_e);
        CurvedRayIntegrals approx{};
        CHECK(tables.lookup_scaled(theta, h_e, k, k, approx));
        worst_shift_abs = std::max(worst_shift_abs, std::abs(approx.shift_m - direct.shift_m));
        if (std::abs(direct.shift_m) > 0.01) {
          worst_shift_rel = std::max(worst_shift_rel, std::abs(approx.shift_m - direct.shift_m) /
                                                          std::abs(direct.shift_m));
        }
        if (std::abs(direct.excess_bent_ns) > 1.) {
          worst_bent_rel =
              std::max(worst_bent_rel, std::abs(approx.excess_bent_ns - direct.excess_bent_ns) /
                                           std::abs(direct.excess_bent_ns));
        }
        if (std::abs(direct.excess_straight_ns) > 1e-3) {
          worst_straight_rel = std::max(
              worst_straight_rel, std::abs(approx.excess_straight_ns - direct.excess_straight_ns) /
                                      std::abs(direct.excess_straight_ns));
        }
        ++count;
      }
    }
  }
  // Table interpolation alone (k = 1, same grid): table vs kernel on the unscaled profile.
  double worst_table_shift_rel = 0., worst_table_bent_rel = 0.;
  for (int deg = 0; deg <= 85; deg += 5) {
    for (const double h_e : {10000., 30000.}) {
      const double theta = radians(deg);
      const CurvedRayIntegrals direct =
          curved_ray_integrals(profile, kEarthRadius, kObservation, theta, h_e);
      CurvedRayIntegrals looked{};
      CHECK(tables.lookup(theta, h_e, looked));
      if (std::abs(direct.shift_m) > 0.01) {
        worst_table_shift_rel =
            std::max(worst_table_shift_rel,
                     std::abs(looked.shift_m - direct.shift_m) / std::abs(direct.shift_m));
      }
      worst_table_bent_rel =
          std::max(worst_table_bent_rel, std::abs(looked.excess_bent_ns - direct.excess_bent_ns) /
                                             std::abs(direct.excess_bent_ns));
    }
  }
  std::fprintf(stderr, "info: table interpolation alone (k=1): shift rel %.3e, bent rel %.3e\n",
               worst_table_shift_rel, worst_table_bent_rel);
  std::fprintf(stderr,
               "info: spherical scaled (%d cases, k=0.97,1.03, 0-85 deg, 10/30 km): shift abs "
               "%.3e m, shift rel %.3e, bent rel %.3e, straight rel %.3e\n",
               count, worst_shift_abs, worst_shift_rel, worst_bent_rel, worst_straight_rel);
  CHECK(worst_shift_rel <= 1e-3);
  CHECK(worst_bent_rel <= 1e-3);
  CHECK(worst_straight_rel <= 1e-5);

  // Unequal scales: geometry with k, time with k_group, against the quadrature with that split.
  double worst_group_bent = 0., worst_group_straight = 0.;
  double worst_group_bent_rel = 0.;
  for (const auto [k, k_group] : {std::pair{0.97, 0.98}, std::pair{1.03, 1.02}}) {
    for (int deg = 0; deg <= 85; deg += 5) {
      const double theta = radians(deg);
      for (const double h_e : {10000., 30000.}) {
        const CurvedRayIntegrals ref =
            reference_scaled(profile, kObservation, theta, h_e, k, k_group, 20000);
        CurvedRayIntegrals approx{};
        CHECK(tables.lookup_scaled(theta, h_e, k, k_group, approx));
        worst_group_bent =
            std::max(worst_group_bent, std::abs(approx.excess_bent_ns - ref.excess_bent_ns));
        if (std::abs(ref.excess_bent_ns) > 1.) {
          worst_group_bent_rel =
              std::max(worst_group_bent_rel, std::abs(approx.excess_bent_ns - ref.excess_bent_ns) /
                                                 std::abs(ref.excess_bent_ns));
        }
        if (std::abs(ref.excess_straight_ns) > 1e-3) {
          worst_group_straight = std::max(
              worst_group_straight, std::abs(approx.excess_straight_ns - ref.excess_straight_ns) /
                                        std::abs(ref.excess_straight_ns));
        }
      }
    }
  }
  std::fprintf(stderr,
               "info: spherical scaled, unequal scales: bent abs %.3e ns, bent rel %.3e, "
               "straight rel %.3e\n",
               worst_group_bent, worst_group_bent_rel, worst_group_straight);
  CHECK(worst_group_bent_rel <= 1e-3);
  CHECK(worst_group_straight <= 1e-5);

  // 300 nm at the top of the refraction grid needs k outside the range covered above.
  // Compare with independent direct quadrature at the same phase and group indices.
  {
    constexpr double wavelength_nm = 300.;
    const double k = refractivity_ratio(wavelength_nm, 500.);
    const double k_group = group_refractivity_ratio(wavelength_nm, 500.);
    const double theta = radians(85.);
    constexpr double h_e = 60000.;
    const CurvedRayIntegrals ref =
        reference_scaled(profile, kObservation, theta, h_e, k, k_group, 40000);
    const CurvedRefractionTables accurate_tables(profile, kEarthRadius, kObservation, true, 2);
    CurvedRayIntegrals approx{};
    CHECK(accurate_tables.lookup_scaled(theta, h_e, k, k_group, approx));
    CHECK(std::abs(approx.shift_m - ref.shift_m) / std::abs(ref.shift_m) <= 1e-3);
    CHECK(std::abs(approx.excess_bent_ns - ref.excess_bent_ns) / std::abs(ref.excess_bent_ns) <=
          1e-3);

    // Cover both ends of the supported scale range and the zenith transition.
    for (const auto [lambda, reference] :
         {std::pair{200., 2000.}, std::pair{200., 500.}, std::pair{225., 500.},
          std::pair{250., 500.}, std::pair{300., 500.}, std::pair{350., 500.},
          std::pair{400., 500.}, std::pair{500., 500.}, std::pair{700., 500.},
          std::pair{1000., 500.}, std::pair{1500., 500.}, std::pair{2000., 500.},
          std::pair{2000., 200.}}) {
      const double phase = refractivity_ratio(lambda, reference);
      const double group = group_refractivity_ratio(lambda, reference);
      for (double deg : {75., 77.5, 80., 85., 89.}) {
        const auto truth =
            reference_scaled(profile, kObservation, radians(deg), h_e, phase, group, 40000);
        CurvedRayIntegrals actual{};
        CHECK(accurate_tables.lookup_scaled(radians(deg), h_e, phase, group, actual));
        CHECK(std::abs(actual.shift_m - truth.shift_m) <=
              std::max(0.005, 1e-3 * std::abs(truth.shift_m)));
        CHECK(std::abs(actual.excess_bent_ns - truth.excess_bent_ns) <=
              std::max(0.02, 1e-3 * std::abs(truth.excess_bent_ns)));
      }
    }

    // Approaching either zenith transition from both sides must not jump.
    for (double boundary : {75., 80.}) {
      CurvedRayIntegrals before{}, after{};
      constexpr double eps_deg = 1e-6;
      CHECK(accurate_tables.lookup_scaled(radians(boundary - eps_deg), h_e, k, k_group, before));
      CHECK(accurate_tables.lookup_scaled(radians(boundary + eps_deg), h_e, k, k_group, after));
      CHECK(std::abs(after.shift_m - before.shift_m) < 1e-3);
      CHECK(std::abs(after.excess_bent_ns - before.excess_bent_ns) < 1e-2);
    }

    // The scaled lookup must remain continuous around the physical wavelength scale.
    CurvedRayIntegrals lower{}, upper{};
    constexpr double dk = 1e-6;
    CHECK(accurate_tables.lookup_scaled(theta, h_e, k - dk, k_group - dk, lower));
    CHECK(accurate_tables.lookup_scaled(theta, h_e, k + dk, k_group + dk, upper));
    CHECK(std::abs(upper.shift_m - lower.shift_m) / std::abs(ref.shift_m) < 1e-3);
    CHECK(std::abs(upper.excess_bent_ns - lower.excess_bent_ns) / std::abs(ref.excess_bent_ns) <
          1e-3);
  }

  // Rejections: bad scales and bad angles or heights.
  CurvedRayIntegrals out{};
  CHECK(!tables.lookup_scaled(0.5, 20000., 0., 1., out));
  CHECK(!tables.lookup_scaled(0.5, 20000., 1., kNan, out));
  CHECK(!tables.lookup_scaled(kNan, 20000., 1., 1., out));
  CHECK(!tables.lookup_scaled(0.5, kObservation - 1., 1., 1., out));
}

} // namespace

int main() {
  return run_checks([] {
    const AtmosphereProfile profile = AtmosphereProfile::from_file(kFixture);
    test_ciddor(profile);
    test_plane_scaled(profile);
    test_spherical_scaled(profile);
  });
}
