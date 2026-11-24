#include "check.hpp"

#include "obdeect/sources/density.hpp"
#include "obdeect/sources/streams.hpp"

#include "obdeect/stream/random.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace obdeect::sources;

namespace {

constexpr std::size_t n = 200000;

// Uniform deviates of photon i for a test stream (keys as in obdeect/sources/source.hpp).
Uniforms draw(std::uint64_t i, Stream stream = Stream::pulse) {
  return obdeect::stream::uniform4(
      {12345, stage_sources, 0, 0, i, static_cast<std::uint32_t>(stream)});
}

// Sample mean and variance (unbiased) of f over n draws.
template <class F> std::pair<double, double> moments(F f) {
  double s1 = 0.0, s2 = 0.0;
  for (std::uint64_t i = 0; i < n; ++i) {
    const double x = f(i);
    s1 += x;
    s2 += x * x;
  }
  const double mean = s1 / n;
  return {mean, (s2 - n * mean * mean) / (n - 1)};
}

// Four-sigma checks of the sample mean and variance against analytic values. The
// variance of a sample variance is (kurtosis - 1) variance^2 / n; the kurtosis of the
// shapes tested here is below 10, which bounds the sigma of the check. The shape itself is
// checked by check_cdf (the moments alone do not test it).
template <class F> void check_moments(const char *name, double mean, double variance, F f) {
  const auto [m, v] = moments(f);
  const double mean_sigma = std::sqrt(variance / n);
  const double variance_sigma = variance * std::sqrt(10.0 / n);
  std::fprintf(stderr, "%-22s mean %.6g (expected %.6g), variance %.6g (expected %.6g)\n", name, m,
               mean, v, variance);
  CHECK(std::abs(m - mean) <= 4.0 * mean_sigma + 1e-12);
  CHECK(std::abs(v - variance) <= 4.0 * variance_sigma + 1e-12);
}

// Kolmogorov-Smirnov test of n draws against an analytic cdf: sqrt(n) D_n < 2.0 (p about 7e-4 for
// a correct sampler). The draws are deterministic, so the check is a fixed test of the shape.
template <class Sample, class Cdf> void check_cdf(const char *name, Sample sample, Cdf cdf) {
  std::vector<double> x;
  x.reserve(n);
  for (std::uint64_t i = 0; i < n; ++i)
    x.push_back(sample(i));
  std::sort(x.begin(), x.end());
  double d = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    const double f = cdf(x[k]);
    d = std::max(d,
                 std::max((static_cast<double>(k) + 1.0) / n - f, f - static_cast<double>(k) / n));
  }
  const double statistic = std::sqrt(static_cast<double>(n)) * d;
  std::fprintf(stderr, "%-22s KS sqrt(n) D %.4f\n", name, statistic);
  CHECK(statistic < 2.0);
}

// Mean and variance of a density over u by midpoint integration (for tables).
std::pair<double, double> numeric_moments_u(const Angular &a) {
  const double total = integral(a);
  double m1 = 0.0, m2 = 0.0;
  const int steps = 400000;
  for (int k = 0; k < steps; ++k) {
    const double u = (k + 0.5) * 2.0 / steps;
    const double w = std::visit([u](const auto &x) { return x.intensity(u); }, a) * (2.0 / steps);
    m1 += u * w;
    m2 += u * u * w;
  }
  m1 /= total;
  m2 /= total;
  return {m1, m2 - m1 * m1};
}

// Integral of f over [0, b] by Simpson's rule with 2e5 panels (reference quadrature).
template <class F> double simpson(F f, double b) {
  const int panels = 200000;
  const double h = b / panels;
  double sum = f(0.0) + f(b);
  for (int k = 1; k < panels; ++k)
    sum += (k % 2 == 1 ? 4.0 : 2.0) * f(k * h);
  return sum * h / 3.0;
}

void quantile_reference() {
  CHECK_NEAR(normal_quantile(0.5), 0.0, 1e-15);
  CHECK_NEAR(normal_quantile(0.975), 1.959963984540054, 1e-12);
  CHECK_NEAR(normal_quantile(1e-6), -4.753424308822899, 1e-10);
  CHECK_NEAR(normal_quantile(1.0 - 1e-6), 4.753424308822899, 1e-10);
}

void tabulated_density() {
  // The same density as obdeect test_kernels.cpp (densities): 2x on [0, 1], 2 on [1, 3].
  const LinearDensity d({0.0, 1.0, 3.0}, {0.0, 2.0, 2.0});
  CHECK_NEAR(d.integral(), 5.0, 1e-15);
  CHECK_NEAR(d.maximum(), 2.0, 0.0);
  CHECK_NEAR(d.at(0.5), 1.0, 1e-15);
  CHECK_NEAR(d.at(3.5), 0.0, 0.0);
  CHECK_NEAR(d.sample(0.1), std::sqrt(0.5), 1e-15);
  CHECK_NEAR(d.sample(0.6), 2.0, 1e-15);
  CHECK_NEAR(d.sample(1.0), 3.0, 1e-15);
  CHECK_NEAR(d.cdf(1.0), 1.0, 1e-15);
  CHECK_NEAR(d.cdf(2.0), 3.0, 1e-15);
  for (double u : {0.05, 0.3, 0.77, 0.999})
    CHECK_NEAR(d.cdf(d.sample(u)) / d.integral(), u, 1e-12);
  CHECK_NEAR(LinearDensity().sample(0.3), 0.0, 0.0); // Empty: point mass at zero.
  bool rejected = false;
  try {
    const LinearDensity bad({1.0, 0.0}, {1.0, 1.0});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
  // Tabulated density sampling reproduces its own mean: triangle on [0, 2] peaking at 1.
  const LinearDensity triangle({0.0, 1.0, 2.0}, {0.0, 1.0, 0.0});
  check_moments("triangle (table)", 1.0, 1.0 / 6.0,
                [&](std::uint64_t i) { return Tabulated{triangle}.sample(draw(i)); });
  // Triangle cdf: x^2 / 2 on [0, 1], 1 - (2 - x)^2 / 2 on [1, 2].
  check_cdf(
      "triangle (table) cdf", [&](std::uint64_t i) { return Tabulated{triangle}.sample(draw(i)); },
      [](double x) { return x <= 1.0 ? 0.5 * x * x : 1.0 - 0.5 * (2.0 - x) * (2.0 - x); });
}

void one_dimensional_shapes() {
  check_moments("fixed", 2.5, 0.0, [](std::uint64_t i) { return Fixed{2.5}.sample(draw(i)); });
  const UniformRange flat{-1.25, 1.25};
  check_moments("uniform", 0.0, 2.5 * 2.5 / 12.0,
                [&](std::uint64_t i) { return flat.sample(draw(i)); });
  check_cdf(
      "uniform cdf", [&](std::uint64_t i) { return flat.sample(draw(i)); },
      [&](double x) { return std::clamp((x - flat.low) / (flat.high - flat.low), 0.0, 1.0); });
  const GaussianShape gauss{3.0, 0.7};
  check_moments("gaussian", 3.0, 0.49, [&](std::uint64_t i) { return gauss.sample(draw(i)); });
  check_cdf(
      "gaussian cdf", [&](std::uint64_t i) { return gauss.sample(draw(i)); },
      [&](double x) { return 0.5 * std::erfc(-(x - 3.0) / (0.7 * std::sqrt(2.0))); });
  const Exponential exp_shape{1.5, 2.0};
  check_moments("exponential", 1.5, 4.0,
                [&](std::uint64_t i) { return exp_shape.sample(draw(i)); });
  // Leading edge at mean - tau: x - mean = tau (E - 1) with E exponential of unit mean.
  check_cdf(
      "exponential cdf", [&](std::uint64_t i) { return exp_shape.sample(draw(i)); },
      [&](double x) {
        const double e = (x - exp_shape.mean_value) / exp_shape.tau + 1.0;
        return e <= 0.0 ? 0.0 : 1.0 - std::exp(-e);
      });
  const ComplexPulse complex{2.0, 1.0, 1.0, 0.5};
  check_moments("complex pulse", 0.5, complex.variance(),
                [&](std::uint64_t i) { return complex.sample(draw(i)); });
  // Cherenkov spectrum: density 1/lambda^2 on [300, 600] nm (option "Cherenkov:").
  const InverseSquare cherenkov{300.0, 600.0};
  const double mean = std::log(2.0) / (1.0 / 300.0 - 1.0 / 600.0);
  CHECK_NEAR(cherenkov.mean(), mean, 1e-9);
  check_moments("inverse square", mean, cherenkov.variance(),
                [&](std::uint64_t i) { return cherenkov.sample(draw(i, Stream::wavelength)); });
  // cdf of the density 1/lambda^2: (1/low - 1/lambda) / (1/low - 1/high).
  check_cdf(
      "inverse square cdf",
      [&](std::uint64_t i) { return cherenkov.sample(draw(i, Stream::wavelength)); },
      [](double x) {
        return std::clamp((1.0 / 300.0 - 1.0 / x) / (1.0 / 300.0 - 1.0 / 600.0), 0.0, 1.0);
      });
  // Every draw lies inside the range.
  for (std::uint64_t i = 0; i < 20000; ++i) {
    const double wl = cherenkov.sample(draw(i, Stream::wavelength));
    CHECK(wl >= 300.0 && wl <= 600.0);
  }
}

void angular_profiles() {
  // Isotropic: uniform in u on [0, 2]; mean u = 1.
  const Angular iso = Isotropic{};
  check_moments("isotropic u", 1.0, 1.0 / 3.0,
                [&](std::uint64_t i) { return angle_quantile(iso, draw(i)[0]); });
  check_cdf(
      "isotropic u cdf", [&](std::uint64_t i) { return angle_quantile(iso, draw(i)[0]); },
      [](double u) { return std::clamp(u / 2.0, 0.0, 1.0); });
  // Lambertian (intensity cos theta): u has density 1 - u on [0, 1]; mean 1/3.
  const Angular lambert = Lambertian{};
  check_moments("lambertian u", 1.0 / 3.0, 1.0 / 18.0,
                [&](std::uint64_t i) { return angle_quantile(lambert, draw(i)[0]); });
  check_cdf(
      "lambertian u cdf", [&](std::uint64_t i) { return angle_quantile(lambert, draw(i)[0]); },
      [](double u) {
        const double v = std::clamp(u, 0.0, 1.0);
        return 2.0 * v - v * v;
      });
  for (double p : {0.1, 0.5, 0.9})
    CHECK_NEAR(std::visit([p](const auto &a) { return a.cdf(a.quantile(p)); }, lambert), p, 1e-12);
  // Filled cone of half angle 10 deg: uniform on [0, u0].
  const double u0 = 1.0 - std::cos(10.0 * pi / 180.0);
  const Angular cone = FilledCone{u0};
  check_moments("filled cone u", 0.5 * u0, u0 * u0 / 12.0,
                [&](std::uint64_t i) { return angle_quantile(cone, draw(i)[0]); });
  check_cdf(
      "filled cone u cdf", [&](std::uint64_t i) { return angle_quantile(cone, draw(i)[0]); },
      [u0](double u) { return std::clamp(u / u0, 0.0, 1.0); });
  // Gaussian (table) and Rayleigh (table): moments against midpoint integration of the table.
  for (const auto &table : {gaussian_angular_table(10.0 * pi / 180.0), rayleigh_angular_table()}) {
    const Angular a = TabulatedAngular{table};
    const auto [mean_u, var_u] = numeric_moments_u(a);
    check_moments("tabulated u", mean_u, var_u,
                  [&](std::uint64_t i) { return angle_quantile(a, draw(i)[0]); });
  }
  // Gaussian table against its analytic cdf: the light of exp(-theta^2 / 2 sigma^2) inside the
  // cone of angle theta, weighted by the solid angle sin(theta) (truncated at 5 sigma, as the
  // table).
  const double sigma = 10.0 * pi / 180.0;
  const double theta_max = std::min(5.0 * sigma, pi);
  const auto light = [sigma](double theta) {
    return std::exp(-theta * theta / (2.0 * sigma * sigma)) * std::sin(theta);
  };
  // Cumulative trapezoid of the light on 200000 panels (error ~1e-10), interpolated in theta.
  const int panels = 200000;
  const double step = theta_max / panels;
  std::vector<double> cumulative(panels + 1, 0.0);
  for (int k = 0; k < panels; ++k)
    cumulative[k + 1] = cumulative[k] + 0.5 * (light(k * step) + light((k + 1) * step)) * step;
  const double total = cumulative.back();
  const auto analytic_cdf = [&](double theta) {
    const double t = std::min(theta, theta_max) / step;
    const auto k = std::min(static_cast<int>(t), panels - 1);
    const double frac = t - k;
    return (cumulative[k] + frac * (cumulative[k + 1] - cumulative[k])) / total;
  };
  const Angular gauss_table = TabulatedAngular{gaussian_angular_table(sigma)};
  check_cdf(
      "gaussian table cdf",
      [&](std::uint64_t i) { return angle_quantile(gauss_table, draw(i)[0]); },
      [&](double u) { return analytic_cdf(std::acos(std::clamp(1.0 - u, -1.0, 1.0))); });
  // Table values: the density at the table's own u equals exp(-theta^2 / 2 sigma^2).
  const auto g = gaussian_angular_table(sigma);
  for (double theta_deg : {2.0, 5.0, 10.0, 20.0}) {
    const double theta = theta_deg * pi / 180.0;
    CHECK_NEAR(g.at(1.0 - std::cos(theta)), std::exp(-theta * theta / (2.0 * sigma * sigma)), 1e-3);
  }
  CHECK_NEAR(g.at(0.0), 1.0, 1e-15);
  // Widths whose square underflows (or are not finite) are rejected, not tabulated as NaN.
  for (double bad_sigma : {1e-200, 0.0, -1.0, std::numeric_limits<double>::infinity()}) {
    bool rejected_width = false;
    try {
      (void)gaussian_angular_table(bad_sigma);
    } catch (const std::invalid_argument &) {
      rejected_width = true;
    }
    CHECK(rejected_width);
  }
}

} // namespace

int main() {
  return run_checks([] {
    quantile_reference();
    tabulated_density();
    one_dimensional_shapes();
    angular_profiles();
  });
}
