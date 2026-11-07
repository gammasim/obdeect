#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

namespace dlamb {

// Four uniform deviates in (0, 1) of one photon and stream (gstrickd::uniform4).
using Uniforms = std::array<double, 4>;

inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double two_pi = 2.0 * pi;

// Standard normal quantile for p in (0, 1): Acklam's rational approximation, refined
// by one Halley step (relative error about 1e-15).
[[nodiscard]] inline double normal_quantile(double p) {
  constexpr double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,  -2.759285104469687e+02,
                          1.383577518672690e+02,  -3.066479806614716e+01, 2.506628277459239e+00};
  constexpr double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                          6.680131188771972e+01, -1.328068155288572e+01};
  constexpr double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                          -2.549732539343734e+00, 4.374664141464968e+00,  2.938163982698783e+00};
  constexpr double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                          3.754408661907416e+00};
  constexpr double p_low = 0.02425;
  double x;
  if (p < p_low) {
    const double q = std::sqrt(-2.0 * std::log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  } else if (p <= 1.0 - p_low) {
    const double q = p - 0.5, r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
  } else {
    const double q = std::sqrt(-2.0 * std::log1p(-p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }
  const double e = 0.5 * std::erfc(-x / std::sqrt(2.0)) - p;
  const double step = e * std::sqrt(two_pi) * std::exp(0.5 * x * x);
  return x - step / (1.0 + 0.5 * x * step);
}

// Unnormalised density, linear between nodes and zero outside, sampled by inverting
// its cumulative distribution (port of obdeect LinearDensity). An empty density is a
// point mass at zero.
class LinearDensity {
public:
  LinearDensity() = default;
  LinearDensity(std::vector<double> x, std::vector<double> y) : x_(std::move(x)), y_(std::move(y)) {
    if (x_.size() < 2 || x_.size() != y_.size())
      throw std::invalid_argument("a density needs at least two nodes and one value per node");
    cdf_.assign(x_.size(), 0.0);
    for (std::size_t i = 0; i < x_.size(); ++i) {
      if (!std::isfinite(x_[i]) || !std::isfinite(y_[i]) || y_[i] < 0.0 ||
          (i > 0 && !(x_[i] > x_[i - 1])))
        throw std::invalid_argument("density nodes must increase and values be finite and >= 0");
      if (i > 0)
        cdf_[i] = cdf_[i - 1] + 0.5 * (y_[i] + y_[i - 1]) * (x_[i] - x_[i - 1]);
    }
    if (!(cdf_.back() > 0.0))
      throw std::invalid_argument("a density must have a positive integral");
  }

  [[nodiscard]] bool empty() const { return x_.empty(); }
  [[nodiscard]] double lower() const { return x_.empty() ? 0.0 : x_.front(); }
  [[nodiscard]] double upper() const { return x_.empty() ? 0.0 : x_.back(); }
  [[nodiscard]] double integral() const { return cdf_.empty() ? 0.0 : cdf_.back(); }
  [[nodiscard]] double maximum() const {
    return y_.empty() ? 0.0 : *std::max_element(y_.begin(), y_.end());
  }

  [[nodiscard]] double at(double x) const {
    if (x_.empty() || x < x_.front() || x > x_.back())
      return 0.0;
    const auto i = segment(x_, x);
    return y_[i] + (y_[i + 1] - y_[i]) * (x - x_[i]) / (x_[i + 1] - x_[i]);
  }

  // Integral of the density below x (unnormalised).
  [[nodiscard]] double cdf(double x) const {
    if (x_.empty() || x <= x_.front())
      return 0.0;
    if (x >= x_.back())
      return cdf_.back();
    const auto i = segment(x_, x);
    const double s = x - x_[i], k = (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
    return cdf_[i] + y_[i] * s + 0.5 * k * s * s;
  }

  // The value x with cdf(x) = u * integral(), for u in (0, 1).
  [[nodiscard]] double sample(double u) const {
    if (x_.empty())
      return 0.0;
    const double area = u * cdf_.back();
    const auto i = segment(cdf_, area);
    // Solve y_i s + k s^2 / 2 = a for s in the segment, stable for any slope k.
    const double a = area - cdf_[i], k = (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
    const double root = std::sqrt(std::max(0.0, y_[i] * y_[i] + 2.0 * k * a));
    const double s = y_[i] + root > 0.0 ? 2.0 * a / (y_[i] + root) : 0.0;
    return std::min(x_[i] + s, x_[i + 1]);
  }

private:
  // Segment [i, i + 1] of the nodes containing `value` of an increasing column
  // (x_ or cdf_), clamped to the table.
  [[nodiscard]] static std::size_t segment(const std::vector<double> &column, double value) {
    const auto upper = std::upper_bound(column.begin(), column.end(), value) - column.begin();
    return static_cast<std::size_t>(
        std::clamp<std::ptrdiff_t>(upper - 1, 0, static_cast<std::ptrdiff_t>(column.size()) - 2));
  }

  std::vector<double> x_, y_, cdf_;
};

// ---------------------------------------------------------------------------
// One-dimensional time and wavelength shapes. Each struct samples from the uniform
// deviates of one stream: sample(u) uses u[0] (the composite pulse uses u[0..2]).
// Moments are the analytic values of the shape, for tests and documentation.
// ---------------------------------------------------------------------------

// Delta function at `value`.
struct Fixed {
  double value{};
  [[nodiscard]] double sample(const Uniforms &) const { return value; }
  [[nodiscard]] double mean() const { return value; }
  [[nodiscard]] double variance() const { return 0.0; }
};

// Flat (top-hat) distribution on [low, high].
struct UniformRange {
  double low{};
  double high{};
  [[nodiscard]] double sample(const Uniforms &u) const { return low + (high - low) * u[0]; }
  [[nodiscard]] double mean() const { return 0.5 * (low + high); }
  [[nodiscard]] double variance() const { return (high - low) * (high - low) / 12.0; }
};

// Normal distribution with mean and standard deviation sigma.
struct GaussianShape {
  double mean_value{};
  double sigma{};
  [[nodiscard]] double sample(const Uniforms &u) const {
    return mean_value + sigma * normal_quantile(u[0]);
  }
  [[nodiscard]] double mean() const { return mean_value; }
  [[nodiscard]] double variance() const { return sigma * sigma; }
};

// Exponential decay with decay time tau, mean `mean_value` (leading edge at mean - tau),
// as LightEmission "Exp:" (t = mean + tau (E - 1), E ~ Exp(1)).
struct Exponential {
  double mean_value{};
  double tau{};
  [[nodiscard]] double sample(const Uniforms &u) const {
    return mean_value + tau * (-std::log1p(-u[0]) - 1.0);
  }
  [[nodiscard]] double mean() const { return mean_value; }
  [[nodiscard]] double variance() const { return tau * tau; }
};

// Sum of a flat pulse of full width flat_width, a Gaussian of r.m.s. gauss_rms and an
// exponential of decay time decay (LightEmission "Complex:"), all centred on delay.
struct ComplexPulse {
  double flat_width{};
  double gauss_rms{};
  double decay{};
  double delay{};
  [[nodiscard]] double sample(const Uniforms &u) const {
    return delay + flat_width * (u[0] - 0.5) + gauss_rms * normal_quantile(u[1]) +
           decay * (-std::log1p(-u[2]) - 1.0);
  }
  [[nodiscard]] double mean() const { return delay; }
  [[nodiscard]] double variance() const {
    return flat_width * flat_width / 12.0 + gauss_rms * gauss_rms + decay * decay;
  }
};

// Density proportional to 1/lambda^2 on [low, high] (LightEmission "Cherenkov:").
struct InverseSquare {
  double low{};
  double high{};
  [[nodiscard]] double sample(const Uniforms &u) const {
    const double a = 1.0 / low, b = 1.0 / high;
    return 1.0 / (a - u[0] * (a - b));
  }
  [[nodiscard]] double mean() const { return std::log(high / low) / (1.0 / low - 1.0 / high); }
  [[nodiscard]] double variance() const {
    const double norm = 1.0 / low - 1.0 / high;
    const double second = (high - low) / norm;
    return second - mean() * mean();
  }
};

// Tabulated density over x (wavelength [nm] or time [ns]).
struct Tabulated {
  LinearDensity density;
  [[nodiscard]] double sample(const Uniforms &u) const { return density.sample(u[0]); }
};

using Spectrum = std::variant<Fixed, UniformRange, InverseSquare, Tabulated>;
using Pulse =
    std::variant<Fixed, UniformRange, GaussianShape, Exponential, ComplexPulse, Tabulated>;

// ---------------------------------------------------------------------------
// Angular profiles. The argument u = 1 - cos(theta) is measured from the pointing
// axis, so a density in u is a density per solid angle (dOmega = 2 pi du). Intensity
// is the relative emission per solid angle; cdf(u) and quantile(p) refer to the
// emitted fraction.
// ---------------------------------------------------------------------------

// Isotropic emission into 4 pi.
struct Isotropic {
  [[nodiscard]] double intensity(double u) const { return u >= 0.0 && u <= 2.0 ? 1.0 : 0.0; }
  [[nodiscard]] double maximum() const { return 1.0; }
  [[nodiscard]] double integral() const { return 2.0; }
  [[nodiscard]] double cdf(double u) const { return std::clamp(0.5 * u, 0.0, 1.0); }
  [[nodiscard]] double quantile(double p) const { return 2.0 * p; }
};

// All photons at one angle with u0 = 1 - cos(theta) (LightEmission "Cone:", "Parallel"
// for u0 = 0). It has no density, so it cannot be used with a disc target.
struct FixedAngle {
  double u0{};
  // No density: the intensity members are zero and the disc target rejects this profile.
  [[nodiscard]] double intensity(double) const { return 0.0; }
  [[nodiscard]] double maximum() const { return 0.0; }
  [[nodiscard]] double integral() const { return 0.0; }
  [[nodiscard]] double cdf(double u) const { return u >= u0 ? 1.0 : 0.0; }
  [[nodiscard]] double quantile(double) const { return u0; }
};

// Uniform intensity within u <= u0 (LightEmission "FilledCone:").
struct FilledCone {
  double u0{};
  [[nodiscard]] double intensity(double u) const { return u >= 0.0 && u <= u0 ? 1.0 : 0.0; }
  [[nodiscard]] double maximum() const { return 1.0; }
  [[nodiscard]] double integral() const { return u0; }
  [[nodiscard]] double cdf(double u) const { return std::clamp(u / u0, 0.0, 1.0); }
  [[nodiscard]] double quantile(double p) const { return p * u0; }
};

// Lambertian emission, intensity proportional to cos(theta) over the hemisphere.
// In u: (1 - u) on [0, 1].
struct Lambertian {
  [[nodiscard]] double intensity(double u) const { return u >= 0.0 && u <= 1.0 ? 1.0 - u : 0.0; }
  [[nodiscard]] double maximum() const { return 1.0; }
  [[nodiscard]] double integral() const { return 0.5; }
  [[nodiscard]] double cdf(double u) const {
    const double v = std::clamp(u, 0.0, 1.0);
    return 2.0 * v - v * v;
  }
  [[nodiscard]] double quantile(double p) const { return 1.0 - std::sqrt(1.0 - p); }
};

// Tabulated intensity over u (Gaussian, Rayleigh and user tables).
struct TabulatedAngular {
  LinearDensity density;
  [[nodiscard]] double intensity(double u) const { return density.at(u); }
  [[nodiscard]] double maximum() const { return density.maximum(); }
  [[nodiscard]] double integral() const { return density.integral(); }
  [[nodiscard]] double cdf(double u) const {
    return density.empty() ? 0.0 : density.cdf(u) / density.integral();
  }
  [[nodiscard]] double quantile(double p) const { return density.sample(p); }
};

using Angular = std::variant<Isotropic, FixedAngle, FilledCone, Lambertian, TabulatedAngular>;

[[nodiscard]] inline double intensity_maximum(const Angular &a) {
  return std::visit([](const auto &x) { return x.maximum(); }, a);
}

// The u with the emitted fraction p below it.
[[nodiscard]] inline double angle_quantile(const Angular &a, double p) {
  return std::visit([p](const auto &x) { return x.quantile(p); }, a);
}

// Integral of the intensity over u (the solid angle integral is two pi times this).
[[nodiscard]] inline double integral(const Angular &a) {
  return std::visit([](const auto &x) { return x.integral(); }, a);
}

// Gaussian intensity exp(-theta^2 / (2 sigma^2)) (LightEmission "Gauss-rms:"), tabulated
// on a grid of 0.05 sigma up to 5 sigma (or 180 deg).
[[nodiscard]] inline LinearDensity gaussian_angular_table(double sigma_rad) {
  if (!(sigma_rad > 0.0) || !std::isfinite(sigma_rad))
    throw std::invalid_argument("a Gaussian angular distribution needs a positive width");
  const double theta_max = std::min(5.0 * sigma_rad, pi);
  const double step = 0.05 * sigma_rad;
  std::vector<double> u, value;
  for (double theta = 0.0; theta < theta_max; theta += step) {
    u.push_back(1.0 - std::cos(theta));
    value.push_back(std::exp(-0.5 * theta * theta / (sigma_rad * sigma_rad)));
  }
  u.push_back(1.0 - std::cos(theta_max));
  value.push_back(std::exp(-0.5 * theta_max * theta_max / (sigma_rad * sigma_rad)));
  return {u, value};
}

// Rayleigh-scattering intensity 0.5 (1 + cos^2 theta) (LightEmission "Rayleigh"),
// tabulated in 1 degree steps.
[[nodiscard]] inline LinearDensity rayleigh_angular_table() {
  std::vector<double> u, value;
  for (int deg = 0; deg <= 180; ++deg) {
    const double c = std::cos(deg * pi / 180.0);
    u.push_back(1.0 - c);
    value.push_back(0.5 * (1.0 + c * c));
  }
  return {u, value};
}

} // namespace dlamb
