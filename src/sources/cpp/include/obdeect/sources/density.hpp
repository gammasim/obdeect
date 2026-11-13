#pragma once

// Densities, pulse and spectrum shapes and angular profiles (docs/spec/distributions.md).

#include "obdeect/sources/vec3.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

namespace obdeect::sources {

using Uniforms = std::array<double, 4>;

// FWHM to sigma of a Gaussian.
inline constexpr double fwhm_per_sigma = 2.35482004503095;

// Inverse of the standard normal distribution function (relative error ~1e-15).
inline double normal_quantile(double p) {
  if (!(p > 0.0))
    return p == 0.0 ? -std::numeric_limits<double>::infinity()
                    : std::numeric_limits<double>::quiet_NaN();
  if (!(p < 1.0))
    return p == 1.0 ? std::numeric_limits<double>::infinity()
                    : std::numeric_limits<double>::quiet_NaN();
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
    const double q = p - 0.5;
    const double r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
  } else {
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }
  // Two Halley steps bring the approximation to machine precision.
  for (int step = 0; step < 2; ++step) {
    const double e = 0.5 * std::erfc(-x / std::sqrt(2.0)) - p;
    const double u = e * std::sqrt(two_pi) * std::exp(0.5 * x * x);
    x -= u / (1.0 + 0.5 * x * u);
  }
  return x;
}

// Unnormalised density, piecewise linear between nodes, zero outside [x_0, x_n].
class LinearDensity {
public:
  LinearDensity() = default;
  LinearDensity(std::vector<double> x, std::vector<double> y) : x_(std::move(x)), y_(std::move(y)) {
    if (x_.size() < 2 || x_.size() != y_.size())
      throw std::invalid_argument("a tabulated density needs at least two nodes of equal count");
    for (std::size_t i = 0; i < x_.size(); ++i) {
      if (!std::isfinite(x_[i]) || !std::isfinite(y_[i]))
        throw std::invalid_argument("a tabulated density needs finite values");
      if (y_[i] < 0.0)
        throw std::invalid_argument("a tabulated density needs non-negative values");
      if (i > 0 && !(x_[i] > x_[i - 1]))
        throw std::invalid_argument("the abscissae of a tabulated density must increase");
    }
    cumulative_.assign(x_.size(), 0.0);
    for (std::size_t i = 1; i < x_.size(); ++i) {
      cumulative_[i] = cumulative_[i - 1] + 0.5 * (y_[i] + y_[i - 1]) * (x_[i] - x_[i - 1]);
      maximum_ = std::max(maximum_, y_[i]);
    }
    maximum_ = std::max(maximum_, y_[0]);
    if (!std::isfinite(cumulative_.back()) || !(cumulative_.back() > 0.0))
      throw std::invalid_argument("a tabulated density needs a positive, finite integral");
  }

  [[nodiscard]] bool empty() const { return x_.empty(); }
  [[nodiscard]] double integral() const { return cumulative_.empty() ? 0.0 : cumulative_.back(); }
  [[nodiscard]] double maximum() const { return maximum_; }
  [[nodiscard]] double lower() const { return x_.empty() ? 0.0 : x_.front(); }
  [[nodiscard]] double upper() const { return x_.empty() ? 0.0 : x_.back(); }
  [[nodiscard]] const std::vector<double> &nodes_x() const { return x_; }
  [[nodiscard]] const std::vector<double> &nodes_y() const { return y_; }

  [[nodiscard]] double at(double x) const {
    if (x_.empty() || x < x_.front() || x > x_.back())
      return 0.0;
    const std::size_t i = segment_of(x);
    const double s = x - x_[i];
    return y_[i] + (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]) * s;
  }

  // Integral from the first node to x.
  [[nodiscard]] double cdf(double x) const {
    if (x_.empty() || x <= x_.front())
      return 0.0;
    if (x >= x_.back())
      return cumulative_.back();
    const std::size_t i = segment_of(x);
    const double s = x - x_[i];
    const double slope = (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
    return cumulative_[i] + y_[i] * s + 0.5 * slope * s * s;
  }

  // The x with cdf(x) = p * integral, for p in (0, 1).
  [[nodiscard]] double sample(double p) const {
    if (x_.empty())
      return 0.0;
    const double target = std::clamp(p, 0.0, 1.0) * cumulative_.back();
    // Segment i with cumulative_[i] < target <= cumulative_[i + 1].
    auto it = std::lower_bound(cumulative_.begin(), cumulative_.end(), target);
    std::size_t i =
        it == cumulative_.begin() ? 0 : static_cast<std::size_t>(it - cumulative_.begin()) - 1;
    i = std::min(i, x_.size() - 2);
    const double width = x_[i + 1] - x_[i];
    const double slope = (y_[i + 1] - y_[i]) / width;
    const double r = std::max(0.0, target - cumulative_[i]);
    const double root = std::sqrt(std::max(0.0, y_[i] * y_[i] + 2.0 * slope * r));
    const double denominator = y_[i] + root;
    const double s = denominator > 0.0 ? 2.0 * r / denominator : 0.0;
    return x_[i] + std::clamp(s, 0.0, width);
  }

private:
  [[nodiscard]] std::size_t segment_of(double x) const {
    auto it = std::upper_bound(x_.begin(), x_.end(), x);
    const std::size_t i = it == x_.begin() ? 0 : static_cast<std::size_t>(it - x_.begin()) - 1;
    return std::min(i, x_.size() - 2);
  }

  std::vector<double> x_, y_, cumulative_;
  double maximum_{};
};

// One-dimensional shapes. sample(u) uses u[0] (ComplexPulse u[0..2]).

struct Fixed {
  double value{};
  [[nodiscard]] double sample(const Uniforms &) const { return value; }
  [[nodiscard]] double mean() const { return value; }
  [[nodiscard]] double variance() const { return 0.0; }
};

struct UniformRange {
  double low{}, high{};
  [[nodiscard]] double sample(const Uniforms &u) const { return low + (high - low) * u[0]; }
  [[nodiscard]] double mean() const { return 0.5 * (low + high); }
  [[nodiscard]] double variance() const { return (high - low) * (high - low) / 12.0; }
};

struct GaussianShape {
  double mean_value{}, sigma{};
  [[nodiscard]] double sample(const Uniforms &u) const {
    return mean_value + sigma * normal_quantile(u[0]);
  }
  [[nodiscard]] double mean() const { return mean_value; }
  [[nodiscard]] double variance() const { return sigma * sigma; }
};

// Leading edge at mean_value - tau; mean mean_value.
struct Exponential {
  double mean_value{}, tau{};
  [[nodiscard]] double sample(const Uniforms &u) const {
    return mean_value + tau * (-std::log(1.0 - u[0]) - 1.0);
  }
  [[nodiscard]] double mean() const { return mean_value; }
  [[nodiscard]] double variance() const { return tau * tau; }
};

// Flat top, Gaussian smearing and exponential decay, summed.
struct ComplexPulse {
  double flat_width{}, gauss_rms{}, decay{}, delay{};
  [[nodiscard]] double sample(const Uniforms &u) const {
    return delay + flat_width * (u[0] - 0.5) + gauss_rms * normal_quantile(u[1]) +
           decay * (-std::log(1.0 - u[2]) - 1.0);
  }
  [[nodiscard]] double mean() const { return delay; }
  [[nodiscard]] double variance() const {
    return flat_width * flat_width / 12.0 + gauss_rms * gauss_rms + decay * decay;
  }
};

// Density proportional to 1 / lambda^2 on [low, high].
struct InverseSquare {
  double low{}, high{};
  [[nodiscard]] double sample(const Uniforms &u) const {
    const double a = 1.0 / low, b = 1.0 / high;
    return 1.0 / (a - u[0] * (a - b));
  }
  [[nodiscard]] double mean() const { return std::log(high / low) / (1.0 / low - 1.0 / high); }
  [[nodiscard]] double variance() const {
    const double m = mean();
    return (high - low) / (1.0 / low - 1.0 / high) - m * m;
  }
};

struct Tabulated {
  LinearDensity density;
  [[nodiscard]] double sample(const Uniforms &u) const { return density.sample(u[0]); }
};

using Spectrum = std::variant<Fixed, UniformRange, InverseSquare, Tabulated>;
using Pulse =
    std::variant<Fixed, UniformRange, GaussianShape, Exponential, ComplexPulse, Tabulated>;

template <class Shape> double sample(const Shape &shape, const Uniforms &u) {
  return std::visit([&u](const auto &s) { return s.sample(u); }, shape);
}

// Angular profiles in u = 1 - cos(theta), u in [0, 2].

struct Isotropic {
  [[nodiscard]] double intensity(double u) const { return u >= 0.0 && u <= 2.0 ? 1.0 : 0.0; }
  [[nodiscard]] double maximum() const { return 1.0; }
  [[nodiscard]] double integral() const { return 2.0; }
  [[nodiscard]] double cdf(double u) const { return std::clamp(u / 2.0, 0.0, 1.0); }
  [[nodiscard]] double quantile(double p) const { return 2.0 * p; }
};

// All photons at exactly one angle: no density.
struct FixedAngle {
  double u0{};
  [[nodiscard]] double intensity(double) const { return 0.0; }
  [[nodiscard]] double maximum() const { return 0.0; }
  [[nodiscard]] double integral() const { return 0.0; }
  [[nodiscard]] double cdf(double u) const { return u >= u0 ? 1.0 : 0.0; }
  [[nodiscard]] double quantile(double) const { return u0; }
};

struct FilledCone {
  double u0{};
  [[nodiscard]] double intensity(double u) const { return u >= 0.0 && u <= u0 ? 1.0 : 0.0; }
  [[nodiscard]] double maximum() const { return 1.0; }
  [[nodiscard]] double integral() const { return u0; }
  [[nodiscard]] double cdf(double u) const { return std::clamp(u / u0, 0.0, 1.0); }
  [[nodiscard]] double quantile(double p) const { return p * u0; }
};

// Intensity proportional to cos(theta) = 1 - u.
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

struct TabulatedAngular {
  LinearDensity density;
  [[nodiscard]] double intensity(double u) const { return density.at(u); }
  [[nodiscard]] double maximum() const { return density.maximum(); }
  [[nodiscard]] double integral() const { return density.integral(); }
  [[nodiscard]] double cdf(double u) const { return density.cdf(u) / density.integral(); }
  [[nodiscard]] double quantile(double p) const { return density.sample(p); }
};

using Angular = std::variant<Isotropic, FixedAngle, FilledCone, Lambertian, TabulatedAngular>;

inline double intensity_maximum(const Angular &a) {
  return std::visit([](const auto &x) { return x.maximum(); }, a);
}
inline double integral(const Angular &a) {
  return std::visit([](const auto &x) { return x.integral(); }, a);
}
inline double angle_quantile(const Angular &a, double p) {
  return std::visit([p](const auto &x) { return x.quantile(p); }, a);
}

// u = 1 - cos(theta), evaluated without cancellation for small angles.
inline double u_of_angle(double theta) {
  const double s = std::sin(0.5 * theta);
  return 2.0 * s * s;
}

// Gaussian angular intensity exp(-theta^2 / (2 sigma^2)); sigma is the 2-D r.m.s. in angle.
inline LinearDensity gaussian_angular_table(double sigma_rad) {
  if (!std::isfinite(sigma_rad) || sigma_rad < 1e-12)
    throw std::invalid_argument("the Gaussian angular width must be finite and at least 1e-12 rad");
  const double theta_max = std::min(5.0 * sigma_rad, pi);
  const double step = 0.05 * sigma_rad;
  std::vector<double> u, y;
  for (double theta = 0.0; theta < theta_max; theta += step) {
    const double ui = u_of_angle(theta);
    if (!u.empty() && !(ui > u.back()))
      continue;
    u.push_back(ui);
    y.push_back(std::exp(-theta * theta / (2.0 * sigma_rad * sigma_rad)));
  }
  const double u_last = u_of_angle(theta_max);
  const double y_last = std::exp(-theta_max * theta_max / (2.0 * sigma_rad * sigma_rad));
  if (!u.empty() && !(u_last > u.back())) {
    u.back() = u_last;
    y.back() = y_last;
  } else {
    u.push_back(u_last);
    y.push_back(y_last);
  }
  return LinearDensity(std::move(u), std::move(y));
}

// Rayleigh scattering intensity 0.5 (1 + cos^2 theta), nodes every degree.
inline LinearDensity rayleigh_angular_table() {
  std::vector<double> u, y;
  for (int degrees = 0; degrees <= 180; ++degrees) {
    const double theta = degrees * pi / 180.0;
    u.push_back(1.0 - std::cos(theta));
    y.push_back(0.5 * (1.0 + std::cos(theta) * std::cos(theta)));
  }
  u.front() = 0.0;
  u.back() = 2.0;
  return LinearDensity(std::move(u), std::move(y));
}

} // namespace obdeect::sources
