// The shared 1/lambda^2 wavelength draw: range, moments, continuity, quantiles.

#include "check.hpp"

#include "obdeect/stream/random.hpp"
#include "obdeect/stream/wavelength.hpp"

#include <cmath>
#include <cstdint>

int main() {
  return run_checks([] {
    using obdeect::stream::inverse_square_wavelength_nm;
    const double lo = 300.0, hi = 600.0;
    CHECK_NEAR(inverse_square_wavelength_nm(0.0, lo, hi), lo, 1e-9);
    // u = 1/2 gives the harmonic midpoint of the range.
    CHECK_NEAR(inverse_square_wavelength_nm(0.5, lo, hi), 2.0 * lo * hi / (lo + hi), 1e-9);

    // Counter-based draws: mean and variance of density (1/l^2)/(1/lo - 1/hi) within 5 sigma.
    const double norm = 1.0 / lo - 1.0 / hi;
    const double mean = std::log(hi / lo) / norm;
    const double variance = (hi - lo) / norm - mean * mean;
    const std::uint64_t n = 400000;
    double sum = 0.0, sum2 = 0.0;
    bool inside = true, continuous = false;
    for (std::uint64_t i = 0; i < n; ++i) {
      const obdeect::stream::DrawKey key{1, 2, 3, 4, i, 5};
      const double w = inverse_square_wavelength_nm(obdeect::stream::uniform4(key)[0], lo, hi);
      inside = inside && w >= lo && w <= hi;
      continuous = continuous || w != std::floor(w);
      sum += w;
      sum2 += w * w;
    }
    const double sample_mean = sum / double(n);
    const double sample_var = sum2 / double(n) - sample_mean * sample_mean;
    CHECK(inside);
    CHECK(continuous); // not rounded to integer nm
    CHECK(std::fabs(sample_mean - mean) < 5.0 * std::sqrt(variance / double(n)));
    CHECK(std::fabs(sample_var - variance) < 0.02 * variance);
    // Share below 400 nm: F(400) = (1/lo - 1/400) / norm.
    std::uint64_t below = 0;
    for (std::uint64_t i = 0; i < n; ++i) {
      const obdeect::stream::DrawKey key{1, 2, 3, 4, i, 5};
      below += inverse_square_wavelength_nm(obdeect::stream::uniform4(key)[0], lo, hi) < 400.0;
    }
    const double p = (1.0 / lo - 1.0 / 400.0) / norm;
    CHECK(std::fabs(double(below) / double(n) - p) < 5.0 * std::sqrt(p * (1.0 - p) / double(n)));
  });
}
