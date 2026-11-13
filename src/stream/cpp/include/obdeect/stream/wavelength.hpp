#pragma once

// Wavelength draw for photons without a wavelength: density proportional to 1/lambda^2 on
// [lo, hi] nm by inverse CDF, continuous (not rounded to integer nm). One implementation for
// propagation, ray tracing and sources.

namespace obdeect::stream {

// Inverse CDF of the 1/lambda^2 density on [lo, hi] for u in [0, 1).
[[nodiscard]] inline double inverse_square_wavelength_nm(double u, double lo_nm, double hi_nm) {
  const double a = 1.0 / lo_nm;
  return 1.0 / (a - u * (a - 1.0 / hi_nm));
}

} // namespace obdeect::stream
