#pragma once

// Wavelength draw for photons without a wavelength: density proportional to 1/lambda^2 on
// [lo, hi] nm by inverse CDF, continuous (not rounded to integer nm). One implementation for
// propagation, ray tracing and sources.

#include <cmath>
#include <stdexcept>

namespace obdeect::stream {

// Inverse CDF of the 1/lambda^2 density on [lo, hi] for u in [0, 1). Throws std::invalid_argument
// unless 0 < lo_nm <= hi_nm and hi_nm is finite (a few comparisons per draw). lo_nm == hi_nm is a
// valid degenerate range (propagation's WavelengthRange allows it): every draw is lo_nm.
[[nodiscard]] inline double inverse_square_wavelength_nm(double u, double lo_nm, double hi_nm) {
  if (!(lo_nm > 0.0 && hi_nm >= lo_nm && std::isfinite(hi_nm))) {
    throw std::invalid_argument("wavelength range must satisfy 0 < lo <= hi (finite), in nm");
  }
  const double a = 1.0 / lo_nm;
  return 1.0 / (a - u * (a - 1.0 / hi_nm));
}

} // namespace obdeect::stream
