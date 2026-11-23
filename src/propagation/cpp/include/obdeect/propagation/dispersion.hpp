#pragma once

// Wavelength dependence (dispersion) of air refractivity. See docs/propagation/refraction.md
// ("Dispersion").
//
// Ciddor (1996) standard dry air, vacuum wavelength lambda in um, s2 = 1 / lambda^2:
//   (n - 1) 1e8 = k1 / (k0 - s2) + k3 / (k2 - s2)
// with k0 = 238.0185, k1 = 5792105, k2 = 57.362, k3 = 167917. The group index n_g = n - lambda
// dn/dlambda follows from d/dlambda = -2 s2 / lambda d/ds2 (lambda in um): (n_g - 1) 1e8 = N + 2 s2
// dN/ds2 with dN/ds2 = k1 / (k0 - s2)^2 + k3 / (k2 - s2)^2.
//
// The profile of n-1 is given at a reference wavelength. The ratios below scale it to wavelength
// lambda: n(lambda) - 1 = k (n(ref) - 1) for the phase index, and n_g(lambda) - 1 = k_group
// (n(ref) - 1) for the group index. Only the ratio of the air composition is used, so the
// profile's absolute values carry the air conditions.

#include <cmath>
#include <stdexcept>

namespace obdeect::propagation {

namespace detail {

inline constexpr double ciddor_k0 = 238.0185;
inline constexpr double ciddor_k1 = 5792105.;
inline constexpr double ciddor_k2 = 57.362;
inline constexpr double ciddor_k3 = 167917.;
inline constexpr double dispersion_min_nm = 200.;
inline constexpr double dispersion_max_nm = 2000.;

} // namespace detail

// Validity range of the dispersion formulas [nm].
inline constexpr double refractivity_min_nm = detail::dispersion_min_nm;
inline constexpr double refractivity_max_nm = detail::dispersion_max_nm;

// (n - 1) 1e8 of standard dry air at vacuum wavelength lambda_nm [nm]. Throws std::invalid_argument
// outside 200 to 2000 nm or for non-finite input.
[[nodiscard]] inline double refractivity_e8(double lambda_nm) {
  if (!std::isfinite(lambda_nm) || lambda_nm < detail::dispersion_min_nm ||
      lambda_nm > detail::dispersion_max_nm) {
    throw std::invalid_argument("wavelength must be finite and within 200 to 2000 nm");
  }
  const double lambda_um = 1e-3 * lambda_nm;
  const double s2 = 1. / (lambda_um * lambda_um);
  return detail::ciddor_k1 / (detail::ciddor_k0 - s2) +
         detail::ciddor_k3 / (detail::ciddor_k2 - s2);
}

// (n_g - 1) 1e8 of standard dry air at vacuum wavelength lambda_nm [nm]; same validation.
[[nodiscard]] inline double group_refractivity_e8(double lambda_nm) {
  if (!std::isfinite(lambda_nm) || lambda_nm < detail::dispersion_min_nm ||
      lambda_nm > detail::dispersion_max_nm) {
    throw std::invalid_argument("wavelength must be finite and within 200 to 2000 nm");
  }
  const double lambda_um = 1e-3 * lambda_nm;
  const double s2 = 1. / (lambda_um * lambda_um);
  const double a = 1. / (detail::ciddor_k0 - s2);
  const double b = 1. / (detail::ciddor_k2 - s2);
  const double n_minus_one = detail::ciddor_k1 * a + detail::ciddor_k3 * b;
  const double d_ds2 = detail::ciddor_k1 * a * a + detail::ciddor_k3 * b * b;
  return n_minus_one + 2. * s2 * d_ds2;
}

// Phase factor k = (n-1)(lambda) / (n-1)(reference): n(lambda) - 1 = k (n - 1) of the profile.
[[nodiscard]] inline double refractivity_ratio(double lambda_nm, double reference_nm) {
  return refractivity_e8(lambda_nm) / refractivity_e8(reference_nm);
}

// Group factor k_group = (n_g - 1)(lambda) / (n - 1)(reference) for the travel time.
[[nodiscard]] inline double group_refractivity_ratio(double lambda_nm, double reference_nm) {
  return group_refractivity_e8(lambda_nm) / refractivity_e8(reference_nm);
}

} // namespace obdeect::propagation
