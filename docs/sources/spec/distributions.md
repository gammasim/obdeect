# Distributions

Header `obdeect/sources/density.hpp`, namespace `obdeect sources`. `Uniforms = std::array<double,4>`; `sample(u)`
of 1-D shapes uses `u[0]` (complex pulse `u[0..2]`).

## LinearDensity

Unnormalised density, piecewise linear between nodes `(x_i, y_i)`, zero outside
`[x_0, x_n]`. Constructor error (`invalid_argument`): fewer than 2 nodes, size mismatch, a
non-finite value, `y < 0`, non-increasing `x`, or integral `<= 0`. Default-constructed =
empty = point mass at 0 (`sample` returns 0, `cdf` 0, `integral` 0).

- `cdf(x)`: integral from `x_0` to `x` (trapezoid per segment; within a segment
  `C_i + y_i s + k s^2/2`, `s = x - x_i`, `k` the slope); 0 below, total above.
- `sample(p)` for `p in (0,1)`: the `x` with `cdf(x) = p * integral` (solve the segment
  quadratic stably for any slope, clamp to the segment).
- `at(x)`, `maximum()`, `integral()`, `lower()`, `upper()`.
- Example: nodes `{0,1,3}`, values `{0,2,2}`: integral 5, `sample(0.1) = sqrt(0.5)`,
  `sample(0.6) = 2`, `sample(1.0) = 3`, `cdf(2) = 3`.

## Normal quantile

`normal_quantile(p)`: inverse standard-normal CDF to relative error ~1e-15 (reference:
`0.975 -> 1.959963984540054`, `1e-6 -> -4.753424308822899`, `0.5 -> 0`). Used by Gaussian
shapes: `x = mean + sigma * normal_quantile(u)`.

## 1-D shapes (pulses in ns, spectra in nm)

| Shape | Parameters | Sample from `u` | Mean / variance |
| --- | --- | --- | --- |
| Fixed | value | `value` | value / 0 |
| UniformRange | low, high | `low + (high-low) u0` | (low+high)/2 / (high-low)^2/12 |
| GaussianShape | mean, sigma | `mean + sigma * normal_quantile(u0)` | mean / sigma^2 |
| Exponential | mean, tau | `mean + tau * (-ln(1-u0) - 1)` (leading edge at mean - tau) | mean / tau^2 |
| ComplexPulse | flat, rms, decay, delay | `delay + flat (u0-0.5) + rms * normal_quantile(u1) + decay (-ln(1-u2) - 1)` | delay / flat^2/12 + rms^2 + decay^2 |
| InverseSquare | low, high | `1 / (1/low - u0 (1/low - 1/high))` (density ~ 1/lambda^2) | `ln(high/low)/(1/low-1/high)` / `(high-low)/(1/low-1/high) - mean^2` |
| Tabulated | LinearDensity | `density.sample(u0)` | from table |

`Spectrum = Fixed | UniformRange | InverseSquare | Tabulated`;
`Pulse = Fixed | UniformRange | GaussianShape | Exponential | ComplexPulse | Tabulated`. All
members `sample`, `mean`, `variance` (not Tabulated) are public. Moments are tested to four
sigma (`test_density`).

## Angular profiles

Argument `u = 1 - cos(theta)`. Each profile has `intensity(u)` (relative per solid angle),
`maximum()`, `integral()` (over `u`), `cdf(u)` (emitted fraction below `u`), `quantile(p)`.

| Profile | intensity | integral | cdf | quantile |
| --- | --- | --- | --- | --- |
| Isotropic | 1 on [0,2] | 2 | clamp(u/2,0,1) | 2p |
| FixedAngle(u0) | none (0) | 0 | 1 if u >= u0 else 0 | u0 |
| FilledCone(u0) | 1 on [0,u0] | u0 | clamp(u/u0,0,1) | p u0 |
| Lambertian | 1-u on [0,1] | 1/2 | 2v - v^2, v = clamp(u,0,1) | 1 - sqrt(1-p) |
| TabulatedAngular(LinearDensity over u) | table | table | cdf/integral | density.sample(p) |

`Angular` is the variant of these five. Helpers `intensity_maximum`, `angle_quantile`,
`integral(Angular)`.

- `gaussian_angular_table(sigma_rad)`: intensity `exp(-theta^2/(2 sigma^2))`, nodes at
  `theta = 0, 0.05 sigma, 0.10 sigma, ...` below `theta_max = min(5 sigma, pi)`
  (nodes computed by repeated addition of `0.05 sigma` starting at 0, while `theta < theta_max`)
  plus a last node at `theta_max`; `u = 1 - cos(theta)`. Error if `sigma < 1e-12` rad or
  non-finite (message mentions 1e-12 rad). `sigma` is the 2-D r.m.s. in angle: the fitted
  `sqrt(mean(theta^2)/2)` equals sigma (test: sigma 5 deg within 2%). Intensity at `u(sigma)` is
  `exp(-1/2)` to 1e-3 and at `u = 0` is 1.
- `rayleigh_angular_table()`: intensity `0.5 (1 + cos^2 theta)`, nodes every 1 deg from 0 to
  180 deg (181 nodes), `u = 1 - cos(theta)`.
