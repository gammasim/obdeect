# Refraction and timing

Air bends light and slows it down. This page describes how obdeect computes the bending of a ray and
the air delay from an atmospheric profile. The atmosphere stage uses it for light sources and for
photons emitted in the air. Photons of a run with `refraction_applied` = 1 are not refracted again.

The bending and air delay described on this page are plane-parallel (layers of constant height, no
Earth curvature). This is the expansion used in mode `none`. In mode `spherical` the same
quantities come from integrals over spherical shells ([](curvature.md#refraction-in-the-spherical-mode)).
In mode `corsika` CORSIKA has already applied the curved refraction.

## Input table

Lines starting with `#` and blank lines are skipped. The first four columns are used:

```text
# Alt [km]   rho [g/cm^3]   thick [g/cm^2]   n-1
    0.000    1.21171e-03     1.04473e+03    2.80503e-04
    1.000    1.10668e-03     9.28757e+02    2.55945e-04
```

Rules (violations raise `ValueError`):

- 5 to 120 levels, altitudes strictly ascending and within ±1000 km.
- density and n-1 positive; density and thickness strictly decreasing with height.
- top thickness may be zero; the top density must be positive.
- the thickness column must agree with the integral of the density over the table to 5 %.

## Interpolation

Heights in m a.s.l. (cm internally).

1. Density and n-1 are interpolated by a natural cubic spline of their logarithm through the levels.
2. The spline is tabulated on a fine grid: every level is a grid point, cells at most 3 m long.
3. Inside a cell density and n-1 are exponential (linear in the logarithm), so integrals of n-1 over
   height are exact for the tabulated profile.
4. Thickness between two levels is the upper level's table value plus the interval's table
   thickness times the fraction of the interval's mass (from the interpolated density) above the
   height. Level values are reproduced exactly.
5. Below the table the bottom values are used; above it density, thickness and n-1 are 0.

```python
profile = obdeect.propagation.AtmosphereProfile("atmprof_ecmwf_north_winter_fixed.dat")
profile.density([0.0, 2180.0])  # [1.21171e-03, 9.75149e-04] g/cm^3
profile.group_index(2180.0)  # 1.000225191893913
```

## Ray bending

Plane-parallel layers, n(h) = 1 + nu(h). A photon emitted at height `h_e` with zenith angle
`theta_e` (direction cosines `(cx, cy)`, downward) follows Snell's law

```text
n(h) sin(theta(h)) = n(h_e) sin(theta_e)
```

to the observation level, where `sin(theta_obs) = sin(theta_e) n(h_e) / n(obs)`. Over the path

```text
horizontal offset   dx = integral of tan(theta(h)) dh  -  (h_e - obs) tan(theta_e)
travel time         T  = (1/c) integral of n(h) / cos(theta(h)) dh
```

`dx` is the displacement against the straight line along the emission direction (negative: the ray
ends nearer the vertical of the emission point). `T` includes the vacuum path along the emission
direction, `(h_e - obs) / (c cos(theta_e))`; `T` minus that vacuum path is the air excess. The
refractive index n is used for the travel time (no dispersion term).

The integrals are evaluated as a series. With `delta = n_e/n - 1` (about 1e-4) the integrands are
power series in `s_e delta`, with coefficients from the Taylor series of `s/sqrt(1-s^2)` and
`1/sqrt(1-s^2)`. `delta^k` is a polynomial in `nu_e` and `nu`, so each term is a product of powers
of `nu_e` and the moments `M_b(h) = integral of nu^b dh`, tabulated once per profile and observation
level. Six terms are summed; the truncation error is below 1e-6 m for zenith angles up to 82 degrees
(tested). No loop over layers, no allocation, and two table lookups at most.

The new direction is `(cx, cy) sin(theta_obs)/sin(theta_e)`; the position moves by
`dx (cx, cy)/sin(theta_e)`.

`bend` returns `False` and leaves the photon unchanged if an input (height, direction, position,
time) is not finite, the direction is not inside the unit circle (upward or horizontal), the emission
is below the observation level, or the ray is too close to the horizon for the series
(`sin(theta_e) (nu_obs - nu_e) > 0.1 cos^2(theta_e)`, about 87 degrees; `bend_series_limit_exceeded`
in `refraction.hpp` is the one definition of this condition). Public height lookups, the timing
helpers and the observation level throw for non-finite values.

Example (ECMWF table, observation level 2180 m, photon direction purely in x):

| emission height | zenith | displacement | travel time minus vacuum straight path |
|---|---|---|---|
| 10 km | 20 deg | -0.18 m | +4.07 ns |
| 10 km | 60 deg | -3.04 m | -0.73 ns |
| 20 km | 60 deg | -9.32 m | -15.40 ns |
| 40 km | 60 deg | -12.59 m | -24.02 ns |

Vertical photons have no displacement. Their travel-time excess over vacuum equals
`straight_path_air_delay_ns` (4.02 ns from 10 km, 5.77 ns from 20 km).

At large zenith angles the excess can be negative: the ray bends towards the vertical, so its path
is shorter than the straight line, which outweighs the slower speed in air. The series agrees with a numerical integration of the exact integrands to 1e-6 m and 1e-6 ns on
the ECMWF table.

## Spherical mode

Mode `spherical` of the [curvature](curvature.md) option does not use the expansion above. It bends
the ray in spherical shells, with Bouguer's law $n r \sin\theta = K$, and adds the air delay along
the bent ray. The equations, the lookup tables and their accuracy are on the
[curvature page](curvature.md#refraction-in-the-spherical-mode). The plane-parallel bending above
is the limit of that mode for an infinite Earth radius.

## Dispersion

The profile's n-1 is given at a reference wavelength `lambda_ref`. Air refracts shorter wavelengths
more, so the phase index and the group index at wavelength `lambda` differ from the reference
profile by a factor. The ratios come from the Ciddor (1996) standard dry-air formula, with vacuum
wavelength in um and `s2 = 1 / lambda^2`:

```text
(n - 1) 1e8   = k1 / (k0 - s2) + k3 / (k2 - s2)            k0 = 238.0185, k1 = 5792105
                                                            k2 = 57.362,  k3 = 167917
(n_g - 1) 1e8 = (n - 1) 1e8 + 2 s2 [k1 / (k0 - s2)^2 + k3 / (k2 - s2)^2]
```

The group index follows from `n_g = n - lambda dn/dlambda` with `d/dlambda = -2 s2 / lambda d/ds2`.
Both are valid for 200 nm to 2000 nm; outside that range they throw `ValueError` (`std::invalid_argument`
in C++).

```text
k       = (n-1)(lambda) / (n-1)(lambda_ref)      phase factor:  n(lambda) - 1 = k (n - 1)
k_group = (n_g-1)(lambda) / (n-1)(lambda_ref)    group factor:  n_g(lambda) - 1 = k_group (n - 1)
```

Example, reference 500 nm: `k` = 1.045147, 1.013578, 0.988639 and `k_group` = 1.197287, 1.090693,
1.011885 at 300, 400 and 700 nm. At 300 nm the group correction is large (about 15 % of n-1).

Use: `obdeect propagation run --profile atmprof.dat --dispersion-reference-nm 500` (Python
`dispersion_reference_nm=500.0`). Each sky-frame photon that the stage bends (`refraction_applied = 0`)
is scaled to its own wavelength, and the wavelength is assigned before the geometry for these runs.
A run with `refraction_applied = 1` is rejected: the photons are already refracted and the straight
air delay is not scaled, so that `obdeect propagation to-sky` stays the exact inverse. The profile's
values are the index at `lambda_ref`; the ratios come from Ciddor (1996).

Plane-parallel bending (`bend` in `refraction.hpp`, overload with `k`, `k_group`): the series of
the [ray bending](#ray-bending) is evaluated with `nu -> k nu`. Each moment `M_b` is multiplied by
`k^b`, and `nu_e`, `nu_obs` by `k`. The geometry (position and direction) is exact for the scaled
profile up to the series truncation. The travel time is the series time with `k` plus the group
term

```text
T = T(k) + (k_group - k) straight_path_air_delay_ns(profile, obs, h_e, zenith_e)
```

The group term is first order: it treats `1 / cos(theta(h))` along the bent ray as `1 / cos(theta_e)`.
The error is of second order in the bending. The 7-argument `bend` calls the scaled overload with
`k = k_group = 1`, which gives the same result bit for bit.

Spherical mode uses the following first-order scaling below 75 degrees. With the
unscaled integrals `shift(1)`, `bent(1)` and `straight(1)` from the table:

```text
shift          = k shift(1)
bent           = k bent(1) + (k_group - k) straight(1)
straight       = k_group straight(1)
```

Between 75 and 80 degrees the result blends smoothly into a wavelength-dependent table. Above
80 degrees the table evaluates the ray geometry with the phase index and the travel time with
the group index. It interpolates on nine refractivity-scale nodes covering all combinations of
wavelength and reference wavelength from 200 to 2000 nm. The table is built on first use of a
high-zenith photon and cached per telescope altitude level for the lifetime of the propagator.
The direct quadrature fallback uses the same transition. Each cached altitude level requires
additional memory and construction time; telescopes at different levels require separate tables.
No dispersion tables are built when `--dispersion-reference-nm` is absent.

Measured accuracy (`test_dispersion.cpp`, fixture ECMWF north winter, observation level 2180 m,
reference numerical integrals with 20 000 Simpson intervals):

| case | error |
|---|---|
| plane, `k = k_group`, 10-30 km, 0-75 deg, 126 photons | displacement 1.7e-10 m, time 1.8e-9 ns |
| plane, `k_group != k` (Ciddor 300/400/700 nm, synthetic 0.97/0.98 and 1.03/1.02) | displacement 2.0e-10 m, time 6.3e-3 ns (worst, 300 nm) |
| plane, group term alone against `(k_group - k) straight` | 2.5e-11 ns |
| spherical, table interpolation alone (`k = 1`) | shift 8.4e-7 rel, bent 9.4e-6 rel |
| spherical, `k = 0.97, 1.03`, 0-85 deg, 10 and 30 km, against direct quadrature of the scaled profile | shift 7.0e-4 rel (1.3 m absolute at 85 deg, 30 km), bent 7.0e-4 rel, straight 1.6e-7 rel |
| spherical, `(k, k_group) = (0.97, 0.98)`, `(1.03, 1.02)` | bent 7.0e-4 rel, straight 1.6e-7 rel |

For the first-order calculation, the spherical error is second order in the scale: near the horizon the bent delay is not linear in
k, and the linear extrapolation is off by about 7e-4 of the total. Vertical photons are exact (bent =
straight = k straight(1)). The plane-parallel error is set by the group approximation alone and grows with `|k_group - k|`.

## Timing

- `ambient_group_index(profile, obs)`: refractive index `n = 1 + (n-1)` at the observation level.
  ln(n-1) is interpolated with a natural cubic spline in height (m), tabulated on 40 000 equidistant
  heights between the lowest and highest level, linear between them. Below the table the first value
  is used; at or above the top level `n = 1`. The height must be below the second-highest level
  (else `ValueError`).
- `propagation_time_ns(path, n) = path * n / c`, `c = 0.299792458 m/ns`.
- `straight_path_air_delay_ns(profile, h1, h2, zenith)`: `(1/c) integral of (n-1) dh / cos(zenith)`
  between two heights, exact for the interpolated profile.
- `straight_path_air_delay_from_cos_ns(profile, h1, h2, cos_zenith)`: the same from `cos(zenith)`
  (> 0). The scattering kernel uses it, since `acos` rounds to pi/2 for `|cz|` below about 1e-16.

## Units

| quantity | unit |
|---|---|
| height, position | m (public API); cm internally (`*_cm` methods) |
| density | g/cm³ |
| thickness | g/cm² |
| n-1 | dimensionless |
| direction | cosines (cx, cy); cz = sqrt(1 - cx² - cy²) > 0 |
| time | ns, added to the emission time |

Refraction tables depend on the observation level: build one `RefractionTables` per level.
`RefractionTables.etadsn` is the Gladstone-Dale ratio `(n-1)/density` (about 0.23 cm³/g) at the
first level more than 1.5 km above the observation level.

## Near the horizon

`bend_quadrature` (`obdeect/propagation/refraction_horizon.hpp`) evaluates the two integrals above
directly, from Snell's law, with Gauss-Legendre quadrature (8 nodes). Its outputs, units and
rejection rules are those of `bend`. The integrand is taken as a difference from the emission
value (`tan theta(h) - tan theta_e`, `n/cos theta(h) - 1/cos theta_e`), so no cancellation is
lost. Sub-intervals are at most 250 m, starting at 1 m near the emission. For `cos theta_e < 0.2`
the first sub-intervals use `x = u^2`, which removes the square-root behaviour of `1/cos theta`
at the emission level. It returns `False` only if `sin theta` reaches one at some height (total
internal reflection).

The stage uses it for the photons that `bend` rejects only at the series limit, so no photon is
dropped near the horizon:

```cpp
if (!bend(tables, h, cx, cy, x, y, t, k, k_group) &&
    bend_series_limit_exceeded(tables, h, cx, cy, k)) {
  bend_quadrature(tables, h, cx, cy, x, y, t, k, k_group); // exact, photon unchanged on failure
}
```

`bend_series_limit_exceeded` returns true only for valid inputs that violate the condition in
`bend`, so other failures (invalid input, emission below the observation level) are not sent to
the quadrature. The dispersion factors `k` and `k_group` (see [Dispersion](#dispersion)) enter both
functions with the same phase scaling. Quadrature integrates the group index along the bent ray;
the series uses the straight-path approximation above. Example (ECMWF table, 88 degrees, 20 km): displacement -27 705 m, air excess
-92 200 ns; the photon is bent towards the vertical by about 5 %.

Accuracy (measured with `test_horizon`, observation level 2180 m): for zenith 0 to 82 degrees at
10, 20 and 40 km the quadrature agrees with `bend` to 3e-8 m and 9e-8 ns. For 87 to 89 degrees,
halving the largest sub-interval (250 to 125 m) changes the displacement by at most 4e-6 m and the
time by at most 1.4e-5 ns. The cost is 8 profile lookups per sub-interval, about 1300 lookups for
a 40 km emission.
