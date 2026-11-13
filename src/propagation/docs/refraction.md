# Refraction and timing

obdeect propagation reads an atmospheric table in the `atmprof` format and computes interpolation, ray
bending and air delay from it. Light sources and photons emitted in the air use this path;
photons of a run with `refraction_applied` = 1 are not refracted again.

## Input table

Lines starting with `#` and blank lines are skipped. The first four columns are used:

```text
# Alt [km]   rho [g/cm^3]   thick [g/cm^2]   n-1
    0.000    1.21171e-03     1.04473e+03    2.80503e-04
    1.000    1.10668e-03     9.28757e+02    2.55945e-04
```

Rules (violations raise `ValueError`):

- 5 to 120 levels, altitudes strictly ascending and within +-1000 km.
- density and n-1 positive; density and thickness strictly decreasing with height.
- top thickness may be zero; the top density must be positive.
- the thickness column must agree with the integral of the density over the table to 5 %.

## Interpolation

Heights are m a.s.l. (cm internally).

1. Density and n-1 are interpolated by a natural cubic spline of their logarithm through the
   levels.
2. The spline is tabulated on a fine grid: every level is a grid point and the cells are at
   most 3 m long.
3. Inside a cell density and n-1 are exponential (linear in the logarithm), so integrals of
   n-1 over height are exact for the tabulated profile.
4. Thickness between two levels is the table value of the upper level plus the table
   thickness of the interval times the fraction of the interval's mass (from the interpolated
   density) above the height. Level values are reproduced exactly.
5. Below the table the bottom values are used; above it density, thickness and n-1 are 0.

```python
profile = obdeect.propagation.AtmosphereProfile("atmprof_ecmwf_north_winter_fixed.dat")
profile.density([0.0, 2180.0])        # [1.21171e-03, 9.75149e-04] g/cm^3
profile.group_index(2180.0)           # 1.000225191893913
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

`dx` is the displacement against the straight line along the emission direction (it is
negative: the ray ends nearer the vertical of the emission point). `T` includes the vacuum
path along the emission direction, `(h_e - obs) / (c cos(theta_e))`, and the difference
`T - vacuum path` is the air excess. The refractive index n is used for the travel time (no
dispersion term).

The integrals are evaluated with a series. With `delta = n_e/n - 1` (about 1e-4) the
integrands are power series in `s_e delta`, whose coefficients follow from the Taylor series of
`s/sqrt(1-s^2)` and `1/sqrt(1-s^2)`. `delta^k` is a polynomial in `nu_e` and `nu`, so each term
is a product of powers of `nu_e` and the moments `M_b(h) = integral of nu^b dh`, which are
tabulated once per profile and observation level. Six terms are summed; the truncation error is
below 1e-6 m for zenith angles up to 82 degrees (tested). A photon needs no loop over layers, no
allocation and no random access beyond two table lookups.

The new direction is `(cx, cy) sin(theta_obs)/sin(theta_e)`; the position moves by
`dx (cx, cy)/sin(theta_e)`.

`bend` returns `False` and leaves the photon unchanged if an input (height, direction,
position, time) is not finite, the direction is not inside the unit circle (upward or
horizontal), the emission is below the observation level, or the ray is too close to the
horizon for the series (`sin(theta_e) (nu_obs - nu_e) > 0.1 cos^2(theta_e)`, about 87 degrees).
Public height lookups, the timing helpers and the observation level throw for non-finite values.

Example (ECMWF table, observation level 2180 m, photon direction purely in x):

| emission height | zenith | displacement | travel time minus vacuum straight path |
|---|---|---|---|
| 10 km | 20 deg | -0.18 m | +4.07 ns |
| 10 km | 60 deg | -3.04 m | -0.73 ns |
| 20 km | 60 deg | -9.32 m | -15.40 ns |
| 40 km | 60 deg | -12.59 m | -24.02 ns |

Vertical photons have no displacement. Their travel-time excess over vacuum equals
`straight_path_air_delay_ns` (4.02 ns from 10 km, 5.77 ns from 20 km).

The excess time can be negative at large zenith angles: the ray is bent towards the vertical,
so its path is shorter than the straight line along the emission direction, which outweighs
the slower speed in air. Cross-check: the tests integrate the exact integrands numerically and
agree with the series to 1e-6 m and 1e-6 ns on the ECMWF table.

## Timing

- `ambient_group_index(profile, obs)`: refractive index `n = 1 + (n-1)` at the observation
  level. ln(n-1) is interpolated with a natural cubic spline in height (m), tabulated on 40 000
  equidistant heights between the lowest and highest level, and linear between them. Below the
  table the first value is used, at or above the top level `n = 1`. The height must be below the
  second-highest level (else `ValueError`).
- `propagation_time_ns(path, n) = path * n / c`, `c = 0.299792458 m/ns`.
- `straight_path_air_delay_ns(profile, h1, h2, zenith)`: `(1/c) integral of (n-1) dh / cos(zenith)`
  between two heights, exact for the interpolated profile.

## Units and conventions

| quantity | unit |
|---|---|
| height, position | m (public API); cm internally (`*_cm` methods) |
| density | g/cm^3 |
| thickness | g/cm^2 |
| n-1 | dimensionless |
| direction | cosines (cx, cy); cz = sqrt(1 - cx^2 - cy^2) > 0 |
| time | ns, added to the emission time |

Refraction tables depend on the observation level; build one `RefractionTables` per level.
`RefractionTables.etadsn` is the Gladstone-Dale ratio `(n-1)/density` (about 0.23 cm^3/g) at
the first level more than 1.5 km above the observation level.

## Tests

Analytic checks in `cpp/tests/test_physics.cpp`: constant n gives a straight line and
`T = L n / c`; exponential atmosphere with analytic integrals of density, thickness and n-1;
the small-angle limit; rotation symmetry; monotonic growth with zenith angle and height;
agreement with direct numerical integration of Snell's law (synthetic and ECMWF profile).

Example profile: `cpp/tests/data/atmprof_ecmwf_north_winter_fixed.dat`.
