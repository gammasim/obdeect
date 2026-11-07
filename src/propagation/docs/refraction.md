# Refraction and timing

dlufd reads an atmospheric table in the CORSIKA / sim_telarray `atmprof` format and
reproduces the interpolation and ray-bending approximations of CORSIKA's `atmo.c`
(IACT/atmo, K. Bernlöhr, LGPL 2.1+). Light sources and photons emitted in the air use this
path; CORSIKA air-shower photons are already refracted by CORSIKA.

## Input table

Lines starting with `#` and blank lines are skipped. The first four columns are used:

```text
# Alt [km]   rho [g/cm^3]   thick [g/cm^2]   n-1
    0.000    1.21171e-03     1.04473e+03    2.80503e-04
    1.000    1.10668e-03     9.28757e+02    2.55945e-04
```

Rules (violations raise `ValueError`):

- 5 to 120 levels, altitudes strictly ascending.
- density and n-1 positive; density and thickness strictly decreasing with height.
- top thickness may be zero; the top density must be positive.
- the density and thickness columns must agree to 5 % (integrated from the levels).

## Interpolation (atmo.c)

Heights are m a.s.l. Internally atmo.c uses cm.

1. Each column is interpolated by a natural cubic spline in ln(value) against altitude
   (`rpolator.c`, `set_1d_cubic_params`). Thickness uses the levels without the top one.
2. The spline is sampled on 40 000 equidistant heights from the bottom to the top of the
   table (`init_fast_interpolation`, `atmo.c:553`).
3. Lookups interpolate linearly on that grid (`rhofx_` `atmo.c:1393`, `refim1x_` `:1513`,
   `thickx_` `:1446`).
4. Above the last but one level, thickness is the exponential top layer
   `t(h) = rho0 cfac hscale (exp(-h/hscale) - exp(-top/hscale))`, and the grid n-1 there is
   `etadsn * rho`.

With the ECMWF north winter table (`cpp/tests/data/atmprof_ecmwf_north_winter_fixed.dat`):

```python
profile = dlufd.AtmosphereProfile("atmprof_ecmwf_north_winter_fixed.dat")
profile.density([0.0, 2180.0])        # [1.21171e-03, 9.75149e-04] g/cm^3
profile.group_index(2180.0)           # 1.000225191893913
```

## Ray bending (atmo.c `raybnd_`)

For a photon emitted at height `h` with direction cosines `(cx, cy)` (downward), CORSIKA
computes the zenith angles at emission and at the observation level

```text
sin(theta_obs) = sin(theta_em) * n(h) / n(obs)
```

and, with `c = cos(theta_em + 0.28 (theta_obs - theta_em))`, `s` the matching sine:

```text
displacement  d = -(H^2) s / c^3                        H from the density table (cm)
travel time   T = T0(h) - (T_a^2) s^2 / c^3 + (dz + de) / (c_vac cos(theta_em))
```

- `T0(h)` (ns) is the altitude table of the vertical air-time excess: air time minus its straight-line estimate (`init_refraction_tables`).
- `H` (cm) and `T_a` (ns) are density tables from the same 45 degree integration (`trace_ray_planar`, 2.5 m step).
- `dz = |h - obs|`, `de = etadsn (thick(obs) - thick(h))`, `etadsn` is `(n-1)/rho` at the first level more than 1.5 km above the observation level.
- The new direction is `(cx, cy) sin(theta_obs)/sin(theta_em)`; the position moves by `d (cx, cy)/sin(theta_em)`.

Photons with a direction outside the unit circle (upward or horizontal) are returned
unchanged with `bent = False`.

Example (ECMWF table, observation level 2180 m, photon direction purely in x):

| emission height | zenith | displacement | travel time minus vacuum straight path |
|---|---|---|---|
| 10 km | 20 deg | -0.18 m | +4.07 ns |
| 10 km | 60 deg | -3.04 m | -0.73 ns |
| 20 km | 60 deg | -9.32 m | -15.40 ns |
| 40 km | 60 deg | -12.59 m | -24.02 ns |

Vertical photons have no displacement. Their travel-time excess over vacuum equals
`straight_path_air_delay_ns` exactly (4.02 ns from 10 km, 5.77 ns from 20 km).

Caution: the 60 degree time term is CORSIKA's second-order fit to the 45 degree integration,
extrapolated. At 20 km and above it gives negative excess times, which the air delay
(a few ns) does not explain. dlufd reproduces CORSIKA here; the term has not been checked
against full ray tracing, so timing at large zenith angles needs a separate validation
before it is used for physics.

## Timing

Following obdeect (`propagation.py`, `trace.hpp`):

- `ambient_group_index(profile, obs)`: natural cubic spline of ln(n-1) in height (m), on
  40 000 equidistant heights, linear in ln(n-1); returns `1 + (n-1)` at the observation
  level. obdeect names this the group index; it is the refractive index at that level.
- `propagation_time_ns(path, n) = path * n / c`, `c = 0.299792458 m/ns`.
- `straight_path_air_delay_ns(profile, h1, h2, zenith)`: `(1/c) ∫ (n-1) dh / cos(zenith)`
  between two heights, integrated exactly over the piecewise-linear grid.

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

## Accuracy

Checked against CORSIKA's `atmo.c` (built with `-DDLUFD_CORSIKA_REFERENCE=<corsika7/bernlohr>`,
`ctest`, test `corsika_reference`): density, thickness and n-1 agree to 1e-12 relative at
all tested heights, and bent positions, directions and times agree to 1e-9 m, 1e-12 and
1e-6 ns. Lookups reproduce the table nodes to about 2e-5 relative (density), not exactly,
because the grid is 3 m. The coarse 1 km synthetic table in the tests reproduces the
analytic thickness to 0.4 % in its top 2 km.

## References

- K. Bernlöhr, `corsika7/bernlohr/atmo.c` (IACT/atmo package for CORSIKA), functions
  `init_refraction_tables`, `trace_ray_planar`, `raybnd_`, `rhofx_`, `thickx_`, `refim1x_`.
- `rpolator.c`, `set_1d_cubic_params` and `rpol_cspline` (cubic spline interpolation).
- sim_telarray `cfg/CTA/atmprof_*.dat` (ECMWF profiles; column description in the file header).
- obdeect `python/obdeect/propagation.py` (`ambient_group_index`) and `cpp/include/obdeect/trace.hpp`.
