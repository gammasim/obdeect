# Background source (night sky)

Header `obdeect/sources/background.hpp`. Photons start on the sphere of one telescope on the hemisphere
facing their own direction, travel inwards, with directions uniform in solid angle inside a
cone about a viewing axis.

`BackgroundSource{centre, radius_m, pointing (default (0,0,-1)), viewcone_rad (default pi/2),
spectrum (Fixed 400), pulse_width_ns=0, pulse_delay_ns=0}`.

- `background_count(mean_bunches, seed, event, telescope)`: rounding rule of conventions.
- `sample_background(s, seed, event, telescope, index)`:
  - direction stream: `cos_theta = 1 - u0 (1 - cos(viewcone))`, `phi = 2 pi u1`; direction `d =
    pointing_unit cos_theta + (e1 cos phi + e2 sin phi) sin_theta` (basis of `pointing`).
  - target stream: `ri = 1 - sqrt(u0)`, `phi' = 2 pi u1`, `zi = sqrt(max(0, 1 - ri^2))`; the
    local point `(ri cos phi', ri sin phi', zi)` is rotated so that its z axis is `-d`:
    `za = acos(clamp(-d.z,-1,1))`, `az = atan2(-d.y, -d.x)`, `unit = rotate_z(rotate_y(local, za),
    az)`. Position `centre + unit * radius`. Then `dot(unit, d) <= 0` for every photon
    (tested), and `E[zi] = mean over u of sqrt(1 - (1-sqrt(u))^2)`.
  - time `= pulse_delay_ns - zi radius/c + w`, with `w = pulse_width_ns/2.35482004503095 *
    (pulse_u0 - 0.5)` if `pulse_width_ns > 0`, else 0 (pulse stream).
  - wavelength from the spectrum (wavelength stream).
- `emit_background(s, seed, event, telescope, first, rows, weight, chunk)`: requires `radius > 0`
  and `0 < viewcone <= pi` (else `invalid_argument`) and `rows <= 65536`; frame sky; every row
  weight `float(weight)` (the bunch size); columns `source_columns`.
