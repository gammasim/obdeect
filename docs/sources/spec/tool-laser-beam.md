# obdeect sources laser-beam (implementation)

User options and examples: [laser-beam](../programs/laser-beam.md). Writer B, producer
`obdeect sources laser-beam`. Only primary photons are generated; scattering is in
[obdeect propagation](../../propagation/index.md).

## Validation

- `--step` metres, `> 0` and finite. `--primary-per-m` `>= 0`, finite (a sampling choice).
- `--brightness` is rejected: `--brightness is replaced by --primary-per-m`.
- `--laser-position x,y,z`; `x,y` allowed (z = 0), else `--laser-position needs x,y or x,y,z`.
- `--fov`: `0 < fov < 90`, else `laser-beam: field of view must lie in (0, 90) deg`.
- `--angular`, `--angular-distribution`: stored verbatim, not interpreted.
- Astronomical azimuth `a` converts to the mathematical azimuth by `-a` if `a <= 0`, else
  `360 - a` (both `--laser-phi` and `--telescope-phi`). The C++ option struct holds the
  mathematical values (defaults 180 and 0).

## Geometry

- Beam direction `dir = direction_from_altitude((90 - laser_theta) deg, laser_phi deg)`,
  so `dir.z = cos(theta)`.
- Telescope-axis frame: `p = rotate_y(rotate_z(laser_position, -phi_t), -theta_t)`,
  `q = rotate_y(rotate_z(dir, -phi_t), -theta_t)`.
- Cone of half angle `h = fov/2`, `T = tan(h)`. Points `p + t q` on the cone satisfy
  `a t^2 + b t + c = 0` with `a = qx^2 + qy^2 - qz^2 T^2`, `b = 2 (qx px + qy py - qz pz T^2)`,
  `c = px^2 + py^2 - pz^2 T^2`.
- Roots `t1 <= t2` (numerically stable, relative accuracy 1e-12). If `a = 0` and `b != 0`, the
  root is `-c/b`.
- Interval: `a > 0`: `[max(t1, 0), t2]`; `a = 0`, `b > 0`, `c < 0`: `[0, -c/b]`.
- Rejected: negative discriminant or no forward interval (`laser-beam: laser beam does not
  cross the field of view`); `a < 0` or an unbounded linear case.
- Beam start (sky frame, relative to the telescope): `laser_position + dir x interval_start`.
  The nappe is not checked.

## Emissions

- Steps `s_k` by repeated addition in double: `s = 0; while s < L: emit k; s += step`
  (about `ceil(L/step)` steps).
- Step `k`: emitter at `start + dir s_k`, pointing `dir`, FixedAngle(0), frame sky,
  `Source.photons = primary_per_m x step / b`, `b = max(1, bunchsize)`, seed `seed + k`,
  bunch weight `b`, time offset `s_k / c` (float32 after the pulse).
- One case, `events` events. Expected total weight `primary_per_m x L` within `4 b sqrt(steps) / 2`.

## Run record

- One telescope id 1 at (0, 0, 0), radius `--telescope-radius`; wavelength range from the spectrum.
- Extra keys: `scattering_angular` (string = angular text), `primary_photons_per_m` (f64),
  `beam_step_m`, `beam_path_m` (= L), `fov_deg`, `telescope_zenith_deg`,
  `telescope_azimuth_deg` (f64; the last two are the stored mathematical-sense values).

## Acceptance (`test_laser-beam`)

- Geometry at default; step count; time offsets `k step / c` to 1e-9.
- `cz` of all photons `= dir.z`, wavelength 355, `w = 3`.
- One event; `scattering_angular == "Rayleigh"`; `beam_path_m == L`.
- Gaussian pulse residual moments (sigma 0.5 ns, fov 20, step 0.5, bunch 1, `1e3` per m).
- fov in `{-15, 0, 90, 180, nan, inf}` rejected; `{15, 89.9}` accepted.
- Geometry: both ends have angle `h` to the axis (cosine to 1e-9); misses, behind-laser
  intersections and unbounded intervals are rejected.
