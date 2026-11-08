# dlamb-laser-beam

Primary photons of a laser beam inside the field-of-view cone of a telescope, sky frame.
Scattering is not generated here. Writer B, producer `dlamb-laser-beam`.

## Options

| Option (aliases) | Default | Rules |
| --- | --- | --- |
| `--events n` | 1 | integer `0..2^32-1` |
| `--step m` (`--steplength`) | 0.1 | metres, `> 0` and finite |
| `--bunchsize f` | 3 | below 1 used as 1 |
| `--primary-per-m f` | `1e4` | primary photons per metre; `>= 0`, finite (a sampling choice) |
| `--brightness` | | rejected: `--brightness is replaced by --primary-per-m` |
| `--spectrum` | `355` | number 200..2000 nm, `Cherenkov:a:b` or table |
| `--lightpulse` (`--pulse`) | `Simple:0` | bare number = Simple width |
| `--angular`, `--angular-distribution` | `Rayleigh` | any text, stored verbatim, not interpreted |
| `--run` | 1 | written run number is `run` if `> 0`, else 1 |
| `--altitude` | 2150 | |
| `--atmosphere` | 1 | recorded as `atmosphere_id` |
| `--laser-position x,y,z` | `1000,0,0` | metres relative to the telescope; `x,y` allowed (z = 0); else `--laser-position needs x,y or x,y,z` |
| `--laser-theta` | 60 | zenith angle of the beam [deg] |
| `--laser-phi` | 180 | azimuth, astronomical sense |
| `--telescope-theta`, `--telescope-phi` | 30, 0 | telescope axis |
| `--telescope-radius` | 12.5 | m |
| `--fov` | 15 | full opening angle [deg], `0 < fov < 90`, else `laser-beam: field of view must lie in (0, 90) deg` |
| `--seed`, `-o`, `--help` | | |

Astronomical azimuth `a` converts to the mathematical azimuth by `-a` if `a <= 0` else `360 -
a` (both `--laser-phi` and `--telescope-phi`; the C++ option struct holds mathematical values,
defaults 180 and 0).

## Geometry

Beam direction `dir = direction_from_altitude((90 - laser_theta) deg, laser_phi deg)` (unit,
so `dir.z = cos(theta)`). Telescope-axis frame: `p = rotate_y(rotate_z(laser_position,
-phi_t), -theta_t)`, `q = rotate_y(rotate_z(dir, -phi_t), -theta_t)`. Cone of half angle
`h = fov/2`, `T = tan(h)`: points `p + t q` on the cone satisfy `a t^2 + b t + c = 0` with
`a = qx^2 + qy^2 - qz^2 T^2`, `b = 2 (qx px + qy py - qz pz T^2)`, `c = px^2 + py^2 - pz^2
T^2`. Roots `t1 <= t2` (numerically stable evaluation, relative accuracy 1e-12; if `a = 0`
and `b != 0` both roots are `-c/b`; `a = 0 = b`, or negative discriminant, is an error
`laser-beam: laser beam does not cross the field of view`). Beam start (sky frame, relative to
the telescope) `laser_position + dir t1`, inside-cone path length `L = t2 - t1`. The nappe is
not checked. Test: both ends have angle `h` to the axis (cosine to 1e-9); a beam that misses
the cone (laser at (0,1000,0), theta 90) is rejected.

## Emissions

Steps `s_k` by repeated addition in double: `s = 0; while s < L: emit k; s += step; k += 1`
(about `ceil(L/step)` steps, tested to within 1). Step `k`: Emitter at `start + dir s_k`,
pointing `dir`, FixedAngle(0), frame sky, spectrum, pulse, `Source.photons =
primary_per_m * step / b` with `b = max(1, bunchsize)`, seed `seed + k`, bunch weight `b`,
time offset `s_k / c` (added in float32 after the pulse). One case, `events` events. Expected
total weight `primary_per_m * L` within `4 b sqrt(steps)/2`.

## Run record

Writer B: one telescope id 1 at (0,0,0), radius `--telescope-radius`; wavelength range from the
spectrum; `atmosphere_id` as given; extra keys: `scattering_angular` (string = angular
text), `primary_photons_per_m` (f64), `beam_step_m`, `beam_path_m` (= L), `fov_deg`,
`telescope_zenith_deg`, `telescope_azimuth_deg` (f64; the last two are the stored
mathematical-sense values).

## Acceptance (`test_laser-beam`)

Geometry at default; step count; time offsets `k * step / c` to 1e-9 (steps use `k*step` within
1e-9 of the repeated sum); `cz` of all photons `= dir.z`, wavelength 355, `w = 3`; stream has
one event and the three keys `scattering_angular == "Rayleigh"`, `beam_path_m == L`; Gaussian
pulse residual moments (sigma 0.5 ns, fov 20, step 0.5, bunch 1, `1e3` per m); fov in
`{-15, 0, 90, 180, nan, inf}` rejected, `{15, 89.9}` accepted.
