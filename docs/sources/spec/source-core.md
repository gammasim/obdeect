# Source core

User view: [sources](../sources.md). This page lists types, validation and sampling.

Headers `obdeect/sources/source.hpp`, `streams.hpp`, `vec3.hpp`. Namespace `obdeect sources`.

## Types

- `Vec3{x,y,z}` double; `dot`, `cross`, `norm`, `normalised`, `is_finite`; `orthonormal_basis(n)`
  returns `(e1, e2)`: `e1 = normalise(n x (1,0,0))` if `|n.x| < 0.9` else `normalise(n x
  (0,1,0))`, `e2 = n x e1` (right-handed). The basis choice fixes the azimuth zero point; it is
  part of the contract only for tests that compare azimuths with this definition.
- `Emitter{position, pointing (unit), angular}` (default pointing (0,0,+1), Isotropic).
- `Star{towards_source (unit, z>0), distance_m (default +inf), disk_radius_m, disk_z_m, start_z_m}`.
- Targets: `DiscTarget{centre, normal (unit), radius_m}`, `SphereTarget{centre, radius_m}`,
  `Target = monostate | DiscTarget | SphereTarget`.
- `Source{frame (default sky), geometry (Emitter|Star), spectrum (Fixed 400), pulse (Fixed 0),
  photons, target}`.
- `Photon{position, direction, time_ns, weight=1, wavelength_nm}`.
- `Stream` enum 1..5 (see conventions), `stage_sources`, `photon_stream_format_major = 1`,
  `source_columns` = x y z cx cy cz t w wavelength (obdeect stream bit mask).

## Validation (`validate(source)`, `invalid_argument`)

- `photons` finite, `>= 0`, `<= 2^53`.
- Frame is sky, T or C (not ground, not any other code).
- Spectrum: Fixed `value` positive-within-limits; UniformRange `0 < low <= high` both
  positive-within-limits; InverseSquare `low < high` strictly, both positive-within-limits;
  Tabulated non-empty with `lower` positive-within-limits and `upper` within limits.
- Pulse: Fixed within limits; UniformRange `low <= high` within limits; Gaussian mean and
  `sigma >= 0` finite within limits; Exponential mean and `tau >= 0`; ComplexPulse all three
  widths `>= 0` and finite within limits and `delay` within limits (an infinite width is an
  error); Tabulated non-empty, `lower`/`upper` within limits.
- Star: no target; `towards_source` finite, `|norm - 1| <= 1e-9`, `z > 0`; distance `+inf`
  (plane wave) or positive-within-limits (so 1e200 is an error); `disk_radius_m`
  positive-within-limits; `disk_z_m`, `start_z_m` within limits.
- Emitter: `position` within limits; `pointing` finite, unit to 1e-9; angular: FixedAngle
  `0 <= u0 <= 2`; FilledCone `0 < u0 <= 2`; TabulatedAngular non-empty with `lower >= 0` and
  `upper <= 2 + 1e-12`.
- DiscTarget: centre within limits, normal unit to 1e-9, radius positive-within-limits; the
  angular profile must not be FixedAngle; `|(centre - position) . normal| > 0` (emitter not in
  the disc plane).
- SphereTarget: centre within limits, radius positive-within-limits; angular must be
  Isotropic; `|centre - position| > radius` (emitter outside).

## Counts

`expected_count(source)` (validates first) =
- star: `photons` (mean photons on the entrance disc);
- emitter without target: `photons`;
- disc target: `photons * candidates_per_photon`, `candidates_per_photon = I_max /
  (2 pi * integral_u) * pi R^2 / h^2`, `h = |(centre - position) . normal|`;
  example: isotropic (I_max 1, integral 2), `R = 2`, `h = 10`: `R^2/(4h^2) = 1/100`;
- sphere target: `photons * (1 - cos_max) / 2`, `cos_max = sqrt(1 - (R/d)^2)`, `d = |centre -
  position|`.

Throws if the result is not finite or exceeds `2^53` (message mentions 2^53; a disc at
`h = 1e-170` overflows). `photon_count(source, seed, event, telescope)` rounds the mean by the
rule in conventions (stream count, photon index 0).

## Emission of photon `i` (`sample_photon`, `emit`)

Draws: pulse stream -> time `pulse.sample`; wavelength stream -> `spectrum.sample`; then
geometry draws. `emit(source, seed, event, telescope, first_index, rows, chunk)` validates,
rejects `rows > 65536` (`invalid_argument`), sets `chunk.telescope_id`, `frame`, `first_index`,
columns `source_columns`, and converts to float32.

### Emitter without target

`u_angle = quantile(direction_u[0])` (full range, `umax = 2`); azimuth `phi = 2 pi
direction_u[1]`; `direction = pointing cos_theta + (e1 cos phi + e2 sin phi) sin_theta`, `cos_theta
= 1 - u_angle`, `(e1,e2) = orthonormal_basis(pointing)`; position = emitter position; weight 1.
FixedAngle gives exactly the angle (`dot(direction, pointing) = cos(angle)` to 1e-6).

### Sphere target (isotropic only)

Cone about `n = normalise(centre - position)` with `cos_max` as above; `cos_theta = 1 -
direction_u[0] (1 - cos_max)`, azimuth `2 pi direction_u[1]`. Weight 1. Every photon passes
within `R` of the centre.

### Disc target

Target stream: `r = R sqrt(u0)`, `phi = 2 pi u1`, point `P = centre + e1 r cos phi + e2 r sin
phi` with `(e1,e2) = orthonormal_basis(normal)`. `direction = (P - position)/|P - position|`,
`cos = clamp(dot(direction, pointing), -1, 1)`, weight
`I(1 - cos) / I_max * (h/|P - position|)^3` (zero outside the profile). Sum of weights is the
expected number of photons reaching the disc. Test: isotropic, `h = 10`, `R = 2`: mean weight
`= (solid angle / 4 pi) * 100`, solid angle `2 pi (1 - 10/sqrt(104))`; a 5 deg FilledCone gives
`w > 0` inside and `w == 0` outside (margin 1e-3 deg).

### Star (plane wave, `distance_m = +inf`)

Target stream: `r = R_disk sqrt(u0)`, `phi = 2 pi u1`; `target = (r cos phi, r sin phi,
disk_z)`; `direction = -towards_source`. The ray starts where it crosses the plane
`z = start_z`: `start = target + direction * ((start_z - target.z)/direction.z)`. Time added:
`dot(start - centre, direction)/c` with `centre = (0,0,disk_z)` (wavefront through the disc
centre at t = 0), on top of the pulse time. Weight 1.

### Star (point source, finite distance)

`position = centre + towards_source * distance`. `direction = normalise(target - position)`.
If `position.z <= start_z` the photon starts at `position`; otherwise at the crossing of the
start plane: `target + direction * ((start_z - target.z)/direction.z)`. Time added:
`(|start - position| - distance)/c`. Test: source at 20 m, direction (0.6, 0, 0.8), disc centre
height 1.5, start plane 30: start at (12, 0, 17.5), `t = -20/c`.

The frame is stored unchanged in the chunk. A plane-wave star or a point source above the start plane starts at `z = start_z` (in its frame, e.g. T); a point source at or below the start plane starts at its own position.
