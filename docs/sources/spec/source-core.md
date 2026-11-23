# Source core

User view: [source model](../sources.md). This page lists types, validation, counts and
sampling. Headers `obdeect/sources/source.hpp` (types, validation, counts) and
`obdeect/sources/emit.hpp` (sampling). Namespace `obdeect::sources`.

## Types

- `Vec3{x, y, z}`; `dot`, `cross`, `norm`, `normalised`, `is_finite`; `orthonormal_basis(n)`
  gives `(e1, e2)`, right-handed.
- `Point{position, axis (unit, default +z), angular (default Isotropic)}`.
- `Line{start, direction (unit), length_m, beta (default 1), angular (default FixedAngle 0),
  profile (LinearDensity over distance, empty = uniform)}`.
- `PlaneWave{towards_source (unit; z = 0 horizontal), distance_m (default +inf), disc_radius_m,
  disc_z_m, start_z_m (default 0), start_distance_m (optional)}`: at most one of the two starts
  is given; neither means `start_z_m = 0`.
- `SkyBackground{centre, radius_m, axis (unit, default -z), viewcone_deg (default 90)}`.
- `Geometry = Point | Line | PlaneWave | SkyBackground`.
- Targets: `DiscTarget{centre, normal, radius_m}`, `SphereTarget{centre, radius_m}`,
  `AutoTarget` (resolved by the run, see [run](run.md)); `Target = monostate | DiscTarget |
  SphereTarget | AutoTarget`.
- `Source{geometry, spectrum (default Fixed 400 nm), pulse (default Fixed 0), photons, target}`.
  A source has no frame; the run sets it.
- `Photon{position, direction, time_ns, weight = 1, wavelength_nm}`.
- Columns of emitted rows: `source_columns` = x y z cx cy cz t w wavelength (no zem).

## Validation

`validate(source)` throws `std::invalid_argument`. It rejects an `AutoTarget` (resolve it first).

- `photons` finite, `>= 0`, `<= 2^53`.
- Spectrum: Fixed positive; UniformRange `0 < low <= high`; InverseSquare `0 < low < high`;
  Tabulated non-empty, lower node positive. All within limits (`|x| <= 1e30`).
- Pulse: within limits; UniformRange `low <= high`; Gauss and Exponential `sigma`/`tau >= 0`;
  Complex widths `>= 0`; Tabulated non-empty.
- Point: position within limits; axis unit; angular in `u = 1 - cos(theta)` in `[0, 2]`.
  Disc target: centre and normal ok, radius positive, not a fixed angle (needs a density), point
  not in the disc plane. Sphere target: isotropic only, point outside the sphere.
- Line: start within limits; direction unit; `length_m >= 0`; `beta` in `[0, 1]`; profile within
  `[0, length]`; no disc target; sphere target only with fixed angle or isotropic.
- PlaneWave: no target; `towards_source` unit; distance infinite or positive; disc radius
  positive; at most one start key (both is an error); with neither, `start_z_m = 0`. `start_z_m`:
  not for z = 0 (a horizontal direction needs `start_distance_m`), plane on the source side of the
  disc. `start_distance_m`: finite, `>= 0`, and
  `>= disc_radius * sqrt(1 - towards_source.z^2)` (the plane lies on the source side of the whole
  disc).
- SkyBackground: no target; radius positive; axis unit; viewcone in `(0, 180]` deg.

## Counts

`expected_count(source)` is the mean number of candidate rows. `photon_count(source, seed,
event, telescope)` rounds it (see [conventions](conventions.md#rounding-to-counts)).

| Case | Mean candidates |
| --- | --- |
| Point, none or fixed-angle | `photons` |
| Point, disc | `photons * I_max / (2 pi integral(I du)) * pi R^2 / h^2`, `h` = distance of the point to the disc plane |
| Point, sphere (isotropic) | `photons * (1 - cos_max) / 2`, `cos_max = sqrt(1 - (R/d)^2)` |
| Line, none | `photons` |
| Line, sphere (isotropic or fixed angle) | `photons * r`, `r` the line reference fraction (below) |
| PlaneWave, SkyBackground | `photons` |

Sphere fraction of a point at distance `d`: `(1 - sqrt(1 - (R/d)^2)) / 2` (1 inside the sphere).
Example: a point 10 m from a sphere of radius 2 m: `(1 - sqrt(0.96)) / 2 = 0.0101`.

Line reference fraction `r` (`line_reference_fraction`): for an isotropic line, the largest sphere
fraction of 4 pi along the line (closed form); for a fixed angle, the largest arc fraction `half / pi`
of a uniform scan of 1025 points along the line. Any positive `r` gives the same expected weights.

## Sampling

`sample_photon(source, seed, event, telescope, index, photon)` draws the time (pulse stream) and
wavelength (spectrum stream) first, then the geometry. Draws a profile ignores are skipped; the
result is the same. `emit(source, seed, event, telescope, first, candidates, chunk)` samples
candidates `first .. first + candidates - 1` and returns the rows written (line misses dropped).

Point, no target: direction = profile quantile of `u[0]` about the axis; azimuth `2 pi u[1]`.
Weight 1.

Point, disc target: point `P` on the disc (target stream, `r = R sqrt(u0)`, `phi = 2 pi u1`);
direction towards `P`; weight `I(1 - cos) / I_max * (h / |P - x|)^3`. The sum of weights
is the expected count on the disc.

Point, sphere target: cone about the direction to the centre, `cos = 1 - u0 (1 - cos_max)`,
azimuth `2 pi u1`. Weight 1.

Line, position: `s = profile.sample(u0)` (target stream), or `length * u0` if uniform. Position
`start + direction s`; time `+ s / (beta c)` (none if `beta = 0`).

Line, no target or a non-sphere: direction from the angular profile about the line direction.
Weight 1.

Line, isotropic sphere target: direction in the cone about the direction to the centre, covering
this point's sphere fraction: `cos = 1 - u0 * 2 * fraction(x)`. Weight `fraction(x) / r`.

Line, fixed-angle sphere target: at each `s`, the azimuths of the cone (angle `theta`, about the line)
that reach the sphere form an arc `centre +- half` (none: the row is dropped). The azimuth is
`centre + half (2 u0 - 1)`; weight `(half / pi) / r`.

PlaneWave, infinite distance, `start_z`: target point `T` on the disc (target stream); direction
`-towards_source`; start where the ray crosses `z = start_z`. Time `+ dot(start - (0,0,disc_z),
direction) / c`. Weight 1.

PlaneWave, finite distance, `start_z`: source at `(0,0,disc_z) + towards_source * distance`;
direction from the source to `T`; start at the source if it lies on the target side of the start
plane (`direction.z * (source.z - start_z) >= 0`), else at the crossing of `start_z`. Time
`+ (|start - source| - distance) / c`.

PlaneWave with `start_distance = L`: the start plane is `(X - C) . towards_source = L`, `C =
(0,0,disc_z)`. Infinite distance: `start = T + towards_source (L - (T - C) . towards_source)`; the
time is `-L / c` (direction `-towards_source`). Finite distance: the source is at depth `distance`
along `towards_source`; the photons start at the source if `distance <= L`, else where their ray
crosses the plane (`start = source + direction (L - distance) / (direction . towards_source)`).
The time rule is that of the finite case above. For a vertical beam, `L = start_z - disc_z` gives
the same photons as `start_z`.

The `start_z` plane lies on the source side of the disc: `(start_z - disc_z) * towards_source.z >=
0` (above the disc for a beam from above, below it for a beam from below). Horizontal
`towards_source` (z = 0) needs `start_distance`.

SkyBackground: direction in the viewcone about the axis, `cos = 1 - u0 (1 - cos_cone)`. The start
is uniform over the disc seen along the direction (target stream), lifted onto the sphere of
radius `R` about the centre. Time `- h R / c` with `h` the height of the start point above
the centre plane.

Output rows are float32; computation is in double.

## Frame

The source sampler does not know the frame; the run writes the chunk and adds `zem` in the sky frame
(see [run](run.md)).
