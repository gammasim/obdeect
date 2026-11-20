# Source model

A source emits photons from a geometry. A run file lists the sources (emissions) of each event;
`obdeect sources run` writes the photons as a photon stream (`.pst`).

## Geometries

| `kind` | Photons start | Use for |
| --- | --- | --- |
| `point` | at a point, into a cone or 4 pi | flashers, LEDs, point sources, beam scans, laser pulses |
| `line` | along a segment | muon tracks, air showers |
| `plane_wave` | on a plane, parallel or from a finite distance | stars, parallel beams |
| `sky_background` | on a sphere, travelling inwards | night-sky background |

Each emission also has a spectrum, a pulse shape, a photon count and a target (see below).
Angular profiles give directions: `isotropic`, `cone`, `filled_cone`, `lambertian`, `gauss`,
`rayleigh` or `table`. Full list: [run file](spec/run-file.md).

## Example

All examples are emissions of this run (header once):

```toml
[run]
observation_level_m = 2150.0

[[telescopes]]
radius_m = 12.5

[[events]]
```

One emission per geometry:

```toml
# point: 1e6 photons from 300 m above the telescope, into 4 pi; the auto target keeps the
# photons that reach the telescope sphere
[[events.emissions]]
kind = "point"
photons = 1e6
position_m = [0.0, 0.0, 300.0]
spectrum = { kind = "fixed", wavelength_nm = 400.0 }
target = "auto"

# line: a 1 km track at 1.2 deg Cherenkov angle, time running with the track (beta = 1)
[[events.emissions]]
kind = "line"
photons = 1e5
start_m = [-500.0, 0.0, 10.0]
direction = [1.0, 0.0, 0.0]
length_m = 1000.0
angular = { kind = "cone", angle_deg = 1.2 }
spectrum = { kind = "fixed", wavelength_nm = 400.0 }
target = "auto"

# plane_wave: parallel light from one direction, through a disc of radius 20 m
[[events.emissions]]
kind = "plane_wave"
photons = 1e6
towards_source = [0.0, 0.6, 0.8]
disc_radius_m = 20.0
spectrum = { kind = "uniform", low_nm = 300.0, high_nm = 600.0 }

# sky_background: diffuse light on a sphere of radius 3 km, from the sky above
[[events.emissions]]
kind = "sky_background"
photons = 1e7
radius_m = 3000.0
viewcone_deg = 30.0
spectrum = { kind = "uniform", low_nm = 300.0, high_nm = 600.0 }
```

Run it with `obdeect sources run run.toml -o out.pst`. Keys with defaults can be left out; the
normalized form shows them: `obdeect sources check run.toml`.

## Angular, spectrum and pulse profiles

```toml
angular  = { kind = "filled_cone", angle_deg = 11.0 }
spectrum = { kind = "inverse_square", low_nm = 300.0, high_nm = 600.0 }  # Cherenkov
pulse    = { kind = "gauss", sigma_ns = 2.5 }
```

Each kind and its keys: [run file](spec/run-file.md).

## Frames and units

- The run sets one frame for all its photons: `sky` (default), `telescope` or `camera`.
  Positions and directions are in that frame; `sky` is the default, with the origin at the
  telescope on the observation level and z up. Photons in the sky frame also get `zem`, the
  height above sea level.
- Lengths in metres, times in ns, wavelengths in nm, angles in degrees.
- Directions are unit vectors of the propagation direction; `cz < 0` is downward
  (see [ray tracing](../ray-tracing/light-sources.md) and the [stream format](../stream/format.md)).

## Photon count and weights

`photons` is the mean number of physical photons of the emission.

- Without a target, every photon is written with weight 1.
- With a target, only photons that can reach it are written. Each has a weight: the sum of
  weights is the expected number of photons on the target.
- `bunch_size = b` writes one row per `b` photons, with weight `b`.

The actual count is `floor(mean)` plus one with the probability of the fraction (not a Poisson draw).

## Targets

| `target` | Meaning |
| --- | --- |
| `"none"` (default) | all photons |
| `"auto"` | chosen from the telescopes (rule below) |
| `{ kind = "disc", centre_m, normal, radius_m }` | photons that reach a disc (weighted) |
| `{ kind = "sphere", centre_m, radius_m }` | photons inside a cone around a sphere: isotropic points; isotropic or fixed-angle lines |

Auto rule. The fiducial sphere has centre `C` (the mean telescope position) and radius `R` (the
largest distance from `C` to a telescope centre plus that telescope's radius). Two telescopes at
x = 0 and x = 100 m with radius 12.5 m give `C` = (50, 0, 0) and `R` = 62.5 m.

- `point` within `1.02 R` of `C`: none.
- `point`, isotropic: sphere `(C, R)`.
- `point`, fixed angle: none.
- `point`, other profiles: disc of radius `R` centred on `C` in the plane z = `C.z`. Error if the
  point lies in that plane.
- `line` with a fixed angle or isotropic profile: sphere `(C, R)`; other lines: none.
- `plane_wave` and `sky_background`: none.

### Line photons with a sphere target

A line with a sphere target gets `photons` x `r` candidates. `r` is the largest fraction of the
sphere seen from a point of the line: closed form for an isotropic line, and for a fixed angle the
largest arc fraction over a 1025-point scan of the line. A candidate has weight `fraction / r`.
For a fixed angle, a candidate whose cone misses the sphere is dropped. The output has fewer rows
than candidates; in the run example above the 1e5-photon line has an expected weight of 3.4e4.

## Time

- `pulse` gives the emission time at the source (a point) or along the line before the
  per-point delay.
- A line adds `s / (beta c)` at distance `s` from its start; `beta = 0` puts all points at the pulse time.
- `time_offset_ns` is added to every row.
- `to_plane_z0 = true` moves every row straight down to z = 0 (error if a photon does not travel down).

## Random numbers

Each emission has its own seed, derived from the run seed (see [run](spec/run.md)). Draws are
counter-based (Philox), so the output does not depend on chunking or threads
([random numbers](../stream/random-numbers.md)).

## Downward laser pulses

The `lidar` and `laser-beam` presets emit a point pulse with a filled-cone divergence.
The pulse contains E lambda / (h c) photons; `--rows` controls sampling. Propagation
models scattering and transmission. Use `--help` for the required geometry and options.
The presets reject a direct beam intersecting a telescope sphere and limit the tilt for
the plane-parallel atmosphere.
