# Source model

A source emits photons from a geometry. A run file lists the sources (emissions) of each event;
`obdeect sources run` writes the photons as a photon stream (`.pst`).

## Geometries

| `kind` | Photons start | Use for |
| --- | --- | --- |
| `point` | at a point, into a cone or 4 pi | flashers, LEDs, point sources, beam scans, laser pulses |
| `line` | along a segment | muon tracks, air showers |
| `plane_wave` | on a plane, parallel or from a finite distance, from above or below | stars, parallel beams |
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
start_z_m = 100.0
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

- The run sets one frame for all its photons: `sky` (default), `telescope`, `camera` or `ground`.
  Positions and directions are in that frame; `sky` is the default, with the origin at the
  telescope on the observation level and z up. Photons in the sky frame also get `zem`, the
  height above sea level.
- `ground` (G) is the sphere representation of the [stream format](../stream/format.md): positions
  are relative to the centre of the one telescope sphere (z up, sky axes), and the rows are the
  points where each photon's path enters the sphere (or the emission point inside it). See
  [ground runs](spec/run.md#ground-runs-sphere-rows). Photons that miss the sphere are dropped and
  counted (`misses_sphere`).
- Lengths in metres, times in ns, wavelengths in nm, angles in degrees.
- Directions are unit vectors of the propagation direction, and sources always write `cz`:
  `cz < 0` is downward, `cz > 0` upward and `cz = 0` horizontal (see
  [ray tracing](../ray-tracing/light-sources.md) and the [stream format](../stream/format.md)).
  Emission produces all directions. Horizontal photons (`cz = 0`) are written in every frame. In
  the sky, T and C frames they are ordinary rows; a consumer that works in the G frame needs the
  sphere representation (`ground` runs), since plane-mode consumers (propagation,
  `obdeect stream to-simtel`) cannot take horizontal or upward photons.

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
| `{ kind = "disc", centre_m, normal, radius_m }` | the photons that reach a disc (each row a hit, weight 1) |
| `{ kind = "sphere", centre_m, radius_m }` | photons inside a cone around a sphere: isotropic points; isotropic or fixed-angle lines |

Auto rule. The fiducial sphere has centre `C` (the mean telescope position) and radius `R` (the
largest distance from `C` to a telescope centre plus that telescope's radius). Two telescopes at
x = 0 and x = 100 m with radius 12.5 m give `C` = (50, 0, 0) and `R` = 62.5 m.

- `point` inside the sphere (distance to `C` at most `R`): none, all its light is inside.
- `point` outside the sphere, isotropic: sphere `(C, R)`. Every point outside the sphere gets the
  sphere, so the expected light changes continuously with the distance.
- `point` outside, fixed angle: none.
- `point` outside, other profiles: disc of radius `R` centred on `C` in the plane z = `C.z`. Error
  if the point lies in that plane.
- `line` with a fixed angle or isotropic profile: sphere `(C, R)`; other lines: none.
- `plane_wave` and `sky_background`: none.

In a `ground` run the auto target is the sphere `(origin, R)` of the one telescope, since the
positions are relative to its centre. Explicit targets are also relative to that centre. A target
only chooses the directions that are sampled: each row still goes through the sphere rule, so a
sampled photon that misses the telescope sphere is counted in `misses_sphere`.

Targets do not assume a direction. A point below a sphere that emits upward into it gets the same
rule (the cone about the sphere seen from the point); a point above a disc target gets the disc
in the same way. The weights are the solid-angle fractions in both cases.

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
- `to_plane_z0 = true` moves every row along its direction to the plane z = 0. A row reaches the
  plane when `z = 0` (no move) or `z * cz < 0` (the plane lies ahead, from above or from below):
  the path is `-z / cz`, `x` and `y` move along it, and the time grows by path / c. Rows that move
  away from the plane (`z * cz >= 0` with `z != 0`, horizontal rows included) are dropped and
  counted as `not_on_plane`; this is not an error. Isotropic light from above a plane therefore
  keeps about half of its photons. The kept rows are numbered consecutively in emission order, and
  the random draws of the kept photons are the same as without the move.
- `misses_sphere` (ground runs) counts photons whose path misses the telescope sphere, in the same
  way as `not_on_plane`: reported by the CLI, the Python summary and the `run_end` record.
- `not_on_plane` counts the dropped photons of the whole run (all events, counted once for all
  telescopes). It is reported by `obdeect sources run` (summary line), by the Python `run.write`
  summary and in the `run_end` record. It is not an efficiency factor: no efficiency status
  record is written, since the record is written before the photons exist. `expected` still counts
  all emitted candidates, so `expected - weight` is the weight of the dropped photons plus the
  weight lost to targets.

## Random numbers

Each emission has its own seed, derived from the run seed (see [run](spec/run.md)). Draws are
counter-based (Philox), so the output does not depend on chunking or threads
([random numbers](../stream/random-numbers.md)).

## Laser pulses

The `lidar` and `laser-beam` presets emit a point pulse with a filled-cone divergence. The beam
may point down, up or horizontally. A horizontal beam starts `--start-distance` along the beam from
the laser; its photons are cz = 0 rows (sphere representation for G consumers). The scattering chain
is plane-parallel and only for beams within `--max-tilt`. The pulse contains E lambda / (h c)
photons; `--rows` controls sampling. Propagation models scattering and transmission. Use
`--help` for the required geometry and options. The presets reject a direct beam intersecting a
telescope sphere and limit the tilt of non-horizontal beams (see [presets](presets.md#laser-beam)).
