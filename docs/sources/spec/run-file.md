# Run file

A run file is TOML (`.toml`) or JSON (`.json`) with the same keys. `obdeect.sources.config`
checks it and fills the defaults ([python API](python-api.md)). Unknown keys and kinds are
errors, and the message names the key path:

```text
events[0].emissions[1].angular.kind: unknown kind 'cones' (one of isotropic, cone, ...)
```

Required keys have no default. "Table" values are arrays of equal length (at least two entries).

## `[run]`

| Key | Type | Default | Notes |
| --- | --- | --- | --- |
| `number` | int >= 1 | 1 | run number |
| `seed` | int 0 to 2^64 - 1 | 1 | run seed |
| `observation_level_m` | number | required | height of the observation level |
| `atmosphere_id` | int | -1 | -1 = none |
| `frame` | `sky`, `telescope`, `camera`, `ground` | `sky` | frame of all photons. `ground` (G) writes sphere rows: positions relative to the centre of the telescope sphere, exactly one telescope, no `to_plane_z0` (see [run.md](run.md)) |
| `metadata` | table | `{}` | values: int, float, string or list of numbers; added to the run record |

## `[[telescopes]]`

At least one. Distinct ids.

| Key | Type | Default | Notes |
| --- | --- | --- | --- |
| `id` | int >= 1 | position in the list, from 1 | |
| `position_m` | 3 numbers | `[0, 0, 0]` | |
| `radius_m` | number > 0 | required | sphere radius |

## `[[events]]`

At least one.

| Key | Type | Default | Notes |
| --- | --- | --- | --- |
| `repeat` | int >= 0 | 1 | number of events written |
| `metadata` | table | `{}` | added to each event record |

## `[[events.emissions]]`

At least one per event. Keys for every kind:

| Key | Type | Default | Notes |
| --- | --- | --- | --- |
| `kind` | `point`, `line`, `plane_wave`, `sky_background` | required | |
| `photons` | number >= 0 | required | mean physical photons |
| `bunch_size` | number >= 1 | 1 | photons per row |
| `time_offset_ns` | number | 0 | |
| `to_plane_z0` | bool | false | move rows with `z * cz < 0` (or `z = 0`) to z = 0; drop and count the others (`not_on_plane`); not allowed in a `ground` run |
| `spectrum` | spectrum table | required | |
| `pulse` | pulse table | `{ kind = "fixed", time_ns = 0 }` | |
| `target` | `"none"`, `"auto"`, or disc or sphere table | `"none"` | |

Keys per kind:

| `kind` | Key | Type | Default |
| --- | --- | --- | --- |
| `point` | `position_m` | 3 numbers | required |
| | `axis` | 3 numbers | `[0, 0, 1]` |
| | `angular` | angular table | `{ kind = "isotropic" }` |
| `line` | `start_m` | 3 numbers | required |
| | `direction` | 3 numbers | required |
| | `length_m` | number | required |
| | `beta` | number in [0, 1] | 1 |
| | `angular` | angular table | `{ kind = "cone", angle_deg = 0 }` |
| | `profile` | `{ distance_m = [...], density = [...] }` | absent (uniform) |
| `plane_wave` | `towards_source` | 3 numbers (unit; photons travel along -towards_source; z = 0 horizontal) | required |
| | `distance_m` | number > 0 | absent (infinite) |
| | `disc_radius_m` | number > 0 | required |
| | `disc_z_m` | number | 0 |
| | `start_z_m` | number | 0 (on the source side of `disc_z_m`; not for horizontal directions). At most one of `start_z_m` and `start_distance_m` |
| | `start_distance_m` | number >= 0 | absent. Optional alternative to `start_z_m`; required for a horizontal `towards_source` |
| `sky_background` | `radius_m` | number > 0 | required |
| | `centre_m` | 3 numbers | `[0, 0, 0]` |
| | `axis` | 3 numbers | `[0, 0, -1]` |
| | `viewcone_deg` | number in (0, 180] | 90 |

Plane-wave starts (at most one key; giving both is an error; a horizontal `towards_source` without
`start_distance_m` is an error, and the message names it):

- `start_z_m` (default 0 when neither key is given, so existing run files are unchanged):
  photons start on the plane z = `start_z_m`, on the source side of `disc_z_m` (above it for
  `towards_source.z > 0`, below it for `< 0`). Not for horizontal directions.
- `start_distance_m`: photons start on the plane perpendicular to `towards_source` at this distance
  from the disc centre `(0, 0, disc_z_m)`, towards the source. Valid for any direction. It must be at
  least `disc_radius_m * sqrt(1 - towards_source.z^2)` (the plane lies on the source side of the
  whole disc). For a vertical beam, `start_distance_m = start_z_m - disc_z_m` gives the same photons.

Target tables: `{ kind = "disc", centre_m, normal, radius_m }` and
`{ kind = "sphere", centre_m, radius_m }`. All their keys are required; `radius_m > 0`.

## Spectrum kinds

| `kind` | Keys (default) |
| --- | --- |
| `fixed` | `wavelength_nm` |
| `uniform` | `low_nm`, `high_nm` |
| `inverse_square` | `low_nm`, `high_nm` (density ~ 1/lambda^2, Cherenkov) |
| `table` | `wavelength_nm` [...], `density` [...] |

## Pulse kinds

| `kind` | Keys (default) |
| --- | --- |
| `fixed` | `time_ns` |
| `tophat` | `width_ns`, `delay_ns` (0) |
| `gauss` | `sigma_ns`, `mean_ns` (0) |
| `exponential` | `decay_ns`, `mean_ns` (0) |
| `complex` | `flat_ns`, `gauss_rms_ns`, `decay_ns`, `delay_ns` (0) |
| `table` | `time_ns` [...], `rate` [...] |

## Angular kinds

| `kind` | Keys (default) |
| --- | --- |
| `isotropic` | none |
| `cone` | `angle_deg` (fixed angle) |
| `filled_cone` | `angle_deg` |
| `lambertian` | none |
| `gauss` | `sigma_deg` |
| `rayleigh` | none |
| `table` | `angle_deg` [...], `value` [...] |

## Example

Two telescopes, one event with a point source and a line, repeated 10 times:

```toml
[run]
seed = 7
observation_level_m = 2150.0

[[telescopes]]
position_m = [0.0, 0.0, 0.0]
radius_m = 12.5

[[telescopes]]
position_m = [100.0, 0.0, 0.0]
radius_m = 12.5

[[events]]
repeat = 10

[[events.emissions]]
kind = "point"
photons = 1e6
position_m = [0.0, 0.0, 300.0]
spectrum = { kind = "fixed", wavelength_nm = 400.0 }
target = "auto"

[[events.emissions]]
kind = "line"
photons = 1e5
start_m = [-500.0, 0.0, 10.0]
direction = [1.0, 0.0, 0.0]
length_m = 1000.0
angular = { kind = "cone", angle_deg = 1.2 }
spectrum = { kind = "fixed", wavelength_nm = 400.0 }
target = "auto"
```
