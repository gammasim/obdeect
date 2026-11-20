# Run

Header `obdeect/sources/run.hpp`. A run is a list of events; each event is a list of emissions
(see [source core](source-core.md)). `write_run(out, run)` writes the photon stream. The
run file (TOML or JSON) is converted to a `Run` by the Python layer ([run file](run-file.md)).

## Types

- `Telescope{id, position, radius_m}`.
- `Emission{source, bunch_size (default 1), time_offset_ns (0), to_plane_z0 (false)}`.
- `Event{repeat (default 1), emissions, metadata}`.
- `Run{number (1), seed (1), observation_level_m, atmosphere_id (-1), frame (sky), telescopes,
  events, metadata, producer}`.
- `RunSummary{events, rows, weight, expected}`: rows written, sum of weights, expected sum of
  weights, all telescopes.

## Stream layout

```text
stream, run, { event, { telescope, photon chunks*, telescope_end }*, event_end }*, run_end, stream_end
```

- Events are numbered from 1 in order, `repeat` times per event entry.
- Each telescope block holds the photons of all emissions of the event, in emission order.
- A chunk holds at most 65536 rows and has `telescope_id`, `frame` and `first_index` (rows written
  before it in the block).
- Photon columns: x y z cx cy cz t w wavelength; plus `zem` in the sky frame.

## Metadata keys

Values are typed as the stream metadata (`i64`, `f64`, `vector<f64>`, string).

`run` record:

| Key | Type | Value |
| --- | --- | --- |
| `run_number` | i64 | `number` |
| `seed` | i64 | `seed` (as signed 64-bit) |
| `observation_level_m` | f64 | |
| `atmosphere_id` | i64 | |
| `refraction_applied`, `curvature_applied` | i64 | 0 |
| `n_telescopes` | i64 | |
| `telescope_x_m`, `_y_m`, `_z_m`, `_r_m` | vector f64 | one entry per telescope |
| `wavelength_min_nm`, `wavelength_max_nm` | f64 | over all emissions (only if there is one) |

Run `metadata` entries are added to the run record. Writer keys take precedence.
The Python layer also stores the normalized run description as `run_description` (JSON text).

`event` record: `event_number` (i64), then the event `metadata`.

`telescope` record: `telescope_id`, `frame` (i64; sky 0, telescope 2, camera 3), `x_m`, `y_m`,
`z_m`, `r_m` (f64).

`telescope_end` record: `telescope_id` (i64), `photons` (f64, sum of the weights written).

## Seeds and draws

- Emission `k` in event order (counting all events, once per emission entry) uses
  `source_seed(run.seed, k)` (splitmix64). Repeats of an event use the same seed.
- Photon draws are keyed on (emission seed, event number, telescope key 0, photon index).
  All telescopes get the same photons; the telescope key is always 0.
- The count of an emission is `photon_count(source, seed, event, 0)`.

## Bunch, time and plane moves

For each emission, candidates are `photons / bunch_size` and each row has weight `w * bunch_size`.
Then, per row:

1. `to_plane_z0`: move straight down to z = 0 (requires `cz < 0` and `z >= 0`, else error);
   time `+ path / c`.
2. `t += time_offset_ns`.
3. Sky frame: `zem = observation_level_m + z`.

## Auto target

`auto` is resolved per emission from the telescopes (`resolve_auto_target`, `fiducial_sphere`). The
fiducial sphere has centre `C`, the mean telescope position, and radius `R` = the largest
`|position - C| + radius` over the telescopes. Rule: see [source model](../sources.md#targets).
Summary: point within `1.02 R` of `C` none; isotropic point sphere `(C, R)`; fixed-angle point none;
other points disc of radius `R` centred on `C` in the plane z = `C.z` (error if the point is in that
plane). Lines with fixed or isotropic profile get the sphere `(C, R)`; other lines none; plane waves
and sky backgrounds none.

## Telescope cache

Photons of one event are the same for all telescopes. If the run has more than one telescope and
the candidates of the event fit in 1 GiB (`40` bytes per candidate: ten float32 columns), the chunks
are generated once and written for each telescope. Otherwise they are generated again per
telescope. The output is the same.

## Summary values

- `expected` = `expected_weight(run)`: expected sum of weights, all events and telescopes. Disc
  targets use a 2000 x 2000 midpoint grid; line sphere targets use 20000 steps along the line.
- `weight` = sum of the written weights; for a run with targets this approaches `expected`.
