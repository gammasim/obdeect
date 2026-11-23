# Run

Header `obdeect/sources/run.hpp`. A run is a list of events; each event is a list of emissions
(see [source core](source-core.md)). `write_run(out, run)` writes the photon stream. The
run file (TOML or JSON) is converted to a `Run` by the Python layer ([run file](run-file.md)).

## Types

- `Telescope{id, position, radius_m}`.
- `Emission{source, bunch_size (default 1), time_offset_ns (0), to_plane_z0 (false)}`.
- `Event{repeat (default 1), emissions, metadata}`.
- `Run{number (1), seed (1), observation_level_m, atmosphere_id (-1), frame (sky), telescopes,
  events, metadata, producer}`. `frame` is sky, telescope, camera or ground (G, sphere rows).
- `RunSummary{events, rows, weight, expected, not_on_plane, misses_sphere}`: rows written, sum of
  weights, expected sum of weights, all telescopes; photons dropped by `to_plane_z0` (once per
  event); photons of a ground run that miss the telescope sphere (once per event).

## Stream layout

```text
stream, run, { event, { telescope, photon chunks*, telescope_end }*, event_end }*, run_end, stream_end
```

- Events are numbered from 1 in order, `repeat` times per event entry.
- Each telescope block holds the photons of all emissions of the event, in emission order.
- A chunk holds at most 65536 rows and has `telescope_id`, `frame` and `first_index` (rows written
  before it in the block).
- Photon columns: x y z cx cy cz t w wavelength; plus `zem` in the sky and ground frames. Ground
  (sphere) rows also have `path`.

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
| `ground_representation` | string | `sphere` in ground runs only (absent: sky, telescope, camera) |

Run `metadata` entries are added to the run record. Writer keys take precedence.
The Python layer also stores the normalized run description as `run_description` (JSON text).

`event` record: `event_number` (i64), then the event `metadata`.

`telescope` record: `telescope_id`, `frame` (i64; sky 0, ground 1, telescope 2, camera 3), `x_m`, `y_m`,
`z_m`, `r_m` (f64).

`telescope_end` record: `telescope_id` (i64), `photons` (f64, sum of the weights written).

`run_end` record: `not_on_plane` (i64), the photons dropped by `to_plane_z0` in the run (see
[bunch, time and plane moves](#bunch-time-and-plane-moves)); `misses_sphere` (i64), the photons
dropped by a ground run because their path misses the telescope sphere (see
[ground runs](#ground-runs-sphere-rows)).

## Seeds and draws

- Emission `k` in event order (counting all events, once per emission entry) uses
  `source_seed(run.seed, k)` (splitmix64). Repeats of an event use the same seed.
- Photon draws are keyed on (emission seed, event number, telescope key 0, photon index).
  All telescopes get the same photons; the telescope key is always 0.
- The count of an emission is `photon_count(source, seed, event, 0)`.

## Bunch, time and plane moves

For each emission, candidates are `photons / bunch_size` and each row has weight `w * bunch_size`.
Then, per row:

1. `to_plane_z0`: a row with `z = 0` stays; a row with `z * cz < 0` moves along its direction to
   z = 0 (`path = -z / cz`, time `+ path / c`). Other rows (`z * cz >= 0`, `z != 0`, horizontal
   included) are dropped and counted in `not_on_plane`. Kept rows keep their order and are
   numbered consecutively (`first_index`); the draws of the kept photons are unchanged.
2. `t += time_offset_ns`.
3. Sky frame: `zem = observation_level_m + z`. Ground runs: see below.

## Ground runs (sphere rows)

A run with `frame = ground` writes the sphere representation of docs/stream/format.md (run key
`ground_representation = "sphere"`). Requirements: exactly one telescope (its sphere is the
origin of the frame) and no `to_plane_z0` emission (error otherwise).

- Source positions and directions are relative to the telescope sphere centre (x, y, z up, sky-frame
  axes). The sphere has radius `R` (the telescope's `radius_m`).
- Per row, with the emission point `P` and direction `d`:
  - `|P| <= R`: the row is `P`, `path = 0`, time unchanged (the photon starts inside the sphere).
  - `|P| > R`: the line `P + s d` meets the sphere at the smallest `s > 0` (roots of
    `s^2 + 2 (P.d) s + |P|^2 - R^2 = 0`). If it hits, the row is the entry point `P + s d` with
    `path = s` and `t += s / c` (vacuum: no atmosphere at this level). If it misses, the photon is
    dropped and counted in `misses_sphere`.
- `zem` is the emission altitude a.s.l.: `observation_level_m + telescope z + P.z`, also for entry
  points (the altitude of the emission, not of the row).
- Columns: x y z cx cy cz t zem w path (+ wavelength). Rows satisfy `validate_sphere_rows`.
- Numbering and draws are as for other runs: the kept rows are numbered consecutively, the draws
  are keyed on the photon index, so the rows do not depend on the chunk size.
- Targets: the sampling target is in the frame of the source (relative to the sphere centre). An
  `auto` target is the sphere `(origin, R)`. A target only chooses the directions that are sampled;
  each row still goes through the sphere rule, so a sampled photon that misses the telescope sphere
  is counted in `misses_sphere`. `expected` includes these photons, `weight` does not.
- Photons written by sky, telescope or camera runs with directions `cz = 0` are unchanged; the
  sphere representation is what horizontal photons need in G.

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
  Dropped photons (`not_on_plane`, `misses_sphere`) are in `expected` but not in `weight`.
