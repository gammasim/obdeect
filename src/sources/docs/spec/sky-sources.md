# Sky-level layer

Header `dlamb/sky_sources.hpp`. Shared by point-source, beam-scan, laser-beam, night-sky.

## Helpers

- `rotate_y(v, a) = (x cos a + z sin a, y, z cos a - x sin a)`; `rotate_z(v, a) = (x cos a - y
  sin a, x sin a + y cos a, z)`.
- `direction_from_altitude(alt, az) = (cos alt cos az, cos alt sin az, sin alt)` (rad).
- `fiducial_super_radius(telescopes) = max over telescopes of (|centre| + radius)`, 0 for none.
- `spectrum_range` as in the stream-writer page.
- `shared_telescope_key = 0`: every sky emission is generated once with telescope key 0; the
  same photons are written for every telescope (telescopes only select them).

## Types

- `SkyTelescope{id=1, centre, radius_m}`.
- `SkyRun{run_number=1 (i64), observation_level_m=2150, atmosphere_id=-1 (i64),
  wavelength_min_nm=0, wavelength_max_nm=0, telescopes, extra (Metadata)}`.
- `SourceEmission{source, seed, bunch=1, time_offset_ns=0}`: `count(event) =
  photon_count(source, seed, event, 0)`; `fill(event, first, rows, chunk)` = `emit(...,
  telescope 0, ...)` then for each row `w = float(w) * float(bunch)` and `t = float(t) +
  float(time_offset_ns)` (float32 arithmetic).
- `SkyPlan<Emission>{run, events_per_case=1, cases}`; `write_plan(out, producer, plan)` writes
  the stream with events numbered 1, 2, ... continuing across cases: each case is written
  `events_per_case` times.

## Stream layout

`run` keys: `run_number` i64, `observation_level_m` f64, `atmosphere_id` i64,
`refraction_applied` i64 = 0, `n_telescopes` i64, `telescope_x_m`/`y_m`/`z_m`/`r_m`
`vector<double>`; `wavelength_min_nm` and `wavelength_max_nm` f64 only if
`wavelength_max_nm > 0`; then all `extra` entries inserted without overwriting existing keys
(extra may add doubles and strings).

Per event (number `n`): `event` {`event_number` i64 = n}. Counts of all emissions of the case
are computed once (`count(n)`); then for each telescope in run order: `telescope` record
(`telescope_id`, `frame` = 0, `x_m`, `y_m`, `z_m`, `r_m`), for each emission in order its rows
in chunks of at most 65536 (`first_index` = rows of earlier emissions + offset in the emission),
`telescope_end` {`telescope_id`, `photons` = sum of float32 `w` of the block}; `event_end` {}.
The telescope id of a chunk is the telescope's id.

An error inside `count` or `fill` (source validation) occurs after the header was written;
partial output is allowed in that case. Errors raised while building the plan write nothing.
