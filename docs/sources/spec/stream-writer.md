# Stream writers

All tools write a photon stream (`obdeect::stream::Writer(out, producer)`; the writer emits the
`stream` record with key `producer`, and `finish()` writes `stream_end`). Metadata values are
int64, double, `vector<double>` or string, exactly as typed below. Chunks hold at most
`65536` rows; a chunk carries `telescope_id`, `frame`, `first_index` and a column set.

Record order (all writers):

```text
stream, run, { event, telescope, photons*, telescope_end, event_end }*, run_end, stream_end
```

Sky-layer writers write, per event, one `telescope` block per telescope; the program writers
(`flasher`, `led`) and muon-track write exactly one telescope block per event.
`first_index` of a chunk is the number of rows written before it in the same telescope block
and event (photon indices run over all emissions of the block in emission order). An emission
with zero photons writes no chunk. In writers A and B every emission is chunked separately (a chunk never spans two emissions, at most 65536 rows each); writer C buffers across segments and flushes full chunks of 65536 rows.

Record order is exact; the split into chunks is not part of the contract.

## Writer A: program writer (flasher, led)

`write_stream(out, RunSetup, events)`; `RunSetup{producer, run_number, observation_level_m,
atmosphere, telescope_id=1, telescope_position=(0,0,0), telescope_radius_m, frame}`.

`run` keys: `run_number` i64, `observation_level_m` f64, `n_telescopes` i64 = 1,
`telescope_x_m`, `telescope_y_m`, `telescope_z_m`, `telescope_r_m` (each `vector<double>` of
one value), `atmosphere_id` i64, `refraction_applied` i64 = 0, and, if at least one emission
exists, `wavelength_min_nm` / `wavelength_max_nm` f64 = min / max over all emissions of the
spectrum range (Fixed: value both; UniformRange, InverseSquare: low, high; Tabulated: first and
last node).

Per event: `event` {`event_number` i64}; `telescope` {`telescope_id` i64, `frame` i64,
`x_m`, `y_m`, `z_m`, `r_m` f64}; photons; `telescope_end` {`telescope_id` i64, `photons` f64 =
sum of the float32 weights written}; `event_end` {}. Then `run_end` {}.

An emission whose frame differs from the run frame is an error ("all emissions of a stream
must use the stream frame"). For each emission with `count = photon_count(source, seed,
event_number, telescope_id)`: rows `0..count-1` are `emit(...)` with `first_index` set to the
rows already written in the block; then, in this order, on the float32 columns:

1. if `to_camera_plane`: for every row require `cz < 0` and `z >= 0` (else
   `invalid_argument`: photon does not travel towards the camera plane); `path = z / -cz`;
   `x += path cx`, `y += path cy`, `z = 0`, `t += path / c` (float32 results);
2. `w = float(w * bunch_size)`, `t = float(t + time_offset_ns)`;
3. if frame is sky: add column `zem = float(observation_level_m + z)` (not written for T and C).

Columns written: x y z cx cy cz t w wavelength (+ zem in sky).

## Writer B: sky layer (point-source, beam-scan, laser-beam, night-sky)

See [sky-sources](sky-sources.md). Same record order; `run` keys differ and the telescope
list may have several telescopes; every chunk is frame sky with `zem = float(observation_level_m
+ z)`; the `telescope_end` `photons` is the sum of the float32 `w` values of the block.

## Writer C: run writer (muon-track)

`RunInfo{run_number, observation_level_m, wavelength_min_nm, wavelength_max_nm, telescopes}`.

`run` keys: `run_number` i64, `observation_level_m` f64, `wavelength_min_nm`,
`wavelength_max_nm` f64 (always written), `refraction_applied` i64 = 0, `n_telescopes` i64,
`telescope_x_m`/`y_m`/`z_m`/`r_m` vectors. No `atmosphere_id`.

`event` keys: `event_number` i64, `primary_id` i64, `energy_gev` f64, `core_x_m` and `core_y_m`
(`vector<double>` with one value: the impact position). `telescope`: as writer A with the
telescope given (frame sky = 0). Photons are buffered and flushed in chunks of 65536 rows
(`first_index` = rows already flushed in the block), columns x y z cx cy cz t w wavelength zem,
`zem = float(observation_level_m + z)`. `telescope_end.photons` = sum of the double
weights pushed (not float-rounded). The writer reports totals (photon rows, weight sum).
Misuse (push outside a block, nested block) is a logic error (`std::logic_error`).

## Frame code in `telescope` records

Integer obdeect stream code: sky 0, T 2, C 3.
