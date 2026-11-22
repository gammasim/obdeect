# Photon-stream format

A photon stream is a sequence of records. Integers and floats are little-endian. Lengths are in
metres, times in nanoseconds, wavelengths in nanometres and altitudes in metres above sea level.

## Record

| Field | Type | Meaning |
|---|---|---|
| magic | u32 | `0x52545347` (`"GSTR"`) |
| type | u16 | record type (table below) |
| flags | u16 | reserved, 0 |
| length | u64 | payload length in bytes |
| payload | bytes | `length` bytes |

Readers skip records of unknown type. A stream starts with one `stream` record and ends with
`stream_end`.

| Type | Name | Payload |
|---|---|---|
| 1 | `stream` | metadata: `format_major` (int), `format_minor` (int), `producer` (string) |
| 2 | `run` | metadata (run keys below) |
| 3 | `event` | metadata (event keys below) |
| 4 | `telescope` | metadata: `telescope_id` (int), `frame` (int), optional `array_index` (int, CORSIKA array reuse), `x_m`, `y_m`, `z_m`, `r_m` |
| 5 | `photons` | photon chunk (below) |
| 6 | `telescope_end` | metadata: `telescope_id` (int), `photons` (float, sum of weights; from CORSIKA: the photon sum of its block header) |
| 7 | `event_end` | metadata (optional keys) |
| 8 | `run_end` | metadata (optional keys) |
| 9 | `blob` | u32 `eventio_type`, u32 `eventio_id`, then the raw EventIO block (header + body) |
| 10 | `efficiency_status` | metadata: efficiency status (below) |
| 15 | `stream_end` | empty |

Nesting: `stream` (`run` (`efficiency_status`* (`event` (`telescope` `photons`*
`telescope_end`)* `event_end`)* `run_end`)*) `stream_end`. `blob` records may appear anywhere
between `run` and `run_end`. An `efficiency_status` record after a `telescope` record applies to
that telescope only.

Blobs keep the input order (CORSIKA: RUNH and TELPOS after `run`; EVTH and TELOFF before the
telescopes; EVTE after them; RUNE before `run_end`). Writers of EventIO output emit them in that
order. CORSIKA split-mode markers (1213, 1214) are blobs between the telescopes of an event;
`obdeect stream to-simtel` does not write them.

The format has a major and a minor version (`obdeect stream --version` prints them). Readers
reject a stream whose `format_major` differs from theirs. Minor versions are compatible within a
major version: a reader also reads a later `format_minor` (the current one is 2), and a column it
does not know fails the column-mask check. A negative `format_minor` is rejected. A stream without
`format_minor` is read as 1.0.

Version history:

- 1.0: photons move downwards; the `cz` column, if present, is negative.
- 1.1: the `cz` column may be positive (upward) or zero (horizontal) in every frame. Without a
  `cz` column a photon moves downwards. Streams without upward photons are the same in 1.0 and 1.1.
- 1.2: the `path` column (distance along the photon path, any frame) and the run key
  `ground_representation` with sphere rows in the G frame. Readers of 1.0 and 1.1 streams reject
  a `path` column and a `sphere` run in them; streams without these are the same in all versions.

## Metadata

A metadata payload is a u32 entry count followed by entries:

| Field | Type |
|---|---|
| key length | u16 |
| key | UTF-8 bytes |
| value type | u8: 1 = int (i64), 2 = float (f64), 3 = string, 4 = float array |
| value | i64 / f64 / (u32 length + bytes) / (u32 count + f64 values) |

Keys are unique within one payload. Unknown keys are kept.

Run keys (from CORSIKA RUNH and TELPOS, or from obdeect sources): `run_number`,
`observation_level_m`, `energy_min_gev`, `energy_max_gev`, `spectral_index`,
`n_telescopes`, `telescope_x_m`, `telescope_y_m`, `telescope_z_m`,
`telescope_r_m` (float arrays), `atmosphere_id`, `wavelength_min_nm`,
`wavelength_max_nm`, `refraction_applied` (int 0/1), `curvature_applied` (int 0/1).
`curvature_applied` is 1 if the producer already applied Earth curvature (CORSIKA CURVED
option, EVTH word 79 = 1 or 2); sources write 0.

Run key `ground_representation` (string, format 1.2): `plane` or `sphere`. Absent means `plane`;
any other value is rejected by readers and writers. It sets how the G-frame photons of the run
are represented (see "Ground frame" below). Producers write `plane` for CORSIKA runs.

Event keys (from CORSIKA EVTH and TELOFF): `event_number`,
`primary_id` (CORSIKA particle code), `energy_gev`, `zenith_rad`,
`azimuth_rad` (CORSIKA convention), `first_interaction_m`, `core_x_m`,
`core_y_m`, `time_offset_ns` (float arrays, one entry per array reuse). One
`event` record per shower; telescope records carry the `array_index`.

## Photon chunk

| Field | Type | Meaning |
|---|---|---|
| telescope_id | u32 | as in the enclosing `telescope` record |
| frame | u32 | 0 = sky, 1 = ground (G), 2 = telescope (T), 3 = camera (C) |
| columns | u32 | bit mask of present columns (table below) |
| rows | u32 | number of photons, at most 2^16 |
| first_index | u64 | index of the first row among all photons of this telescope and event; continues across array reuses (same `telescope_id`, different `array_index`) |
| data | f32 | one contiguous column of `rows` values per set bit, in bit order |

| Bit | Column | Meaning |
|---|---|---|
| 0 | `x` | position x [m] |
| 1 | `y` | position y [m] |
| 2 | `z` | position z [m]: frames sky and T, and G of a sphere run; in G of a plane run and in C the plane is z = 0 |
| 3 | `cx` | direction cosine x |
| 4 | `cy` | direction cosine y |
| 5 | `cz` | direction cosine z: negative downwards, positive upwards (format 1.1), zero horizontal (format 1.1); if absent: cz = −sqrt(1 − cx² − cy²), i.e. downwards |
| 6 | `t` | arrival time [ns] |
| 7 | `zem` | emission altitude [m a.s.l.] |
| 8 | `w` | weight (photons in the bunch) |
| 9 | `wavelength` | wavelength [nm]; 0 = not yet assigned |
| 10 | `p` | acceptance probability already applied (for reweighting) |
| 11 | `path` | distance [m] along the photon's path from its emission point to the row position (format 1.2); finite and >= 0; valid in any frame; written by producers that know it |

Columns `x y cx cy t w` are required; the others are optional. A chunk with 0 rows has column
mask 0. The column mask has 12 bits; bits 12 to 31 are undefined and rejected.

Direction `(cx, cy, cz)` is a unit vector. Without a `cz` column the photon moves downwards. A
photon moving upwards or horizontally needs the `cz` column (format 1.1). Consumers that cannot
handle upward photons must check the sign of `cz` and reject the photon rather than assume it is
downward.

Frames:

- sky: absolute emission position (x north, y west, z up, origin at the array centre on the
  observation level, as CORSIKA) and direction. z is the Cartesian height above the plane through
  the origin. With Earth curvature (run key `curvature_applied`, [](../propagation/curvature.md))
  the Earth centre lies on the vertical through the origin at z = −(R + observation level).
- G, plane (`ground_representation` absent or `plane`): horizontal plane through the telescope
  sphere centre; rows at z = 0 (no `z` column); x and y relative to the centre, CORSIKA axes.
- G, sphere (`ground_representation` `sphere`): 3D rows, relative to the telescope sphere centre
  of radius R (`telescope_r_m` of the run). Each row is the point where the photon's path inside
  the telescope sphere begins, with its direction and the time there: on the sphere (|r| = R) if
  the photon was emitted outside it, the emission point itself (|r| < R) if it was emitted inside.
  A photon emitted outside and crossing the sphere gives one row; a muon track crossing the sphere
  gives rows of both kinds in one event. Columns `x y z cx cy cz t zem w` are required (`z`, `cz`
  included), `path` is optional, and `t` is the time at the row point. Rows satisfy
  `|r| <= R (1 + 1e-6) + 1e-4 m`. Consumers check them with `validate_sphere_rows`
  (C++, `obdeect/stream/photons.hpp`).
- T: telescope-local 3D frame (z along the optical axis towards the sky).
- C: camera-neutral focal-plane frame (x, y on the camera entrance, telescope axes, z = 0).
  Conversion to sim_telarray conventions happens only in `obdeect stream to-simtel`.

Sphere rows are not plane rows: a consumer that reads G rows must check the run's
`ground_representation` and either handle `sphere` runs or reject them. `obdeect stream to-simtel`
converts downward sphere rows on the sphere to plane rows (docs/stream/command-line.md) and rejects
the others; `obdeect stream to-grisu` rejects sphere runs.

The efficiency factors that modules apply are recorded in `efficiency_status` records: [](efficiency-status.md). Random numbers: [](random-numbers.md). Reading and writing from code: [](programming.md).
