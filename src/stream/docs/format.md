# Photon-stream format (version 1.0)

A gstrickd stream is a sequence of records. All integers and floats are
little-endian. Lengths are in metres, times in nanoseconds, wavelengths in
nanometres, altitudes in metres above sea level.

## Record

| Field | Type | Meaning |
|---|---|---|
| magic | u32 | `0x52545347` (`"GSTR"`) |
| type | u16 | record type (table below) |
| flags | u16 | reserved, 0 |
| length | u64 | payload length in bytes |
| payload | bytes | `length` bytes |

Readers skip records of unknown type. A stream starts with one `stream`
record and ends with `stream_end`.

| Type | Name | Payload |
|---|---|---|
| 1 | `stream` | metadata: `format_major` (int), `format_minor` (int), `producer` (string) |
| 2 | `run` | metadata (run keys below) |
| 3 | `event` | metadata (event keys below) |
| 4 | `telescope` | metadata: `telescope_id` (int), `frame` (int), optional `x_m`, `y_m`, `z_m`, `r_m` |
| 5 | `photons` | photon chunk (below) |
| 6 | `telescope_end` | metadata: `telescope_id` (int), `photons` (float, sum of weights) |
| 7 | `event_end` | metadata (optional keys) |
| 8 | `run_end` | metadata (optional keys) |
| 9 | `blob` | u32 `eventio_type`, u32 `eventio_id`, then the raw EventIO block (header + body) |
| 10 | `ledger` | metadata: efficiency ledger (below) |
| 15 | `stream_end` | empty |

Nesting: `stream` (`run` (`ledger`* `blob`* (`event` `blob`* (`telescope`
`photons`* `telescope_end`)* `event_end`)* `run_end`)*) `stream_end`.
A `ledger` record may also follow a `telescope` record and then applies to
that telescope only.

Readers reject a stream whose `format_major` differs from theirs and accept
any `format_minor`.

## Metadata

A metadata payload is a u32 entry count followed by entries:

| Field | Type |
|---|---|
| key length | u16 |
| key | UTF-8 bytes |
| value type | u8: 1 = int (i64), 2 = float (f64), 3 = string, 4 = float array |
| value | i64 / f64 / (u32 length + bytes) / (u32 count + f64 values) |

Keys are unique within one payload. Unknown keys are kept, not dropped.

Run keys (from CORSIKA RUNH and TELPOS, or from dlamb): `run_number`,
`observation_level_m`, `energy_min_gev`, `energy_max_gev`, `spectral_index`,
`n_telescopes`, `telescope_x_m`, `telescope_y_m`, `telescope_z_m`,
`telescope_r_m` (float arrays), `atmosphere_id`, `wavelength_min_nm`,
`wavelength_max_nm`, `refraction_applied` (int 0/1).

Event keys (from CORSIKA EVTH and TELOFF): `event_number`,
`primary_id` (CORSIKA particle code), `energy_gev`, `zenith_rad`,
`azimuth_rad` (CORSIKA convention), `first_interaction_m`, `core_x_m`,
`core_y_m` (float arrays, one per array reuse), `array_index`.

## Photon chunk

| Field | Type | Meaning |
|---|---|---|
| telescope_id | u32 | as in the enclosing `telescope` record |
| frame | u32 | 0 = sky, 1 = ground (G), 2 = telescope (T), 3 = camera (C) |
| columns | u32 | bit mask of present columns (table below) |
| rows | u32 | number of photons, at most 2^16 |
| first_index | u64 | index of the first row among all photons of this telescope and event |
| data | f32 | one contiguous column of `rows` values per set bit, in bit order |

| Bit | Column | Meaning |
|---|---|---|
| 0 | `x` | position x [m] |
| 1 | `y` | position y [m] |
| 2 | `z` | position z [m] (frame sky and T; in G and C the plane is z = 0) |
| 3 | `cx` | direction cosine x |
| 4 | `cy` | direction cosine y |
| 5 | `cz` | direction cosine z (if absent: cz = −sqrt(1 − cx² − cy²)) |
| 6 | `t` | arrival time [ns] |
| 7 | `zem` | emission altitude [m a.s.l.] |
| 8 | `w` | weight (photons in the bunch) |
| 9 | `wavelength` | wavelength [nm]; 0 = not yet assigned |
| 10 | `p` | acceptance probability already applied (for reweighting) |

Columns `x y cx cy t w` are required; the others are optional.

Frames:

- sky: absolute emission position (x north, y west, z up, origin at the array
  centre on the observation level, as CORSIKA) and direction.
- G: horizontal plane through the telescope sphere centre, x and y relative
  to the centre, CORSIKA axes.
- T: telescope-local 3D frame (z along the optical axis towards the sky).
- C: camera-neutral focal-plane frame (x, y on the camera entrance, telescope
  axes, z = 0); conversion to sim_telarray conventions only in `ps2simtel`.

## Efficiency ledger

Metadata of a `ledger` record:

- `factors` (string): comma-separated factor names, e.g.
  `atmosphere,qe,mirror_reflectivity,telescope_transmission`.
- per factor `<name>.state` (string): `applied`, `envelope` (only the
  envelope f^env was applied), or `pending`.
- per factor `<name>.source` (string): provenance (file name, model version).
- per factor with a λ table: `<name>.wavelength_nm` and `<name>.value`
  (float arrays) of the applied (envelope) values.

A module applying factor k with an envelope in the ledger applies
f_k / f_k^env and sets `<name>.state` to `applied`. A module never applies a
factor whose state is `applied`.

## Random numbers

Philox-4x32-10 (`gstrickd/random.hpp`), keyed on (seed, stage id) with
counter (event, telescope, photon index, stream). Results do not depend on
chunk size, photon order or thread count.
