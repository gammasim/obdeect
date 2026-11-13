# obdeect sources flasher (implementation)

User options and examples: [flasher](../programs/flasher.md). Writer A, producer
`obdeect sources flasher`, `telescope_id = 1`, telescope at (0, 0, 0).

## Option parsing

- Arguments are collected as (name, value) pairs; `--help` anywhere prints usage and exits 0.
- `--type` is evaluated first (default dish), then the other options left to right.
- An option of the other type: `<option> applies only to --type dish|secondary`.
- Other errors (all prefixed `flasher: `): `--type must be dish or secondary, not '<x>'`;
  `--frame` with secondary: `--frame applies only to --type dish`; `--frame` value other than
  sky or camera; `--xy` needs exactly two values with dish (`needs x,y with --type dish`);
  `--delays needs 4 or 5 values` (missing T5 is 0); `--events takes one value with --type
  secondary`; empty `--events` or `--photons` list; negative event count; non-finite bunch size.
- Bunch size below 1 is used as 1. Positive required: `--distance`, `--camera-radius`,
  `--mirror-camera-distance`. Non-negative required: `--flasher-xy`, `--events`, `--photons`.
- `--xy` and `--flasher-xy` are both accepted for both types, with the meaning of the type.

## Dish

- Emitter at `(x, y, distance)` cm converted to m, pointing `(0, 0, -1)`, fiducial radius
  `R = camera_radius / 100` m, target `fiducial_target(position, R, angular)`.
- Photon levels `k = 0..L-1`: level `k` has `events[min(k, len-1)]` events, one emission
  per event, `Source.photons = photons[k] / b`. Event numbers `1, 2, ...` are consecutive
  across levels.
- Emission `j` (counted over all events) has seed `derive_seed(seed, j)`, bunch weight `b`,
  time offset 0, `to_camera_plane = (frame == camera)`. Random keys: telescope 1, event number.
- Camera frame: writer A step 1 (see [stream-writer](stream-writer.md)). Path `z / -cz`,
  `t` gains `path / c`.

## Secondary

- Frame T (code 2). Run radius 1.20 m (not an option), fiducial target radius 1.2 m.
- Flasher `i = 0..3`: `sx = +1` for `i < 2`, else -1; `sy = +1` for even `i`, else -1.
  Position (cm) `(sx xy, sy xy, mirror_camera_distance + depth)` converted to m.
- Pointing: altitude `-90 deg + inclination`, azimuth `atan2(y, x)` of the flasher,
  `pointing = (cos alt cos az, cos alt sin az, sin alt)`. Inside the 1.2 m sphere with
  default values, so no target.
- Emission counter `j` runs over the whole run, so seeds are `derive_seed(seed, j)`.
- Fire modes (per event `iev`, 0-based within a level):

| `--fire` | emissions |
| --- | --- |
| `all` | flashers 0..3, offsets `delays[i]` |
| `alternate` | flasher `iev mod 4`, offset `delays[i]` |
| `1`..`4` | that flasher, offset `delays[n-1]` |
| `demo` | flashers 0..3 with offsets `delays[i]`, then 0..3 again with offset `delays[4]` |

- For `demo`, if all five delays are exactly 0 they become `0,10,20,30,50`.
- Time: pulse time plus offset (float32 add). No `zem` column (frame T).

## Output metadata

Writer A. `observation_level_m` = `--altitude` or the type default; `atmosphere_id` =
`--atmosphere`; `run_number` = `--run`; `wavelength_*` from the spectrum range.

## Reference values

- Dish defaults (20 events): rows per event about `2.5e6 umax / 2`, with
  `umax = 1 - sqrt(1 - (R/d)^2)`, `R = 1.2`, `d = 16`. Every photon has `x = y = 0`,
  `z = 16`, `zem = 2166`, `w = 1`, `t = 0`, `wavelength = 400`; `1 + cz` is uniform on
  `[0, umax]`.
- Secondary `--photons 1e5`: exactly `4 x 100000` rows, frame 2; first block of 100000 rows at
  `(0.1462, 0.1462, 0.5353)` m.
- `--fire demo --photons 1e3`: 8000 rows; row 0 `t = 0`, row 1000 `t = 10`, row 2000
  `t = 20`, row 3000 `t = 30`, rows 4000..7999 `t = 50`.
