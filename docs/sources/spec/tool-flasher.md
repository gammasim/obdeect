# obdeect sources flasher

Flasher photon stream. `obdeect sources flasher [--type dish|secondary] [options] [-o FILE]`.
Replaces an earlier pair of flasher programs. Writer A (producer `obdeect sources flasher`),
`telescope_id = 1`, telescope position (0,0,0).

## Options

Arguments are first collected as (name, value) pairs (every option except `-h`/`--help` takes
one value; `--help` anywhere prints usage and exits 0); `--type` is evaluated first (default
dish), then the rest left to right. An option of the other type is an error
`<option> applies only to --type dish|secondary`.

| Option (aliases) | Type | Default | Meaning / rules |
| --- | --- | --- | --- |
| `--type` | both | `dish` | `dish` or `secondary`, else `flasher: --type must be dish or secondary, not '<x>'` |
| `--photons N[,..]` | both | `2.5e6` | photons into 4 pi per flasher, list = photon levels; each finite and `>= 0` |
| `--bunchsize F` (`--bs`, `-b`) | both | `1` | photons per bunch; finite; values below 1 are used as 1 |
| `--spectrum` | both | `400` | wavelength or `Cherenkov:A:B` (cli-conventions) |
| `--lightpulse` (`--pulse`) | both | `Simple:0` | pulse text |
| `--angular-distribution` (`--angular`) | both | `isotropic` | angular text |
| `--run` | both | `1` | run number (i64 in stream) |
| `--altitude` (`-a`) | both | 2150 dish, 2100 secondary | observation level [m] |
| `--atmosphere` (`-A`) | both | `1` | integer, recorded as `atmosphere_id` |
| `--seed` | both | `1` | base seed |
| `-o` | both | stdout | output |
| `--events N[,..]` | dish | `1` | events per photon level, each `>= 0`; a short list repeats its last value (extra values beyond the levels unused); secondary: exactly one value |
| `--frame sky\|camera` | dish | `sky` | `--frame` with secondary: error `--frame applies only to --type dish`; other value: `--frame must be sky or camera` |
| `--xy X,Y` (`--flasher-xy`) | dish | `0,0` | offset [cm]; needs exactly two values (`needs x,y with --type dish`) |
| `--distance D` (`-z`) | dish | `1600` | height above the telescope centre [cm], `> 0` |
| `--camera-radius R` (`-r`) | dish | `120` | fiducial radius [cm], `> 0` |
| `--flasher-xy X` (`--xy`) | secondary | `14.62` | diagonal offset [cm], `>= 0` |
| `--flasher-depth D` (`-z`, `--z`) | secondary | `2.46` | depth behind the camera front [cm] |
| `--flasher-inclination A` (`--inclination`) | secondary | `0` | outward inclination [deg] |
| `--mirror-camera-distance D` (`--mirror-distance`) | secondary | `356.1075 - 305.03801 = 51.06949` | [cm], `> 0` |
| `--delays T1,..,T4[,T5]` | secondary | `0,0,0,0,0` | 4 or 5 values [ns] (else `--delays needs 4 or 5 values`); missing T5 is 0 |
| `--fire` | secondary | `all` | `all`, `alternate`, `demo`, `1`..`4` |

Other validation: `--events` empty list or a negative value; `--photons` empty; non-finite
bunch size; `--type secondary` with `--events` of more than one value
(`--events takes one value with --type secondary`); all errors begin `flasher: `.

## Dish

Emitter at `(x, y, distance)` cm -> m, pointing `(0,0,-1)`, with the parsed angular profile,
fiducial radius `R = camera_radius/100` m, target `fiducial_target(position, R, angular)`.
Bunch `b = max(1, bunchsize)`. Photon levels `k = 0..L-1` (L = number of `--photons` values):
level `k` has `events[min(k, len-1)]` events; one emission per event; `Source.photons =
photons[k] / b`; frame = `--frame`. Event numbers `1, 2, ...` consecutive across levels.
Emission `j` (counting from 0 over all events of all levels) has seed `derive_seed(seed, j)`
and bunch weight `b`, time offset 0, `to_camera_plane = (frame == camera)`. Random keys:
`event` = event number, `telescope` = 1. Run radius `R`; run/telescope frame = `--frame`.

Camera frame: the photons run in a straight line from the emission point to the plane z = 0
(writer A, step 1). Count and directions equal the sky output for the same seed; `t` gains
`path/c`; positions fill a disc of radius `R / sqrt(1 - (R/d)^2)` for a source on the axis at
height `d`. A profile whose photons do not travel towards the plane (for example 4 pi emission
from a source inside the fiducial sphere) is an error at write time.

## Secondary

Frame T (code 2), run radius 1.20 m (not an option), fiducial target rule with radius 1.2 m.
Flasher `i = 0..3` (flasher numbers 1..4): `sx = +1` for `i < 2`, else -1; `sy = +1` for `i`
even, else -1; position (cm) `(sx*xy, sy*xy, mirror_camera_distance + depth)`, converted to
m; pointing altitude `alt = (-90 + inclination)` deg, azimuth `az = atan2(y, x)` (of the
flasher), `pointing = (cos alt cos az, cos alt sin az, sin alt)`. Flashers lie inside the
1.2 m sphere with default values, so no target.

Per event (`events` events, numbered from 1; for each photon level in order, so numbers are
consecutive across levels): fire pattern, one emission per firing with `Source.photons =
photons[level] / b`, time offset `delays[...]`, bunch `b`, seed `derive_seed(seed, j)` with
`j` the running emission counter (counts emissions over the whole run, in order):

| `--fire` | emissions of event `iev` (0-based within level) |
| --- | --- |
| `all` | flashers 0..3, offsets `delays[i]` |
| `alternate` | flasher `iev mod 4`, offset `delays[i]` |
| `1`..`4` | that flasher, offset `delays[n-1]` |
| `demo` | flashers 0..3 with offsets `delays[i]`, then flashers 0..3 again with offset `delays[4]` |

For `demo`, if all five delays are exactly 0 they become `0,10,20,30,50`.

Times: pulse time plus offset (float32 add). No `zem` column (frame T).

## Output metadata

Writer A. `run.observation_level_m` = `--altitude` or the type default; `atmosphere_id` =
`--atmosphere`; `run_number` = `--run`; `wavelength_*` from the spectrum range.

## Examples and reference numbers

- Defaults (dish, 20 events): rows per event ~ `2.5e6 * umax/2`, `umax = 1 - sqrt(1 - (R/d)^2)`
  with `R = 1.2`, `d = 16`; every photon `x = y = 0`, `z = 16`, `zem = 2166`, `w = 1`,
  `t = 0`, `wavelength = 400`; `1 + cz` uniform on `[0, umax]`.
- Secondary `--photons 1e5`: exactly `4 * 100000` rows, frame 2, first block of 100000 rows at
  `(0.1462, 0.1462, 0.5353)` m.
- `--fire demo --photons 1e3`: 8000 rows; row 0 `t=0`, row 1000 `t=10`, 2000 `t=20`, 3000 `t=30`,
  4000..7999 `t=50`.
