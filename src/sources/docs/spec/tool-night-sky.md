# dlamb-night-sky

Night-sky background photons on one telescope sphere in a cone about the downward viewing axis (0,0,-1), sky frame.
Writer B, producer `dlamb-night-sky`; source in [background](background.md).

## Options

| Option (aliases) | Default | Rules |
| --- | --- | --- |
| `-r`, `--radius`, `--telescope-radius` | 15 | m, finite and `> 0` |
| `-v`, `--viewcone` | 10 | half opening angle [deg], `0 < v <= 180` |
| `-a`, `--altitude`, `--obslev`, `--observation-level` | 2150 | finite |
| `-A`, `--atmosphere` | -1 | integer `>= -1` (`--atmosphere: profile number (or -1 for none)`) |
| `--photons` | `1e6` | per event, finite, `0..2^53` |
| `-b`, `--bs`, `--bunch-size` | 2 | finite `> 0`; below 1 used as 1 |
| `--events` | 1 | integer `0..2^32-1` |
| `--run` | 1 | run if `> 0` else 1 |
| `-s`, `--spectrum` | table file `Benn_LaPalma_sky_converted.lis` in the working directory | number 200..2000, `Cherenkov:a:b`, or table (nm, value) |
| `-w`, `--pulse-width` | 0 | ns, finite `>= 0` |
| `-d`, `--pulse-delay` | 0 | ns, finite |
| `--seed` | 1 | |
| `-o`, `--output` | stdout | |

Without `-s` the default table must exist in the working directory, else the error is
`--spectrum: cannot open table Benn_LaPalma_sky_converted.lis`. The table's content (a
night-sky photon spectrum) is external data. Parameters are validated when the plan is built,
also for `--photons 0`.

## Plan

One case, `events` events; one emission: `BackgroundSource{centre (0,0,0), radius, pointing
(0,0,-1), viewcone = v deg, spectrum, pulse width, pulse delay}`; bunch `b = max(1,
bunchsize)`; mean bunches `photons / b`; seed = `--seed`; rows written with `w = b`
(float32); telescope key 0; count by `background_count` (rounding rule). One telescope id 1 at
the origin with radius `r`; `wavelength_min/max` from the spectrum range. Run extra keys:
`viewcone_deg` (f64), `background` (string `night sky`). `atmosphere_id` = `-A`.

`cos(theta)` (angle from the downward axis (0,0,-1); photons travel with `cz <= -cos v`) is uniform in `[cos v, 1]`; all positions lie at
distance `r` from the origin.

## Acceptance (`test_night-sky`)

Count `3e4` for `6e4` photons with bunch 2; directions inside the cone to 1e-6; `|position| =
15` to 1e-3; `w = 2`; cos moments to 4 sigma; mean time `= delay - E[zi] r / c` with pulse
width 4 and delay 10; uniform 300..600 table spectrum moments; tilted viewing axes always
enter the sphere; stream layout (2 events, `2*20000` rows, `telescope_r_m = [15]`); `zem - z =
altitude`; negative photons, radius `<= 0`, viewcone `0` or `200`, negative pulse width, NaN
delay, bunch 0, NaN altitude all throw at plan time, viewcone 180 accepted.
