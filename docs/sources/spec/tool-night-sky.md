# obdeect sources night-sky (implementation)

User options and examples: [night-sky](../programs/night-sky.md). Source:
[background](background.md). Writer B, producer `obdeect sources night-sky`.

## Validation (at plan time, also with `--photons 0`)

- Radius: finite and `> 0`. Viewcone: `0 < v <= 180` deg. Altitude: finite.
- Atmosphere: integer `>= -1` (`--atmosphere: profile number (or -1 for none)`).
- Photons per event: finite, `0..2^53`. Bunch size: finite, `> 0`; below 1 used as 1.
- Events: integer `0..2^32-1`. Pulse width: finite `>= 0`. Pulse delay: finite.
- Default spectrum: `Benn_LaPalma_sky_converted.lis` must exist in the working directory;
  else `--spectrum: cannot open table Benn_LaPalma_sky_converted.lis`.
- Table spectrum: nm and value columns; repeated wavelengths are merged (mean value).

## Plan

- One case, `events` events; one emission `BackgroundSource{centre (0,0,0), radius, pointing
  (0,0,-1), viewcone = v deg, spectrum, pulse width, pulse delay}`.
- Mean bunches `photons / b`; seed = `--seed`; rows written with `w = b` (float32);
  telescope key 0; count by `background_count` (rounding rule).
- One telescope id 1 at the origin with radius `r`; `wavelength_min/max` from the spectrum.
- Run extra keys: `viewcone_deg` (f64), `background` (string `night sky`). `atmosphere_id` = `-A`.
- `cos(theta)` from the downward axis is uniform in `[cos v, 1]`; photons travel with
  `cz <= -cos v`; all positions lie at distance `r` from the origin.

## Acceptance (`test_night-sky`)

- Count `3e4` for `6e4` photons with bunch 2; directions inside the cone to 1e-6; `|position| = 15`
  to 1e-3; `w = 2`; cos moments to 4 sigma.
- Projected disc: `ri^2` uniform in `[0, 1]` (mean 1/2, KS below 1.95/sqrt(n)); `E[zi] = 2/3`.
- Mean time `= delay - E[zi] r / c` with pulse width 4 and delay 10.
- Uniform 300..600 table spectrum moments; tilted viewing axes always enter the sphere.
- Stream layout: 2 events, `2 x 20000` rows, `telescope_r_m = [15]`; `zem - z = altitude`.
- Negative photons, radius `<= 0`, viewcone `0` or `200`, negative pulse width, NaN delay,
  bunch 0, NaN altitude all throw at plan time; viewcone 180 is accepted.
