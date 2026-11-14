# night-sky: specification

Options and examples: [night-sky](../programs/night-sky.md). Source:
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
