# night-sky

Night-sky background photons on one telescope sphere. Tool `obdeect sources night-sky`, programme
`cpp/include/obdeect/sources/programs/night-sky.hpp`, core `cpp/include/obdeect/sources/background.hpp`.
Replaces the earlier nsbls program; options keep their names.

## Options

| Option | Default | obdeect sources field |
| --- | --- | --- |
| `-r`, `--radius` [m] | 15 | `telescope_radius_m` |
| `-v`, `--viewcone` [deg] | 10 | `viewcone_deg` (half opening angle) |
| `--altitude` [m] | 2150 | `altitude_m` |
| `-A`, `--atmosphere` | -1 | `atmosphere`, recorded only |
| `--photons` | 1e6 | `photons` (per event) |
| `-b`, `--bunch-size` | 2 | `bunch_size` (below 1 used as 1) |
| `--events`, `--run` | 1, 1 | `events`, `run_number` |
| `-s`, `--spectrum` | Benn_LaPalma_sky_converted.lis | nm, table file (default: this file in the working directory), or `Cherenkov:a:b` |
| `-w`, `--pulse-width` [ns] | 0 | `pulse_width_ns` |
| `-d`, `--pulse-delay` [ns] | 0 | `pulse_delay_ns` |
| `--seed` | 1 | `seed` |

## Geometry and time

- Count: `floor(photons / bunch)` plus one with the fractional probability; each row `w = bunch`.
- Direction: uniform in cos(theta) inside the cone about (0,0,-1).
- Emission point: `1 - sqrt(U)` projected radius on the unit hemisphere, rotated by the
  view direction (`Rotate_y(za)`, `Rotate_z(az)`) and scaled by the telescope radius.
- Time: `delay - zi R / c` plus the width term, with c the vacuum speed. The width enters
  as `width / 2.35482 x (U - 0.5)`.
- The photon count is spread over the viewing cone.
- The atmosphere number is recorded only.
- A table's repeated wavelengths are merged into one node with the mean value (the table
  interpolation requires strictly increasing nodes).

## Checks

- Directions inside the cone (within 1e-6), |position| = R, mean and variance of cos(theta)
  (4 sigma), mean time against `-E[zi] R / c` (4 sigma), a table spectrum's mean and variance.
- Invalid source parameters are rejected when the plan is built, also with `--photons 0`:
  `--radius -1` and `--viewcone 200` fail with an error.
- A run of the reference program (`--atmosphere 1 --photons 200000`) wrote zero bunches.
  Cause not found yet, so no comparison.
