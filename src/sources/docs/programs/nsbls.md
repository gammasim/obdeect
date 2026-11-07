# nsbls

Night-sky background photons on one telescope sphere (`sim_telarray/LightEmission/nsbls.cc`).
Tool `dlamb-nsbls`, programme `cpp/include/dlamb/programs/nsbls.hpp`, core
`cpp/include/dlamb/background.hpp` (new: no BackgroundLightSource existed in dlamb).

## Options (original units)

| Original | Default | dlamb |
| --- | --- | --- |
| `-r`, `--radius` [m] | 15 | `telescope_radius_m` |
| `-v`, `--viewcone` [deg] | 10 | `viewcone_deg` (half opening angle) |
| `--altitude` [m] | 2150 | `altitude_m` |
| `-A`, `--atmosphere` | -1 | `atmosphere`, recorded (see below) |
| `--photons` | 1e6 | `photons` (per event) |
| `-b`, `--bunch-size` | 2 | `bunch_size` (below 1 used as 1) |
| `--events`, `--run` | 1, 1 | `events`, `run_number` |
| `-s`, `--spectrum` | Benn_LaPalma_sky_converted.lis | nm, table file (default: this file in the working directory), or `Cherenkov:a:b` |
| `-w`, `--pulse-width` [ns] | 0 | `pulse_width_ns` |
| `-d`, `--pulse-delay` [ns] | 0 | `pulse_delay_ns` |
| (new) `--seed` | 1 | `seed` |

## Mapping and original behaviour

- Count: `floor(photons / bunch)` plus one with the fractional probability; each row `w = bunch`.
- Direction: uniform in cos(theta) inside the cone about (0,0,-1) (`IsotropicRandomEmission`).
- Emission point: `1 - sqrt(U)` projected radius on the unit hemisphere, rotated with
  `Rotate_y(za)`, `Rotate_z(az)` as in the original, scaled by the telescope radius.
- Time: `delay - zi R / c` plus the width term. c is the vacuum speed (the original uses
  the air speed at the observation level).
- Quirks kept: the width enters as `width / 2.35482 x (U - 0.5)` although the help calls
  it a flat width; the photon count is spread over the viewing cone, although the
  class comment says 4 pi.
- The original's default atmosphere `-1` makes it open a file named `-1` and fail
  unless `--atmosphere` is given. dlamb records the number only.
- The table in the working directory is read by default. Its repeated wavelengths
  are merged into one node with the mean value (the original accepts them; dlamb's
  table interpolation requires strictly increasing nodes).

## Checks

- Directions inside the cone (within 1e-6), |position| = R, mean and variance of cos(theta)
  (4 sigma), mean time against `-E[zi] R / c` (4 sigma), a table spectrum's mean and variance.
- Original (container, `--atmosphere 1 --photons 200000`): wrote zero bunches in
  this run. Cause not found yet, so no comparison.
