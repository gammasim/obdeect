# night-sky

Night-sky background photons on one telescope sphere, in a cone about the downward axis
(0, 0, -1), sky frame. Source details: [background](../spec/background.md).

```console
$ obdeect sources night-sky --photons 1e6 --viewcone 10 -o night-sky.ps
```

Replaces the earlier nsbls program; options keep their names.

## Options

| Option | Default | obdeect sources field |
| --- | --- | --- |
| `-r`, `--radius` [m] | 15 | `telescope_radius_m` |
| `-v`, `--viewcone` [deg] | 10 | `viewcone_deg` (half opening angle) |
| `--altitude` [m] | 2150 | `altitude_m` |
| `-A`, `--atmosphere` | -1 | `atmosphere` (recorded only) |
| `--photons` | 1e6 | `photons` (per event) |
| `-b`, `--bunch-size` | 2 | `bunch_size` (below 1 used as 1) |
| `--events`, `--run` | 1, 1 | `events`, `run_number` |
| `-s`, `--spectrum` | `Benn_LaPalma_sky_converted.lis` | nm, table file (default: this file in the working directory), or `Cherenkov:a:b` |
| `-w`, `--pulse-width` [ns] | 0 | `pulse_width_ns` |
| `-d`, `--pulse-delay` [ns] | 0 | `pulse_delay_ns` |
| `--seed` | 1 | `seed` |

## Checks

- Directions inside the cone (within 1e-6), |position| = R, mean and variance of cos(theta)
  (4 sigma), entry points uniform over the projected disc (KS), mean time against `-E[zi] R / c` (4 sigma), a table spectrum's mean and variance.
- Invalid source parameters are rejected when the plan is built, also with `--photons 0`:
  `--radius -1` and `--viewcone 200` fail with an error.
- Reference-program run (`--atmosphere 1 --photons 200000`) wrote zero bunches. Cause not
  found yet, so no comparison.
