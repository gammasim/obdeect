# laser-beam

Laser beam crossing the field of view of a telescope. obdeect sources emits the **primary**
photons along the beam inside the field-of-view cone. The scattered light is computed from
them by [obdeect propagation](../../propagation/index.md) (single scattering, plan step 6b).

```console
$ obdeect sources laser-beam --events 1 --primary-per-m 1e4 -o laser.ps
```

Replaces the earlier ls-beam program; options keep their names.

## Options

| Option (aliases) | Default | obdeect sources field |
| --- | --- | --- |
| `--events` | 1 | `events` |
| `--step` [m] | 0.1 | `step_m` |
| `--bunchsize` | 3 | `bunch_size` (below 1 used as 1) |
| `--spectrum` [nm] | 355 | `spectrum` |
| `--lightpulse` | Simple:0 | `pulse` |
| `--laser-position` [m] | 1000,0,0 | `laser_position_m` (`x,y` or `x,y,z`) |
| `--laser-theta`, `--laser-phi` | 60, 180 | `laser_theta_deg`, `laser_phi_deg` (CLI phi astronomical, converted) |
| `--telescope-theta`, `--telescope-phi` | 30, 0 | telescope axis |
| `--telescope-radius` [m] | 12.5 | `telescope_radius_m` |
| `--fov` [deg] | 15 | `fov_deg` (full cone opening angle, in (0, 90)) |
| `--altitude`, `--run`, `--atmosphere` | 2150, 1, 1 | `altitude_m`, `run_number`, `atmosphere` (recorded) |
| `--primary-per-m` | 1e4 | `primary_photons_per_m` (sampling choice, not physics) |
| `--angular-distribution` | Rayleigh | recorded in the run key `scattering_angular` for obdeect propagation |

## Yield per metre

obdeect propagation computes, for primary weight `w` at height h:
`scattered per m = sum(w) x beta_s(h) x (phase function) x (transmission)`.
With uniform primaries, `sum(w) per m = primary_photons_per_m`. For a scattered yield `B` per
metre at the ground, set `primary_photons_per_m = B / beta_s(ground)`, with `beta_s` the
scattering coefficient per m at the laser wavelength. Propagation takes `beta_s(h)` from the
atmosphere profile.

Not available in obdeect propagation: `beta_s(ground)` at 355 nm for the aerosol part (only
molecular scattering follows from the density profile), the aerosol single-scattering albedo
and the phase-function asymmetry g. These come from simulation-models (step 6b).

## Checks

- Default: all step points lie in the cone (angle to axis = fov/2 at the ends, within 1e-9).
- Non-default (fov 20, step 0.5 m, Gauss sigma 0.5 ns): residual time moments match the pulse.
- A beam that misses the cone is rejected.
