# ls-beam

Laser beam crossing the field of view of a telescope (`sim_telarray/LightEmission/ls-beam.cc`).
Tool `dlamb-ls-beam`, programme `cpp/include/dlamb/programs/ls-beam.hpp`.

## What changes

The original emits the light **scattered** from the beam towards the sphere
(`LightSource::Emit`, Rayleigh profile, `brightness x exp(-h/7.5 km)` photons per step).
dlamb emits the **primary** photons along the beam inside the field-of-view cone.
Scattering is dlufd's single-scattering effect (plan step 6b).

| Original | dlamb |
| --- | --- |
| step point: `N = brightness x exp(-h/7500 m) x step` photons into 4 pi, Rayleigh | step point: `N = primary_photons_per_m x step` primary photons along the beam |
| `--brightness` (photons/m into 4 pi, ground) | replaced by `--primary-per-m` (primary photons per m) |
| `--angular-distribution` (Rayleigh) | recorded in the run key `scattering_angular` for dlufd |
| time = path / c_air | path / c (vacuum); group delay is dlufd's |

## Mapping of the yield to dlufd

Original scattered yield per metre at height h:
`B(h) = brightness x exp(-h/7.5 km)`.
dlufd computes, for primary weight `w` at height h:
`scattered per m = sum(w) x beta_s(h) x (phase function) x (transmission)`.
With uniform primaries, `sum(w) per m = primary_photons_per_m` and so
`primary_photons_per_m = brightness / beta_s(ground)` reproduces the ground yield
(brightness of the original, `beta_s` the scattering coefficient per m at the laser
wavelength). The original's `exp(-h/7.5 km)` is replaced by dlufd's `beta_s(h)`
from the atmosphere profile.

Parameters dlufd does not have: `beta_s(ground)` at 355 nm for the aerosol part (only
molecular scattering follows from the density profile), the aerosol single-scattering
albedo and phase-function asymmetry g. These come from simulation-models (step 6b).

## Options

| Original | Default | dlamb |
| --- | --- | --- |
| `--events` | 1 | `events` |
| `--step` [m] | 0.1 | `step_m` (the original's `steplength_cm` = 100 x step) |
| `--bunchsize` | 3 | `bunch_size` (below 1 used as 1) |
| `--spectrum` [nm] | 355 | `spectrum` |
| `--lightpulse` | Simple:0 | `pulse` |
| `--laser-position` [m] | 1000,0,0 | `laser_position_m` |
| `--laser-theta`, `--laser-phi` | 60, 180 | `laser_theta_deg`, `laser_phi_deg` (CLI phi astronomical, converted as the original) |
| `--telescope-theta`, `--telescope-phi` | 30, 0 | same, telescope axis |
| `--telescope-radius` [m] | 12.5 | `telescope_radius_m` |
| `--fov` [deg] | 15 | `fov_deg` (full cone opening angle) |
| `--altitude`, `--run`, `--atmosphere` | 2150, 1, 1 | `altitude_m`, `run_number`, `atmosphere` (recorded) |
| (new) `--primary-per-m` | 1e4 | `primary_photons_per_m` (sampling choice, not physics) |

Notes on the original: the help text says `--altitude` default 2100, the code uses
2150. `--bunches` is in the help but not parsed.

## Geometry

Beam and cone are solved as in `ls-beam.cc` (rotations with `Rotate_z(-phi)`, `Rotate_y(-theta)`,
quadratic in the cone). The emitted segment is the one between the two cone-surface roots,
as in the original; the original does not check the nappe, so for a beam crossing both
nappes the segment includes points below the telescope.

## Checks

- Default: all step points lie in the cone (angle to axis = fov/2 at the ends, within 1e-9).
- Non-default (fov 20, step 0.5 m, Gauss sigma 0.5 ns): residual time moments match the pulse.
- Beam that misses the cone is rejected.
