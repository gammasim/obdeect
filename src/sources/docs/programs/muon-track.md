# muon-track

Cherenkov light of one muon-like track, as a sky-frame photon stream. Replaces the earlier
fake-muon program; options keep their names.

```bash
dlamb-muon-track -x 70 -h 3000 --theta 0 -o muon-track.ps
dlamb-muon-track --help
```

## Geometry and yield

- The track runs from `-h` above the observation level to the impact point
  `(x, y, 0)`, direction `(-sin theta cos phi, -sin theta sin phi, -cos theta)`.
- Photons per metre into 4 pi (Frank-Tamm):
  `2 pi alpha q^2 sin^2(theta_c) (1/wlmin - 1/wlmax)`, with `cos(theta_c) = 1/(n beta)`.
- Emission is at a fixed angle `theta_c` to the track, uniform azimuth, uniform along
  the segment. Wavelengths follow 1/lambda^2 between `wlmin` and `wlmax`, exactly (no
  tabulated grid).
- Emission time is `-D/(beta c)` for a point at distance `D` along the track from the
  observation level, so `t = 0` is the arrival of the particle there.
- Only photons that reach the telescope sphere (radius `-r`, centre at the origin) are
  generated, on straight rays. Weight = fraction of the azimuths that hit (exact), so the
  weight sum is the expected number of photons on the sphere. Photons are single
  (weight 1 or the arc fraction), not bunches.
- A step with `cos(theta) > 1` advances the track (no effect for beta = 1).
- The default impact point (0, 0) gives a track through the sphere centre. Its Cherenkov
  ring lies at about `h tan(theta_c)`, about 73 m for 3 km, so most photons come from the
  lower part of the track (below about 620 m, where the sphere subtends more than theta_c).

## Output

- Frame `sky`: positions relative to the observation level, z up, metres. No refraction
  and no extinction: dlufd applies them on the way to the telescope.
- One run, `--events` events, one telescope record (the sphere). `refraction_applied = 0`.

## Options

| Option | Default | Note |
| --- | --- | --- |
| `-r`, `--radius` | 15 m | sphere radius |
| `-x`, `-y` | 0 | impact point |
| `-h` | 3000 m | start altitude |
| `--theta`, `--phi` | 0 | degrees |
| `-n`, `--n1` | 1.000292 | index, or n-1 |
| `-A const` or `--index-profile FILE` | constant | see below |
| `-a` | 2150 m | written to the run record only |
| `-b` | 5 | accepted, no effect |
| `-l`, `-L` | 240, 1000 nm | 200 to 2000 nm |
| `-q`, `-m`, `-p`, `--beta` | -1, 0.105658 GeV, -, 1 | `-p` sets beta from the final mass (order-independent; not with `--beta`) |
| `--primary` | 6 (mu-) | sets mass and charge (`mass_charge`) |
| `-N`, `--events` | 1 | number of events, 1 to 4294967295 |
| `--run` | 1 | run number (run record) |
| `--enable-scatter` | off | rejected (needs the density profile) |
| `--rscat`, `-I` | 0 | rejected |
| `--seed` | 1 | explicit random seed |

## Refractive index

The default index is constant (`-n`). `--index-profile FILE` (lines `altitude_m
n_minus_1`, from simulation-models) gives an altitude-dependent index: the track is
stepped (10 m steps, index at each mid-point).

## Example

```bash
dlamb-muon-track -x 70 -N 10 --seed 3 -o muon-track.ps 2> summary.txt
dlamb-muon-track -x 70 --index-profile atm_index.txt | ps-cat --summary -
```
