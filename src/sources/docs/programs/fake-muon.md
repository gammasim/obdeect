# fake-muon

Cherenkov light of one muon-like track, as a sky-frame photon stream. Port of
`sim_telarray/LightEmission/fake-muon.cc`.

```bash
dlamb-fake-muon -x 70 -h 3000 --theta 0 -o fake-muon.ps
dlamb-fake-muon --help
```

## Geometry and yield

- The track runs from `-h` above the observation level to the impact point
  `(x, y, 0)`, direction `(-sin theta cos phi, -sin theta sin phi, -cos theta)`.
- Photons per metre into 4 pi (Frank-Tamm, `fake-muon.cc:504`, `:594`):
  `2 pi alpha q^2 sin^2(theta_c) (1/wlmin - 1/wlmax)`, with `cos(theta_c) = 1/(n beta)`.
- Emission is at a fixed angle `theta_c` to the track, uniform azimuth, uniform along
  the segment. Wavelengths follow 1/lambda^2 between `wlmin` and `wlmax`.
- Emission time is `-D/(beta c)` for a point at distance `D` along the track from the
  observation level, so `t = 0` is the arrival of the particle there.
- Only photons that reach the telescope sphere (radius `-r`, centre at the origin) are
  generated. Weight = fraction of the azimuths that hit (exact), so the weight sum is the
  expected number of photons on the sphere.

## Output

- Frame `sky`: positions relative to the observation level, z up, metres. No refraction
  and no extinction: dlufd applies them on the way to the telescope.
- One run, `--events` events, one telescope record (the sphere). `refraction_applied = 0`.

## Options

| Option (original) | dlamb | Default | Note |
| --- | --- | --- | --- |
| `-r`, `--radius` (`:240`) | `-r`, `--radius` | 15 m | sphere radius |
| `-x`, `-y` (`:241`) | `-x`, `-y` | 0 | impact point |
| `-h` (`:245`) | `-h` | 3000 m | start altitude |
| `--theta`, `--phi` (`:249`) | same | 0 | degrees |
| `-n`, `--n1` (`:250`) | `-n` | 1.000292 | index, or n-1 |
| `-A` (`:258`) | `-A const` or `--index-profile FILE` | constant | see below |
| `-a` (`:239`) | `-a` | 2150 m | written to the run record only |
| `-b` (`:247`) | `-b` | 5 | accepted, no effect |
| `-l`, `-L` (`:248`) | `-l`, `-L` | 240, 1000 nm | 200 to 2000 nm |
| `-q`, `-m`, `-p`, `--beta` | same | -1, 0.105658 GeV, -, 1 | `-p` sets beta from the final mass (order-independent; not with `--beta`) |
| `--primary` (`:327`) | `--primary` | 6 (mu-) | sets mass and charge (`mass_charge`) |
| `-N`, `--events` | same | 1 | number of events, 1 to 4294967295 |
| `--run` | same | 1 | run number (run record) |
| `--enable-scatter` (`:275`) | rejected | off | not ported (needs the density profile) |
| `--rscat`, `-I` | rejected | 0 | disabled in the original / not applicable |
| `--seed` | new | 1 | explicit random seed |

## Differences from the original

- Refractive index: the default is constant (`-n`). The original's default atmosphere
  (`-A 1`, tropical) gives an altitude-dependent index. Pass `--index-profile FILE`
  (lines `altitude_m n_minus_1`, from simulation-models) to get the stepped track
  (10 m steps, index at each mid-point, `fake-muon.cc:529`).
- The original skips a step with `cos(theta) > 1` without advancing the track. The port
  advances it (no effect for beta = 1).
- The default impact point (0, 0) gives a track through the sphere centre. Its Cherenkov
  ring lies at about `h tan(theta_c)`, about 73 m for 3 km, so most photons come from the
  lower part of the track (below about 620 m, where the sphere subtends more than theta_c).
- The original spectrum is tabulated on a 5 nm grid; the port uses the exact 1/lambda^2.
- Refraction in the sphere test: the original tests the refracted impact points, the port
  tests straight rays. The size of the effect was not quantified; it is expected to be
  small for near-vertical tracks.
- Bunches: the original emits bunches of `-b` photons. The port emits single photons with
  weight 1 (or the arc fraction), with the same mean count.

## Example

```bash
dlamb-fake-muon -x 70 -N 10 --seed 3 -o fake-muon.ps 2> summary.txt
dlamb-fake-muon -x 70 --index-profile atm_index.txt | ps-cat --summary -
```
