# point-source

Point light source above the telescope, sky-frame photons. Tool `obdeect sources point-source`,
programme `cpp/include/obdeect/sources/programs/point-source.hpp`. Replaces the earlier xyzls program;
options keep their names.

## Options

| Option | Unit | Default | obdeect sources field |
| --- | --- | --- | --- |
| `-R run` | | 1 | `run_number` |
| `-N nevents` | | 1 | `events` (events per case) |
| `-r rtelcam` | cm | 1500 | `telcam_radius_cm` (default telescope) |
| `-x`, `-y` | cm | 0 | `x_cm`, `y_cm` |
| `-z zpos` | cm | 1000e2 | `z_cm`, above the observation level |
| `-n nphot` | | 1e7 | `photons` (list: one case each) |
| `-b bsize` | | 1 | `bunch_size` (below 1 is used as 1) |
| `-s spectrum` | nm | 400 | `spectrum`: fixed (200..2000), `Cherenkov:a:b`, table |
| `-p pulse` | ns | Gauss:0.0 | `pulse`; bare number = Gauss sigma |
| `-a angdist` | deg | isotropic | `angular`: isotropic, Gauss, Rayleigh, Cone, Parallel, FilledCone, table |
| `-h alt_m` | m | 2150 | `altitude_m` |
| `-A atm` | | 1 | `atmosphere`, recorded only |
| `-d cx,cy,cz` | | 0,0,-1 | `direction` (normalised) |
| `--telpos-file`, `--corsika-telpos-file` | cm | none | telescopes `x y z r` |
| `--seed` | | 1 | `seed` |

## Geometry

- Source: `Emitter` at (x, y, z) with the given profile; `photons = N / bunch`,
  each row has `w = bunch`.
- Fiducial sphere: beyond `1.02 x R_max` an isotropic source emits only into the cone
  towards the origin that contains the sphere (`SphereTarget`).
- Cone, Parallel: no restriction.
- Other profiles beyond the sphere are rejected.
- Time: emission time (pulse). Refraction and extinction are not applied; obdeect propagation does them.

## Output

Stream: stream, run (run keys incl. telescope arrays, `refraction_applied = 0`),
then per event one telescope block per telescope, with the same photons for every
telescope (one emission, shared key). Frame sky, `z` above the observation level.

```console
$ obdeect sources point-source -N 50 -o point-source.ps
$ obdeect stream cat --summary point-source.ps
```

## Checks

- Mean count per event equals `1e7 x 0.5 (1 - cos(asin(1.02 R / d)))` (4 sigma).
- Non-default (Cone:30, bunch 2, Exp:2, 300 nm): exact counts, cone angle, weight,
  exponential moments.
- Reference-program run (50 events): 561.8 bunches per event vs 585.0 expected from the
  formula (not explained yet).
