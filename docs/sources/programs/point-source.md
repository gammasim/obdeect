# point-source

Point light source above the telescope, sky-frame photons.

```console
$ obdeect sources point-source -N 50 -o point-source.ps
$ obdeect stream cat --summary point-source.ps
```

- Source at `(x, y, z)` (cm, `z` above the observation level), emitting `-n` photons into 4 pi
  per event. `-n` accepts a list: one case per value.
- Photons are written for every telescope; refraction and extinction are left to
  [obdeect propagation](../../propagation/index.md).

## Options (units cm unless stated)

| Option | Unit | Default | Meaning |
| --- | --- | --- | --- |
| `-R run` | | 1 | run number |
| `-N nevents` | | 1 | events per case |
| `-r rtelcam` | cm | 1500 | fiducial radius of the default telescope |
| `-x`, `-y` | cm | 0 | source position |
| `-z zpos` | cm | `1000e2` | height above the observation level |
| `-n nphot` | | `1e7` | photons into 4 pi per event (list: one case each) |
| `-b bsize` | | 1 | photons per bunch (below 1 used as 1) |
| `-s spectrum` | nm | 400 | fixed (200..2000), `Cherenkov:a:b`, or table file |
| `-p pulse` | ns | `Gauss:0.0` | pulse; bare number = Gauss sigma |
| `-a angdist` | deg | `isotropic` | isotropic, Gauss, Rayleigh, Cone, Parallel, FilledCone, table |
| `-h alt_m` | m | 2150 | observation level |
| `-A atm` | | 1 | atmosphere id (recorded only) |
| `-d cx,cy,cz` | | `0,0,-1` | direction (normalised) |
| `--telpos-file`, `--corsika-telpos-file` | cm | one telescope | telescopes `x y z r` |
| `--seed` | | 1 | seed |

## Checks

- Mean count per event equals `1e7 x 0.5 (1 - cos(asin(1.02 R / d)))` (4 sigma).
- Non-default (Cone:30, bunch 2, Exp:2, 300 nm): exact counts, cone angle, weight,
  exponential moments.
- Reference-program run (50 events): 561.8 bunches per event vs 585.0 expected from the
  formula (not explained yet).
