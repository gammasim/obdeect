# figure-eight

Isotropic point light moved along a figure-eight 1000 m above the observation level,
emitting at each position. Port of sim_telarray `LightEmission/octo.cc`. The original is a
demonstration with no options. The CLI exposes only `--photons`, `--seed` and `-o`; the
other parameters are library-only.

```bash
dlamb-figure-eight --seed 1 -o figure-eight.ps
```

## Geometry and yield

- 41 positions on the circle centred at `(20 m, 0)`, t = 0 to 20 ns, then 40 positions on the
  circle centred at `(-20 m, 0)`, t = 20.5 to 40 ns (`octo.cc:107-124`). Altitude 1000 m.
- At each position: isotropic emission, 400 nm, instantaneous, `8e7` photons into 4 pi.
- Only photons that reach the sphere of radius 12.5 m at the origin are generated. The
  count is `N * (1 - cos(theta_max)) / 2`, `sin(theta_max) = R / distance` (exact). This is
  the same expectation as figure-eight's cone pre-selection followed by the sphere test.
- Each position is a source with its own seed (`source_seed(seed, k)`).

## Output

- Frame `sky`, one event, one telescope record (sphere, `r_m = 12.5`), weight 1 per photon.
- Output photons are single photons. The original's bunches of 2 photons are not kept.

## Options

| Parameter (original) | dlamb | Default |
| --- | --- | --- |
| `alt_m` (`:76`) | observation level, run record only | 2150 m |
| `r_tel_m` (`:77`) | sphere radius | 12.5 m |
| `dist` (`:78`) | altitude of the light | 1000 m |
| `nphot_per_pos` (`:79`) | `--photons` | 8e7 |
| `rc` (`:86`) | ring radius (library option) | 20 m |
| SimpleLightSource default (`IactLightEmission.hh:591`) | wavelength (library option) | 400 nm |
| - | `--seed` | 1 |

sim_telarray's ls-test is `dlamb-figure-eight --photons 4e7`.

## Check against the original

Container `ghcr.io/gammasim/simtools-dev`, `octo` built from `sim_telarray/LightEmission`
(with `FILEOPEN_DEFAULT_PATH` set):

| | Original | dlamb |
| --- | --- | --- |
| Expected (original printout: bunches; dlamb: photons) | 126461 bunches (= 252922 photons) | 252935 photons |
| Stored | 126320 bunches (= 252640 photons) | 252941 photons |

The stored counts agree within 1 sigma of the original's Poisson fluctuation.
