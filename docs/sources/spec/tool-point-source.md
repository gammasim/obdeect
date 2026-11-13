# obdeect sources point-source

Point light source above the observation level, sky frame. Replaces an earlier point-source
program. Writer B, producer `obdeect sources point-source`.

## Options (units cm unless stated)

| Option | Default | Rules |
| --- | --- | --- |
| `-R run` | 1 | integer; the written run number is `run` if `> 0`, else 1 |
| `-N n` | 1 | events per case, integer `0..2^32-1` |
| `-r r` | 1500 | fiducial radius [cm] of the default telescope |
| `-x`, `-y` | 0 | source position [cm] |
| `-z z` | `1000e2` (= 1e5 cm = 1000 m) | height [cm] |
| `-n N[,..]` | `1e7` | photons into 4 pi per event; one case per value |
| `-b b` | 1 | photons per bunch; below 1 used as 1 |
| `-s spec` | `400` | spectrum: number 200..2000, `Cherenkov:a:b`, or table file |
| `-p pulse` | `Gauss:0.0` | bare number = Gauss sigma; see cli-conventions |
| `-a ang` | `isotropic` | angular text, tables allowed |
| `-h alt` | 2150 | observation level [m] |
| `-A n` | 1 | atmosphere id, recorded as `atmosphere_id` |
| `-d cx,cy,cz` | `0,0,-1` | pointing; normalised; all zero is an error (`point-source: direction vector must not be zero`) |
| `--telpos-file f`, `--corsika-telpos-file f` | one telescope | later one wins |
| `--seed n` | 1 | |
| `-o fname` | stdout | |
| `--help` | | usage, exit 0 |

Spectrum, pulse and angular texts are parsed after all options. Options with a trailing
missing value: `<opt>: missing value`.

## Run record

Telescopes: with no telescope file, one telescope id 1 at (0,0,0) with radius `r/100` m. With a
file: ids 1.. in file order, centre and radius divided by 100. `wavelength_min/max` from the
spectrum range (written only if max > 0). `atmosphere_id` = `-A`. `observation_level_m` =
`-h`. No extra keys.

## Cases and emissions

One case per `-n` value; each case is written `-N` times; event numbers consecutive from 1.
Case `c` has one `SourceEmission`: Emitter at `(x, y, z)/100` m with `pointing = normalised(d)` (d as given,
default (0,0,-1)), angular profile
`-a`, spectrum, pulse, frame sky, `Source.photons = photons[c] / b`, `b = max(1,-b)`,
seed `seed + c`, time offset 0. Photons per case must be finite, `>= 0` and `<= 2^53`
(`point-source: photons must be finite, >= 0 and at most 2^53`).

Target: `rmax = fiducial_super_radius(telescopes)`. If `|position| > 1.02 rmax`: Isotropic ->
`SphereTarget(origin, 1.02 rmax)` (so the photon count per event is `photons/b *
(1 - sqrt(1 - (1.02 rmax/d)^2))/2`); FixedAngle -> no target; any other profile -> error
`point-source: a non-isotropic profile beyond 1.02 x the fiducial radius is not supported
(use isotropic or Cone/Parallel)`. Otherwise no target.

Written for every telescope block identical photons (telescope key 0), `w = b` (float32),
`zem = observation_level + z`.

## Acceptance (`test_point-source`)

Default: one SphereTarget; mean count over 4000 events matches `1e7 * 0.5 * (1 -
sqrt(1 - (1.02*15/1000)^2))` to `4*0.5/sqrt(4000)`; every row `w = 1`, `wavelength = 400`,
`t = 0`, `z = 1000`, `cz <= -cos_max`. Non-default: Cone 30 deg, `-b 2`, `2e5` photons, 2 events
=> exactly `2*100000` rows, `w = 2`, `cz = -cos 30deg`, exponential pulse moments. Lambertian at
default distance is rejected; negative photons rejected; `zem - z = altitude`; telescope file
cases (see cli-conventions).
