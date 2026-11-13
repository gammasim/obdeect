# obdeect sources point-source (implementation)

User options and examples: [point-source](../programs/point-source.md). Writer B, producer
`obdeect sources point-source`.

## Validation and defaults

- `-R`: the written run number is `run` if `> 0`, else 1.
- `-b` below 1 is used as 1.
- Direction `-d`: normalised; all zero is an error (`point-source: direction vector must not be zero`).
- Spectrum, pulse and angular texts are parsed after all options. A trailing missing value:
  `<opt>: missing value`.
- Photons per case must be finite, `>= 0` and `<= 2^53`
  (`point-source: photons must be finite, >= 0 and at most 2^53`).
- A bare number in `-p` is a Gauss sigma (mean 0).

## Run record

- Telescopes: without a telescope file, one telescope id 1 at (0, 0, 0) with radius `r/100` m.
  With a file: ids 1.. in file order, centre and radius divided by 100.
- `wavelength_min/max` from the spectrum range (written only if max > 0).
- `atmosphere_id` = `-A`; `observation_level_m` = `-h`. No extra keys.

## Cases and emissions

- One case per `-n` value; each case is written `-N` times; event numbers run from 1.
- Case `c`: one emitter at `(x, y, z)/100` m, `pointing = normalised(d)`, profile `-a`,
  `Source.photons = photons[c] / b`, seed `seed + c`, time offset 0, frame sky.
- Target: `rmax = fiducial_super_radius(telescopes)`. If `|position| > 1.02 rmax`:
  - isotropic: `SphereTarget(origin, 1.02 rmax)`, count per event `photons/b (1 - sqrt(1 - (1.02 rmax/d)^2)) / 2`;
  - FixedAngle: no target;
  - other profiles: error `point-source: a non-isotropic profile beyond 1.02 x the fiducial
    radius is not supported (use isotropic or Cone/Parallel)`.
- Every telescope block gets the same photons (telescope key 0), `w = b` (float32),
  `zem = observation_level + z`.

## Acceptance (`test_point-source`)

- Default: one SphereTarget; mean count over 4000 events matches
  `1e7 x 0.5 (1 - sqrt(1 - (1.02 x 15/1000)^2))` within `4 x 0.5 / sqrt(4000)`; every row
  `w = 1`, `wavelength = 400`, `t = 0`, `z = 1000`, `cz <= -cos_max`.
- Non-default (Cone 30 deg, `-b 2`, `2e5` photons, 2 events): exactly `2 x 100000` rows,
  `w = 2`, `cz = -cos 30 deg`, exponential pulse moments.
- Lambertian at default distance is rejected; negative photons rejected; `zem - z = altitude`;
  telescope-file cases: see [CLI conventions](cli-conventions.md).
