# Conventions

## Units and frames

User-facing units are in [sources](../sources.md#frames-and-units). Here:

- Stream frame codes: sky 0, ground 1, telescope (T) 2, camera (C) 3. Sources emit in sky, T, C
  or G (`ground` runs: positions relative to the one telescope sphere centre, sphere rows).
- Sky frame: origin at the observation level above the telescope, z up. Sky z values are
  heights above the observation level; the column `zem` is `observation_level_m + z`.
- Angular profiles use `u = 1 - cos(theta)`, so `dOmega = 2 pi du`, `u in [0, 2]`.
- Internally angles are rad; `c = 0.299792458 m/ns`. FWHM to sigma: `sigma = FWHM / 2.35482004503095`.

## Random numbers (contract)

The generator, the key fields and the stage ids are described in
[](../../stream/random-numbers.md). All draws are `obdeect::stream::uniform4(DrawKey{seed, stage,
event, telescope, photon, stream})`, returning four doubles in (0, 1) (Philox-4x32-10). Stage
`0x444C0001` (obdeect stream `stage_sources`) for every photon and count draw. Streams of one
photon (each one `uniform4` call, never shared):

| Stream | Id | Used for |
| --- | --- | --- |
| direction | 1 | emission angle (`u[0]`) and azimuth (`u[1]`) of points and lines; sky-background direction |
| target | 2 | point on a disc (`u[0]` radius, `u[1]` azimuth); line position (`u[0]`); plane-wave disc point; sky-background start point |
| wavelength | 3 | wavelength |
| pulse | 4 | emission time |
| count | 5 | photon count, photon index 0 only, `u[0]` |

Photon results depend only on the key, never on chunking, order or threads. Counts, bunch
numbers, weights, record order and event numbering are exact. The key-to-stream assignment above
is the contract.

Seeds: each emission of a run has its own seed, so emissions never share draws:

- `source_seed(base, i)` (splitmix64): `z = base + 0x9E3779B97F4A7C15 * (i + 1)` (mod 2^64);
  `z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9`; `z = (z ^ (z >> 27)) * 0x94D049BB133111EB`;
  result `z ^ (z >> 31)`. `i` counts the emissions of the run in event order ([run](run.md#seeds-and-draws)).

## Rounding to counts

A mean `m >= 0` becomes `floor(m) + (u < m - floor(m) ? 1 : 0)` with `u = uniform4({seed,
0x444C0001, event, telescope, 0, 5})[0]`. Not a Poisson draw.

## float32 rule

Photon columns x y z cx cy cz t w wavelength (and zem) are float32. All computation is in
double; each value is rounded to float32 once when written. Any later arithmetic on a written
column (bunch weight, time offset, plane move, `zem`) is on the float32 values, as stated on the
[run](run.md) page.

## Rejecting input

Invalid input raises `std::invalid_argument` (C++) or `ValueError` (Python). Fields are never
silently clamped, dropped or ignored. The only dropped rows are line photons that miss a sphere
target ([sources](../sources.md#line-photons-with-a-sphere-target)).

## Magnitude limits

`max_magnitude = 1e30`, `max_photon_count = 2^53`. A length, time or wavelength is "within
limits" if it is finite and `|x| <= 1e30`; a "positive-within-limits" value additionally is
`>= 1e-30`.
