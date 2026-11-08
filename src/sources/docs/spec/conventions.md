# Conventions

## Units and frames

- Lengths m, times ns, wavelengths nm, angles deg at every user interface (tool options may
  use cm where stated). Internally angles are rad; speed of light `c = 0.299792458 m/ns`.
- Stream frames (gstrickd codes): sky 0, ground 1, telescope (T) 2, camera (C) 3. Sources may
  emit in sky, T or C only.
- Directions are unit propagation directions; `cz < 0` is downward. Sky frame: origin at the
  observation level above the telescope, z up. All z values in the sky frame are heights above
  the observation level; the column `zem` is `observation_level_m + z`.
- Angular profiles use `u = 1 - cos(theta)`, theta measured from the pointing axis, so
  `dOmega = 2 pi du`, `u in [0, 2]`.
- FWHM to sigma: `sigma = FWHM / 2.35482004503095`.

## Random numbers (contract)

All draws are `gstrickd::uniform4(DrawKey{seed, stage, event, telescope, photon, stream})`,
returning four doubles in (0, 1) (Philox-4x32-10). Stage `0x444C0001` (gstrickd
`stage_sources`) for every photon, count and track draw; stage `0x444C0002` (`stage_programs`)
for program-level choices. Streams of one photon (each one `uniform4` call, never shared):

| Stream | Id | Used for |
| --- | --- | --- |
| direction | 1 | emission angle (`u[0]`) and azimuth (`u[1]`) of emitters without target; track azimuth; background direction |
| target | 2 | point on the target disc (`u[0]` radius, `u[1]` azimuth) or star entrance disc; track position along the segment; background emission point |
| wavelength | 3 | wavelength |
| pulse | 4 | emission time |
| count | 5 | photon count, photon index 0 only, `u[0]` |

Photon results depend only on the key, never on chunking, order or threads. Counts, bunch
numbers, weights, record order and event numbering are exact; positions and times agree
statistically only if the samplers differ, but the key-to-stream assignment above is the
contract (counts use `uniform4` exactly as specified so that they are reproduced bit for bit).

Seeds: each emission of a tool gets its own seed so emissions never share draws. Two rules
exist and are part of the contract:

- `derive_seed(base, i)` (flasher, led) and `source_seed(base, i)` (muon-track) are the same
  function (splitmix64): `z = base + 0x9E3779B97F4A7C15 * (i + 1)` (mod 2^64);
  `z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9`; `z = (z ^ (z >> 27)) * 0x94D049BB133111EB`;
  result `z ^ (z >> 31)`.
- point-source, beam-scan, laser-beam use `seed + i` (mod 2^64); night-sky uses `seed`.

Program-level uniform: `program_uniform(seed, event, draw) = uniform4({seed, 0x444C0002, event,
0, draw, 1})[0]`.

## Rounding to counts

A mean `m >= 0` becomes `floor(m) + (u < m - floor(m) ? 1 : 0)` with `u = uniform4({seed,
0x444C0001, event, telescope, 0, 5})[0]`. Not a Poisson draw.

## float32 rule

Photon columns x y z cx cy cz t w wavelength (and zem) are float32. All computation is in
double; each value is rounded to float32 once when written. Any later arithmetic on a
written column (bunch weight, time offset, camera-plane projection, `zem`) is on the float32
values as stated on the pages below.

## Rejecting input

Invalid input raises `std::invalid_argument` (C++) or `ValueError` (Python); fields are never
silently clamped, dropped or ignored, except where a page says "used as 1" (bunch sizes below
1).

## Magnitude limits

`max_magnitude = 1e30`, `max_photon_count = 2^53`. A length, time or wavelength is "within
limits" if it is finite and `|x| <= 1e30`; a "positive-within-limits" value additionally is
`>= 1e-30`.
