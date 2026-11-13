# obdeect sources beam-scan (implementation)

User options and examples: [beam-scan](../programs/beam-scan.md). Writer B, producer
`obdeect sources beam-scan`.

## Validation

- `-r` must be finite and positive: `beam-scan: radius must be finite and positive`.
- The defaults of `-x`, `-y`, `-h` are derived from 123 cm and do not follow a later `-r`.
- Photons per position: finite, `>= 0`, `<= 2^53`, else `beam-scan: photons per position must
  be finite, >= 0 and at most 2^53`.

## Plan

- Fixed: run number 1, observation level 2150 m, `atmosphere_id` 1, 100 photons per position,
  bunch 1, wavelength 400 nm (`wavelength_min = wavelength_max = 400`), one telescope id 1 at
  (0, 0, 0) with radius `r/100` m, 1 event, one case.
- Scan positions in cm, by repeated addition in double precision (this fixes the count):
  `x = xpos; x < xpos + 0.4 r; x += 0.001 r` at `y = ypos`; then
  `y = ypos; y < ypos + 0.2 r; y += 0.001 r` at `x = xpos`.
- Position `k = 0, 1, ...` (over both lines) is one emission: emitter at `(x, y, height)/100` m,
  pointing `(0, 0, -1)`, FixedAngle(0), Fixed 400 nm, Gauss sigma 0 (`t = 0`),
  `Source.photons = 100`, bunch 1, seed `seed + k`, no target.
- Total rows `100 x positions` (default: more than 500 and fewer than 700 positions), all
  `w = 1`, `cz = -1`, `z = 0.0615`.
- `XylsOptions` has the fields `photons_per_position` and bunch size: with `bunch 2`,
  `photons 6`, each position gives 3 rows of `w = 2`.
