# beam-scan

Parallel downward light beam scanned along two lines over the camera area. Tool
`obdeect sources beam-scan`, programme `cpp/include/obdeect/sources/programs/beam-scan.hpp`. Replaces the earlier
xyls program; options keep their names.

## Options (cm)

| Option | Default | obdeect sources field |
| --- | --- | --- |
| `-r rcam` | 123 | `radius_cm` (also sets the scan length and the telescope radius) |
| `-x xpos` | 0.3 x 123 | `x_cm` |
| `-y ypos` | 0.2 x 123 | `y_cm` |
| `-h height` | 0.05 x 123 | `height_cm` |
| (fixed) | 100 photons, bunch 1, 2150 m, atmosphere 1 | `photons_per_position`, `bunch_size`, `altitude_m`, `atmosphere` |
| `--seed` | 1 | `seed` |

## Geometry

- One `Emitter` per scan position, direction (0, 0, -1), `FixedAngle{0}`, 400 nm,
  delta pulse. Each position is emitted with its own seed.
- Scan: along x at y = y0 for `0.4 r` in steps of `0.001 r`, then along y at x = x0 for
  `0.2 r`.
- Count per position: `100 / bunch` bunches of weight `bunch`.

## Checks

- Default: total rows `100 x positions`; first scan point at (0.369, 0.246) m;
  direction (0, 0, -1); z = 0.0615 m.
- Reference-program run: 60100 bunches, same as obdeect sources.
