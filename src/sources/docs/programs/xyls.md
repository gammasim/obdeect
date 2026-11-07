# xyls

Parallel downward light beam scanned along two lines over the camera area
(`sim_telarray/LightEmission/xyls.cc`, a test program). Tool `dlamb-xyls`, programme
`cpp/include/dlamb/programs/xyls.hpp`.

## Options (original units, cm)

| Original | Default | dlamb |
| --- | --- | --- |
| `-r rcam` | 123 | `radius_cm` (also sets the scan length and the telescope radius) |
| `-x xpos` | 0.3 x 123 | `x_cm` (default from 123 cm, as in the original) |
| `-y ypos` | 0.2 x 123 | `y_cm` |
| `-h height` | 0.05 x 123 | `height_cm` |
| (fixed) | 100 photons, bunch 1, 2150 m, atmosphere 1 | `photons_per_position`, `bunch_size`, `altitude_m`, `atmosphere` |
| (new) `--seed` | 1 | `seed` |

## Mapping

- One `Emitter` per scan position, direction (0, 0, -1), `FixedAngle{0}`, 400 nm,
  delta pulse. Each position is emitted with its own seed.
- Scan loops as in `xyls.cc`: along x at y = y0 for `0.4 r` in steps of `0.001 r`,
  then along y at x = x0 for `0.2 r`.
- Count per position: `100 / bunch` bunches of weight `bunch` (exact for these defaults).

## Checks

- Default: total rows `100 x positions`; first scan point at (0.369, 0.246) m;
  direction (0, 0, -1); z = 0.0615 m.
- Original (container): 60100 bunches, same as dlamb.
