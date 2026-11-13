# beam-scan

Parallel downward light beam scanned along two lines over the camera area.

```console
$ obdeect sources beam-scan -r 123 -o beam-scan.ps
```

- The beam moves along x at `y = y0` for `0.4 r`, then along y at `x = x0` for `0.2 r`, in
  steps of `0.001 r`. Replaces the earlier xyls program; options keep their names.
- 100 photons per position, one photon per bunch, 400 nm, observation level 2150 m.
- The defaults of `-x`, `-y`, `-h` are fractions of 123 cm and do not follow a later `-r`.

## Options (cm)

| Option | Default | Meaning |
| --- | --- | --- |
| `-r rcam` | 123 | camera radius; sets the scan length and the telescope radius |
| `-x xpos` | `0.3 x 123` | scan start x |
| `-y ypos` | `0.2 x 123` | scan start y |
| `-h height` | `0.05 x 123` | beam height |
| `--seed` | 1 | seed (position `k` uses `seed + k`) |
| `-o FILE` | stdout | output file |

## Checks

- Default: total rows `100 x positions`; first scan point at (0.369, 0.246) m; direction
  (0, 0, -1); z = 0.0615 m.
- Reference-program run: 60100 bunches, same as obdeect sources.
