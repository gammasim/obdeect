# led

LED in front of selected camera pixels, to test trigger definitions.

```console
$ obdeect sources led --camera-file cam.dat --use "0 1 2" --photons 500 > led.pst
$ obdeect sources led --camera-file cam.dat --use-trg-all --list-only
```

- Frame `camera` (C). The LED sits at the pixel centre `(x, y, 0)` in the pixel-file
  coordinates, pointing `(0, 0, -1)`. The photons already lie on the camera plane.
- Give exactly one of `--use`, `--use-trg-rand N` or `--use-trg-all`.
- `--list-only` prints the pixel combinations (`(@)` labels) and writes no photons.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--camera-file FILE` (`-c`) | required | camera definition (`Pixel`, `Trigger` lines) |
| `--use LIST` (`-u`) | | pixels: ids, `+` required, `[a,b]` slaves, `%` comment |
| `--use-trg-rand N` | | N random trigger groups |
| `--use-trg-all` | | all trigger groups |
| `--require N` (`--required`) | `0` (all) | pixels per event |
| `--random N` | | N distinct random pixels per event (uses `--require`) |
| `--slaves` | off | also fire pre-sum slaves |
| `--events N` | `1` | events per combination |
| `--photons N` | `500` | photons per pixel |
| `--bunchsize F` (`-b`) | `1` | photons per bunch (at least 1) |
| `--camera-radius R` (`-r`) | `123` | fiducial radius [cm] |
| `--spectrum W`, `--lightpulse S` | `400`, `Simple:0` | as in [flasher](flasher.md) |
| `--run N`, `--altitude M`, `--atmosphere N` | `1`, `2150`, `1` | run keys |
| `--list-only` (`-L`) | off | print `(@)` labels, write no photons |
| `--quiet` (`-q`) | | accepted; the program prints nothing else |
| `--seed N` | `1` | seed of the random choices and photons |
| `-o FILE` | stdout | output file |

## Input rules

- Pixel lists accept ids, `+` (required), `[..]` (slaves) and `%` (comment).
- Trigger lists: the text after `of` is used; `[..]` is removed unless `--slaves`.
- `--use-file` is not supported.
- A non-isotropic pixel outside the fiducial sphere is rejected (no camera-plane disc target).
- Events: all listed pixels (no `--require`), `--random` draws, or one event per combination of
  `--require` pixels in increasing order. Required (`+`) pixels are always lit.
