# led

An LED in front of selected camera pixels, to test trigger definitions. Replaces the earlier
pixled program; options keep their names. Command: `obdeect sources led --camera-file FILE (--use LIST |
--use-trg-rand N | --use-trg-all) [options] > out.gst`. `--list-only` prints the pixel
combinations instead.

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
| `--spectrum W`, `--lightpulse S` | `400`, `Simple:0` | as in flasher |
| `--run N`, `--altitude M`, `--atmosphere N` | `1`, `2150`, `1` | run keys |
| `--list-only` (`-L`) | off | print `(@)` labels, write no photons |
| `--quiet` (`-q`) | | accepted; the program prints nothing else |
| `--seed N` | `1` | seed of the random choices and photons |
| `-o FILE` | stdout | output file |

## Geometry

- Frame: `camera` (C). The LED is at the pixel centre `(x, y, 0)` in the pixel-file
  coordinates, pointing `(0, 0, -1)`, angle `Parallel`. The photons already lie on the
  camera plane, so no projection is needed.
- Count: `photons / bunchsize` per pixel; one emission per lit pixel and per slave.
- Events: all listed pixels (no `--require`), `--random` draws, or one event per
  combination of `--require` pixels in increasing order. Required (`+`) pixels are always
  lit.

## Input rules

- `--use-file` is not supported.
- Trigger lists: the text after `of` is used, `[..]` removed unless `--slaves`.
- Pixel lists accept ids, `+`, `[..]` and `%`.
- Random choices use Philox keyed on the seed.
- A non-isotropic pixel outside the fiducial sphere is rejected (no camera-plane disc
  target).

## Checks

- `ctest`: `test_led` (camera-plane positions and directions, combinations, slaves,
  trigger groups, random pixels, errors).
- Reference-program run (`scripts/compare_led.sh`, 3 events of 2 pixels, 10 photons each):
  20 photons per event in both; direction (0, 0, -1) and time 0 in both.
