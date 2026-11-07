# led

Port of sim_telarray `LightEmission/pixled.cc`: an LED in front of selected camera pixels,
to test trigger definitions. Command: `dlamb-led --camera-file FILE (--use LIST |
--use-trg-rand N | --use-trg-all) [options] > out.gst`. `--list-only` prints the pixel
combinations instead.

## Options

| Option | Default | Original | Meaning |
| --- | --- | --- | --- |
| `--camera-file FILE` (`-c`) | required | `pixled.cc:513, 652-667` | sim_telarray camera definition (`Pixel`, `Trigger` lines) |
| `--use LIST` (`-u`) | | `pixled.cc:583-590` | pixels: ids, `+` required, `[a,b]` slaves, `%` comment |
| `--use-trg-rand N` | | `pixled.cc:599-607` | N random trigger groups |
| `--use-trg-all` | | `pixled.cc:608-617` | all trigger groups |
| `--require N` (`--required`) | `0` (all) | `pixled.cc:618-619` | pixels per event |
| `--random N` | | `pixled.cc:620-624` | N distinct random pixels per event (uses `--require`) |
| `--slaves` | off | `pixled.cc:540-545` | also fire pre-sum slaves |
| `--events N` | `1` | `pixled.cc:629-630` | events per combination |
| `--photons N` | `500` | `pixled.cc:509, 571-572` | photons per pixel (the original `--help` says 2500000) |
| `--bunchsize F` (`-b`) | `1` | `pixled.cc:573-574` | photons per bunch (at least 1) |
| `--camera-radius R` (`-r`) | `123` | `pixled.cc:575-576` | fiducial radius [cm] |
| `--spectrum W`, `--lightpulse S` | `400`, `Simple:0` | `pixled.cc:579-582` | as in flasher |
| `--run N`, `--altitude M`, `--atmosphere N` | `1`, `2150`, `1` | `pixled.cc:625-628` | run keys |
| `--list-only` (`-L`) | off | `pixled.cc:546-551` | print `(@)` labels, write no photons |
| `--quiet` (`-q`) | | `pixled.cc:552-557` | accepted; the program prints nothing else |
| `--seed N` | `1` | new | seed of the random choices and photons |
| `-o FILE` | stdout | `pixled.cc:567-568` | output file |

## Mapping to dlamb

- Frame: `camera` (C). The LED is at the pixel centre `(x, y, 0)` in the pixel-file
  coordinates, pointing `(0, 0, -1)`, angle `Parallel`. The photons already lie on the
  camera plane, so no projection is needed.
- Count: `photons / bunchsize` per pixel; one emission per lit pixel and per slave.
- Events: as in the original. All listed pixels (no `--require`), `--random` draws, or one
  event per combination of `--require` pixels in increasing order (`pixled.cc:424-497`).
  Required (`+`) pixels are always lit.
- Bypass: meant for sim_telarray `Bypass_Optics=2`. sim_telarray (and `ps2simtel`) apply the
  y flip and camera rotation, as in the original.

## Differences to the original

- `--use-file` is not supported (the original's reader is an empty stub and always fails).
- Trigger lists: the text after `of` is used, `[..]` removed unless `--slaves`. The
  original drops the character after a removed `[..]` group; this port does not.
- Pixel lists are parsed with a simplified reader (same ids, `+`, `[..]`, `%`).
- Random choices use Philox keyed on the seed, not the global RNG.
- A non-isotropic pixel outside the fiducial sphere is rejected (no camera-plane disc
  target).

## Checks

- `ctest`: `test_led` (camera-plane positions and directions, combinations, slaves,
  trigger groups, random pixels, errors).
- Cross-check (`scripts/compare_led.sh`, 3 events of 2 pixels, 10 photons each): 20
  photons per event in both; direction (0, 0, -1) and time 0 in both.
