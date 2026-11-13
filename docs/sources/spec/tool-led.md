# obdeect sources led

LED in front of selected camera pixels (frame C). `obdeect sources led --camera-file FILE (--use LIST |
--use-trg-rand N | --use-trg-all) [options]`. Writer A, producer `obdeect sources led`.

## Options

| Option (aliases) | Default | Rules |
| --- | --- | --- |
| `--camera-file` (`-c`) | required | missing: `no camera definition file defined (--camera-file)`; unopenable: `cannot open camera configuration file <path>` |
| `--use` (`-u`) | | pixel list |
| `--use-trg-rand N` | 0 | N random trigger groups |
| `--use-trg-all` | off | all trigger groups |
| `--require` (`--required`) N | 0 | pixels per event, 0 = all, `>= 0` |
| `--random N` | off | sets random mode and `require = N` |
| `--slaves` | off | also fire pre-sum slaves |
| `--events N` | 1 | `>= 0` |
| `--photons N` | 500 | per pixel, `>= 0` |
| `--bunchsize F` (`-b`) | 1 | `> 0`; below 1 used as 1 |
| `--camera-radius R` (`-r`) | 123 | cm, `> 0` |
| `--spectrum`, `--lightpulse` (`--pulse`) | `400`, `Simple:0` | flasher syntax |
| `--run` (`-R`), `--altitude` (`-a`), `--atmosphere` (`-A`), `--seed` | 1, 2150, 1, 1 | as flasher |
| `--list-only` (`-L`) | off | print labels, write no photons |
| `--quiet` (`-q`) | | accepted, no effect |
| `-o FILE` | stdout | |

The angular profile is fixed to `Parallel`. Validation: exactly one of `--use`,
`--use-trg-rand`(>0), `--use-trg-all` (`give exactly one of ...`); at least one pixel in the
camera file; `events >= 0`, `require >= 0`; `photons >= 0`, `bunch > 0`, `radius > 0` (all
messages prefixed `led: `).

## Camera file

Line-oriented, leading blanks ignored, `#` lines and blank lines skipped. Keyword matching is
case-insensitive on the first word for `Pixel`, `Rotate`, `Pixtype` (the latter two skipped);
`Pixel <id> <type> <x_cm> <y_cm> ...` defines a pixel (a line that does not parse as `int int
double double` prints `led: problem with pixel definition: <line>` to stderr and is skipped);
negative id or duplicate id are errors. Any other line containing the text `Trigger`
(case-sensitive) is a trigger group.

Trigger text cleaning (`clean_trigger_list(line, keep_slaves)`): take the first line; if it
contains `of`, start after the first `of` and skip following blanks/tabs; copy the rest, and
unless `keep_slaves` drop every `[...]` group. The pixel file is read after all options, so
`--slaves` decides.

## Pixel list syntax

Items separated by blanks/tabs: `[+]id` optionally followed by `[a,b,..]` (slaves, commas and
blanks separate ids). `%` or newline ends the list. `+` marks a required pixel. Errors:
`led: invalid in list of pixels to use: '<rest>'`, `led: pre-sum pixel group not terminated
in '<text>'`. A repeated id is skipped with a warning on stderr (`duplicated pixel <id> in use
list, used only once`). `--use` text is cleaned like a trigger line first.

## Lists and events

Lists: `--use-trg-all`: one list per trigger group (error if none: `the camera definition has
no trigger groups`); `--use-trg-rand N`: N lists, each a group drawn by `j = floor(G *
program_uniform(seed, 0, draw++))` repeating until `j < G` (draw counter shared and increasing
across the N draws); `--use`: one list. An empty list: `led: failed to obtain a pixel list`.

Per list, `nreq = require` or the list size if `require == 0`. If `nreq` exceeds the list
size, print `led: required <n> pixels but the list has <m>; list skipped` to stderr and skip
the list. Otherwise repeat `events` times:

- `nreq == list size`: one event with all pixels;
- `--random`: one event with `nreq` distinct pixels: draw `p = floor(L * program_uniform(seed,
  event_number, draw++))` (draw counter restarts at 0 per event), ignore `p >= L` and
  already-used ids; fire in draw order;
- else: one event per combination of `nreq` pixels in lexicographic order of list position,
  starting from the `+` pixels (in list order; if there are at least `nreq` of them, that set
  alone is the single combination) and adding non-required pixels in increasing list order until
  `nreq` pixels are reached; note the whole combination set is emitted in
  each of the `events` repetitions.

Event numbers are consecutive from 1 across all lists. A pixel not in the camera file: `led:
pixel ID <n> is not defined in the camera file`.

## Emissions

Each lit pixel (and, with `--slaves`, each of its slaves, without labels) is one emission, in
pixel order, emitter at `(x_cm/100, y_cm/100, 0)` m, pointing `(0,0,-1)`, fiducial rule with
`R = camera_radius/100` (FixedAngle: never a target), `Source.photons = photons / b`, bunch
`b = max(1,bunchsize)`, seed `derive_seed(seed, j)` over a running emission counter,
telescope key 1, frame C, run radius `R`, photons written at `z = 0`, `cx = cy = 0`, `cz = -1`,
`t` from the pulse. Pixel positions use pixel-file coordinates (no y flip).

## `--list-only`

For every event print to stdout a line `(@) ` followed by, for each label in ascending
numeric order, a blank and the id (i.e. `(@)  0 1` with two blanks before each id), and a
newline. Labels are the listed pixels only (slaves never). No stream is written.

## Acceptance (test_led)

Combinations of "0 1 2" with require 2: three events of 2x10 photons, sets {0,1},{0,2},{1,2};
trigger groups with and without slaves (40 vs 30+20 rows; labels do not include slaves);
random mode: 5 events of two distinct pixels; `+0 1` fires both; require above list size gives
no events; unknown pixel id throws; non-contiguous ids keep their own position; duplicate pixel
definition throws; spectrum `400`, `400.5`, `4e2` accepted, `1e5`, `100`, `400x`, `table.dat`
rejected.
