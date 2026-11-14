# led: specification

Options and examples: [led](../programs/led.md). Frame C. Writer A, producer
`obdeect sources led`.

## Validation (prefix `led: `)

- Missing camera file: `no camera definition file defined (--camera-file)`; unopenable:
  `cannot open camera configuration file <path>`.
- Exactly one of `--use`, `--use-trg-rand` (> 0), `--use-trg-all`; else `give exactly one of ...`.
- At least one pixel in the camera file. `events >= 0`, `require >= 0`, `photons >= 0`,
  `bunch > 0`, `radius > 0`.
- `--random N` sets random mode and `require = N`. `--use-trg-rand 0` is the default.
- Angular profile is fixed to `Parallel`.

## Camera file

- Line-oriented; leading blanks ignored; `#` and blank lines skipped.
- Keywords (case-insensitive, first word): `Pixel <id> <type> <x_cm> <y_cm> ...` defines a
  pixel; `Rotate` and `Pixtype` are skipped. A line that does not parse as `int int double
  double` prints `led: problem with pixel definition: <line>` to stderr and is skipped.
- Negative or duplicate pixel id: error.
- Any other line containing `Trigger` (case-sensitive) is a trigger group.
- Trigger text (`clean_trigger_list(line, keep_slaves)`): if the line contains `of`, use the text
  after the first `of`; unless `keep_slaves`, drop every `[...]` group. The pixel file is read
  after all options, so `--slaves` decides.

## Pixel list syntax

- Items separated by blanks or tabs: `[+]id`, optionally followed by `[a,b,..]` (slaves).
  `%` or newline ends the list. `+` marks a required pixel.
- Errors: `led: invalid in list of pixels to use: '<rest>'`; `led: pre-sum pixel group not
  terminated in '<text>'`; a pixel not in the file: `led: pixel ID <n> is not defined in the
  camera file`.
- A repeated id is skipped with `duplicated pixel <id> in use list, used only once` on stderr.

## Lists and events

- `--use-trg-all`: one list per trigger group. No groups: `the camera definition has no
  trigger groups`.
- `--use-trg-rand N`: N lists; group `j = floor(G program_uniform(seed, 0, draw++))`, repeated
  until `j < G`. The draw counter increases across the N draws.
- `--use`: one list. An empty list: `led: failed to obtain a pixel list`.
- Per list, `nreq = require`, or the list size if `require = 0`. If `nreq` exceeds the list
  size: `led: required <n> pixels but the list has <m>; list skipped` (stderr), list skipped.
- Otherwise repeat `events` times:
  - `nreq` equals list size: one event with all pixels.
  - `--random`: one event with `nreq` distinct pixels. Draw `p = floor(L program_uniform(seed,
    event_number, draw++))` (counter restarts per event); ignore `p >= L` and used ids.
  - else: one event per combination of `nreq` pixels in lexicographic order of list position,
    starting from the `+` pixels, then non-required pixels in increasing list order. The whole
    combination set is emitted in each of the `events` repetitions.
- Event numbers are consecutive from 1 across all lists.

## Emissions

- Each lit pixel (and, with `--slaves`, each slave, without labels) is one emission, in pixel
  order. Emitter at `(x_cm/100, y_cm/100, 0)` m, pointing `(0, 0, -1)`, `Source.photons =
  photons / b`, `b = max(1, bunchsize)`, seed `derive_seed(seed, j)` over a running emission
  counter, telescope key 1, frame C, run radius `R = camera_radius / 100`.
- FixedAngle: never a target. Photons are written at `z = 0`, `cx = cy = 0`, `cz = -1`; `t`
  from the pulse. Pixel positions use pixel-file coordinates (no y flip).

## `--list-only`

For every event print `(@) ` followed by, for each label in ascending numeric order, a blank
and the id (`(@)  0 1`), then a newline. Labels are the listed pixels only (never slaves). No
stream is written.
