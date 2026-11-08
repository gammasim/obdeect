# Atmosphere stage

`dlufd-propagate` (C++ tool) and `dlufd.propagate` (Python) run the atmosphere stage on a
gstrickd photon stream: sky photons are brought to the telescope plane, then every chunk
gets the atmospheric effects in a fixed order. Each effect is recorded in the efficiency
ledger of the output.

## Stages

For each telescope and event, without scattering:

| # | Step | Frame | Reference |
|---|---|---|---|
| 1 | Sky photons: straight line to the horizontal plane through the telescope sphere centre; bending (refraction) unless the run says it is applied; air delay; then G | sky → G | [](refraction.md) |
| 2 | Wavelength from 1/lambda^2 for rows without one | G | [](extinction.md) |
| 3 | Extinction p_atm (plane-parallel airmass, zem = emission altitude) | G | [](extinction.md) |
| 4 | Envelope E = product of the envelope factors; thinning or weighting | G | [](extinction.md) |

With `--scattering` the stage runs one more step per event, before the steps above:

| # | Step | Frame | Reference |
|---|---|---|---|
| 0 | Single scattering of the primary photons (sky, downwards): wavelengths assigned first, then scattered photons towards every telescope sphere, written in their telescope's block | sky | [](scattering.md) |
| 1 | Primaries: as above | sky → G | |
| 2, 3 | Primaries: wavelength, extinction p_atm | G | |
| 4 | Primaries: p = p_atm E, thinning or weighting | G | |
| 1s | Scattered photons: sky → G, same run geometry (refraction, air delay) | sky → G | |
| 4s | Scattered photons: p = E exp(-od_focus), thinning or weighting (their own random stream) | G | |

Scattered photons skip steps 2 and 3: their transmission from the scattering point to the telescope
sphere is already in w (`T(P -> C)` of [](scattering.md)), so the atmosphere is not applied twice.
Step 4s adds only the in-telescope part of the extinction, exp(-od_focus(lambda)) of their telescope
(the same `od_focus` as step 4, so the path length applies to scattered photons too). With a path of
0 m the factor is 1. The primaries are scattered before their own extinction and are otherwise
unchanged.

An event is held until its `event_end` record, because the scattered photons of a telescope come
from the primaries of all telescopes of the event. Memory: one event of the input.

Adding an effect means adding a step in `Propagator::photons` (`cpp/include/dlufd/propagator.hpp`)
and an entry in `make_ledger`.

## Frames and positions

- Sky: x (north), y (west), z above the observation level, direction cosines (cx, cy, cz),
  cz < 0 downwards. A cz column is optional; if present its sign decides whether a photon
  moves upwards.
- Telescope i (1-based `telescope_id`) is at the i-th entry of the run arrays
  `telescope_x_m`, `telescope_y_m`, `telescope_z_m` (z above the observation level).
- G: x, y relative to the telescope centre; z and cz are removed.
- Emission altitude a.s.l. = `observation_level_m` + z. Telescope altitude a.s.l. = the same
  plus the telescope z.

A sky photon reaches the plane if it moves downwards and starts above it. Other photons
(upwards, horizontal, emitted below the plane, or not bendable) are dropped and counted.
Their rows stay in the chunk with w = 0, so the random draws of the other photons do not
change.

Time: without a profile, t + s/c. With a profile and no refraction in the stream, the air
delay `straight_path_air_delay_ns` is added. With refraction (`refraction_applied` = 0 and a
profile), the time and direction come from `bend` at the telescope plane.

## Ledger

One `ledger` record follows each `run` record (before the first record that is not a ledger).
Example:

```text
factors                 atmosphere,wavelength,scattering,qe
atmosphere.state        applied
atmosphere.source       atm_trans_2000_1_10_0_0_2000.dat
atmosphere.pathlen_m    [28, 0]
wavelength.state        applied
wavelength.range_nm     [250, 700]
scattering.state        applied
scattering.source       scattering_exponential.dat
scattering.points_per_photon 8
qe.state                envelope
qe.source               qe.dat
qe.wavelength_nm        [200, 800]
qe.value                [0.9, 0.5]
```

- `atmosphere`: applied; source = transmission file; `pathlen_m` = path inside each telescope [m]
  in telescope_id order (entry 0 is telescope 1), as used for the in-telescope extinction.
- `wavelength`: applied; rows without a wavelength get one from the range.
- `scattering`: only with scattering; applied; source = scattering table; points per photon.
- Each envelope factor: state `envelope` (its value is part of the thinning, not applied in
  full), with its table.

The stage refuses an input ledger that shows `atmosphere`, `scattering` or an envelope factor as
applied (`pending` or absent is fine). Incoming ledgers are passed through before the stage's ledger.
Readers merge the ledger records of a stream in order.

## Output

- Thinned (default): one row per surviving photon, w = 1, frame G, `first_index` counts output
  rows per telescope and event.
- Weighted: all rows with w multiplied by p; rows with w = 0 are removed.
- The column p is written only with `--store-p` (`store_p`).
- `run`: `refraction_applied` = 1 when a profile is given (sky photons are refracted).
- `telescope`: frame set to G.
- Scattered photons of telescope t are in its block, after the primaries; a telescope with
  scattered photons but no input block gets a `telescope` and `telescope_end` record of its own
  before the `event_end` record.
- `telescope_end`: `photons` is the sum of w of the output (not the input).
- Other records (event, blob, end records) are passed through.
- The output stream record comes from the writer; the input stream record is not copied.

## Command line

```console
$ dlufd-propagate --transmission atm_trans_2000_1_10_0_0_2000.dat --site-altitude 2000 \
    --profile atmprof.dat --envelope qe=qe.dat --envelope mirror=mirror.dat \
    --wavelength-range 250 700 --seed 1 -i in.ps -o out.ps
dlufd-propagate: photons in <n>, out <m> (weight <w>), dropped <d>
```

| Option | Meaning |
|---|---|
| `--transmission FILE` | transmission table, extinction-table format (required) |
| `--site-altitude M` | observation level [m a.s.l.] (required) |
| `--profile FILE` | atmospheric profile (atmprof format): air delay and refraction |
| `--envelope NAME=FILE` | envelope factor, two-column text (nm, value); repeatable |
| `--wavelength-range LO HI` | wavelengths from 1/lambda^2 [nm]; default: run keys |
| `--seed N` | Philox seed (required) |
| `--mode thinned\|weighted` | output mode (default thinned) |
| `--store-p` | keep column p |
| `--scattering FILE` | single scattering with this scattering table (off by default) |
| `--scattering-points N` | scattering points per photon, 1 to 1024 (default 8; needs `--scattering`) |
| `--pathlen M` | path inside each telescope [m], finite and >= 0 (default 0) |
| `--telescope-pathlen ID=M` | path inside telescope ID (1-based) [m]; overrides `--pathlen`; repeatable |
| `-i FILE\|-`, `-o FILE\|-` | input and output stream (default stdin, stdout) |

Example with the in-telescope path of 28 m for all telescopes, 0 m for telescope 2:

```console
$ dlufd-propagate ... --pathlen 28 --telescope-pathlen 2=0
```

The summary goes to stderr.

## Python

```python
import dlufd

summary = dlufd.propagate(
    "in.ps", "out.ps",
    transmission="atm_trans_2000_1_10_0_0_2000.dat", site_altitude_m=2000.0,
    seed=1, profile="atmprof.dat",
    envelope={"qe": "qe.dat", "mirror": "mirror.dat"},
    wavelength_range=(250.0, 700.0), mode="thinned",
    scattering="scattering_exponential.dat", scattering_points=8,  # optional
    pathlen_m=28.0, telescope_pathlen_m={2: 0.0},                  # optional
)
```

`summary` is a dict with `photons_in`, `photons_out`, `dropped`, `weight_out`.
Files are paths only (no `-`).

## In-telescope path

The extinction inside the telescope covers the path from the telescope sphere plane to the focus.
With path length l (`pathlen_m`), the optical depth is tau(h_T + l) - tau(h_T), with h_T the
telescope altitude. `--pathlen` (`pathlen_m`) sets l for all telescopes, `--telescope-pathlen`
(`telescope_pathlen_m`, keys 1-based telescope ids) overrides it for single telescopes. Default 0,
so existing runs keep their result. Each telescope's l comes from its optics: the focal length for
single-mirror telescopes, `|p0-s0|+|s0-f0|` for dual-mirror ones (primary, secondary, focus).
An id that the run does not have is an error.

## Limits

- Plane-parallel airmass and straight-line geometry (no Earth curvature).
- Scattering needs sky-frame input and `telescope_r_m` in the run keys; it is not available for
  photons already in frame G.
- The in-telescope path is not derived from the optics: it is given per telescope (default 0).
- Wavelengths are integer nm.
- Envelope tables must cover the wavelength range (rounded outwards).
- Photons within about 3 degrees of the horizon are not bent and are dropped (series limit of
  the ray bending, see [](refraction.md)).
