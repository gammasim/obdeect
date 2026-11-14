# Atmosphere stage

`obdeect propagation run` (command line) and `obdeect.propagation.propagate` (Python) apply the
atmosphere to a photon stream. They move photons from their emission point to the telescope,
then apply the atmospheric effects of each chunk in a fixed order. Each effect is recorded in the
output [efficiency status](../stream/efficiency-status.md).

Use them for photons in the sky frame: air showers read with `obdeect stream from-corsika` and
the sky-frame [light sources](../sources/index.md). The output is in the ground frame (G).

## Stages

Per telescope and event, without scattering:

| # | Step | Frame | Reference |
|---|---|---|---|
| 1 | Straight line from the emission point to the horizontal plane through the telescope sphere centre; bending (refraction) unless the run says it is applied; air delay | sky → G | [](refraction.md) |
| 2 | Wavelength from 1/lambda² for rows without one | G | [](extinction.md) |
| 3 | Extinction p_atm (plane-parallel airmass by default, zem = emission altitude; curved slant depth with `--curvature`) | G | [](extinction.md), [](curvature.md) |
| 4 | Efficiency bound E = product of the `--envelope` tables; thinning or weighting | G | [](extinction.md) |

With `--scattering`, one more step runs per event, before the steps above:

| # | Step | Frame | Reference |
|---|---|---|---|
| 0 | Single scattering of the primary photons (downwards): wavelengths first, then scattered photons towards every telescope sphere, in their telescope's block | sky | [](scattering.md) |
| 1–4 | Primaries: as above (p = p_atm E) | sky → G | |
| 1s | Scattered photons: sky → G, same run geometry (refraction, air delay) | sky → G | |
| 4s | Scattered photons: p = E exp(-od_focus), thinning or weighting (own random stream, see [](../stream/random-numbers.md)) | G | |

Scattered photons skip steps 2 and 3. Their transmission from the scattering point to the
telescope sphere is already in w (`T(P → C)` of [](scattering.md)), so the atmosphere is not
applied twice. Step 4s adds only the in-telescope part exp(-od_focus(lambda)) of their telescope
(the same `od_focus` as step 4, so the path length applies to them too; with a path of 0 m the
factor is 1). The primaries are scattered before their own extinction and are otherwise unchanged.

An event is held until its `event_end` record, because the scattered photons of a telescope come
from the primaries of all telescopes of the event. Memory: one event of the input.

To add an effect, add a step in `Propagator::photons` (`cpp/include/obdeect/propagation/propagator.hpp`)
and an entry in the list of efficiency factors built there.

## Frames and positions

- Sky: x (north), y (west), z above the observation level; direction cosines (cx, cy, cz), cz < 0
  downwards. An optional cz column decides whether a photon moves upwards.
- Telescope i (1-based `telescope_id`) is at the i-th entry of `telescope_x_m`, `telescope_y_m`,
  `telescope_z_m` (z above the observation level).
- G: x, y relative to the telescope centre; z and cz are removed.
- Emission altitude a.s.l. = `observation_level_m` + z. Telescope altitude a.s.l. = the same plus
  the telescope z.

A photon in the sky frame reaches the plane if it moves downwards and starts above it. Other photons (upwards,
horizontal, emitted below the plane, or not bendable) are dropped and counted. Their rows stay in
the chunk with w = 0, so the random draws of the other photons do not change.

Time: without a profile, t + s/c. With a profile and no refraction in the stream, the air delay
`straight_path_air_delay_ns` is added. With refraction (`refraction_applied` = 0 and a profile),
time and direction come from `bend` at the telescope plane ([](refraction.md)).

## Efficiency status

One `efficiency_status` record follows each `run` record (before the first record that is not an
efficiency-status record). Example:

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

- `atmosphere`: applied; source = transmission file; `pathlen_m` = path inside each telescope [m],
  in telescope_id order (entry 0 is telescope 1).
- `wavelength`: applied; rows without a wavelength get one from the range.
- `scattering`: only with scattering; applied; source = scattering table; points per photon.
- Each efficiency bound: state `envelope` (its value is part of the thinning, not applied in full),
  with its table.

The stage refuses an input stream whose efficiency status shows `atmosphere`, `scattering` or an
efficiency bound as applied (`pending` or absent is fine). Input records are passed through before
the stage's own record. Readers merge the efficiency-status records of a stream in order. The
record format is described in [](../stream/efficiency-status.md).

## Output

- Thinned (default): one row per surviving photon, w = 1 ([](extinction.md)). `first_index` counts
  output rows per telescope and event.
- Weighted: all rows with w multiplied by p; rows with w = 0 removed.
- Column p is written only with `--store-p` (`store_p`).
- `run`: `refraction_applied` = 1 when a profile is given.
- `telescope`: frame set to G.
- Scattered photons of telescope t are in its block, after the primaries. A telescope with
  scattered photons but no input block gets its own `telescope` and `telescope_end` records
  before `event_end`.
- `telescope_end`: `photons` is the sum of w of the output.
- Other records (event, blob, end records) are passed through. The output stream record comes
  from the writer; the input one is not copied.

## Command line

```console
$ obdeect propagation run --transmission atm_trans_2000_1_10_0_0_2000.dat --site-altitude 2000 \
    --profile atmprof.dat --envelope qe=qe.dat --envelope mirror=mirror.dat \
    --wavelength-range 250 700 --seed 1 -i in.pst -o out.pst
obdeect propagation run: photons in <n>, out <m> (weight <w>), dropped <d>
```

| Option | Meaning |
|---|---|
| `--transmission FILE` | transmission table, extinction-table format (required, except with `--no-extinction`) |
| `--no-extinction` | geometry only, see below |
| `--site-altitude M` | observation level [m a.s.l.] (required) |
| `--profile FILE` | atmospheric profile (atmprof format): air delay and refraction |
| `--envelope NAME=FILE` | efficiency bound of a later factor, two-column text (nm, value); repeatable |
| `--wavelength-range LO HI` | wavelengths from 1/lambda² [nm]; default: run keys |
| `--seed N` | Philox seed (required) |
| `--mode thinned\|weighted` | output mode (default thinned) |
| `--store-p` | keep column p |
| `--scattering FILE` | single scattering with this scattering table (off by default) |
| `--scattering-points N` | scattering points per photon, 1 to 1024 (default 8; needs `--scattering`) |
| `--pathlen M` | path inside each telescope [m], finite and ≥ 0 (default 0) |
| `--telescope-pathlen ID=M` | path inside telescope ID (1-based) [m]; overrides `--pathlen`; repeatable |
| `-i FILE\|-`, `-o FILE\|-` | input and output stream (default stdin, stdout) |

Path 28 m for all telescopes, 0 m for telescope 2:

```console
$ obdeect propagation run ... --pathlen 28 --telescope-pathlen 2=0
```

The summary goes to stderr.

## Python

```python
import obdeect.propagation

summary = obdeect.propagation.propagate(
    "in.pst", "out.pst",
    transmission="atm_trans_2000_1_10_0_0_2000.dat", site_altitude_m=2000.0,
    seed=1, profile="atmprof.dat",
    envelope={"qe": "qe.dat", "mirror": "mirror.dat"},
    wavelength_range=(250.0, 700.0), mode="thinned",
    scattering="scattering_exponential.dat", scattering_points=8,  # optional
    pathlen_m=28.0, telescope_pathlen_m={2: 0.0},                  # optional
)
```

`summary` is a dict with `photons_in`, `photons_out`, `dropped`, `weight_out`. Files are paths only
(no `-`).

## Geometry only

`--no-extinction` (Python: `extinction=False`) moves photons from the sky frame to G and assigns wavelengths, and
nothing else: no extinction, envelope or thinning. Use it when a later program applies the
transmission itself, for example a light source feeding the camera simulation directly.

```console
$ obdeect propagation run --no-extinction --site-altitude 2000 --profile atmprof.dat --seed 1 \
    -i light.pst -o ground.pst
```

- Positions, directions and times are those of the normal run (same refraction and air delay).
- Photons that reach the telescope plane are kept with their weight w; photons that do not (w = 0)
  are removed, as in weighted mode. There is no column p.
- The efficiency status records `atmosphere.state` as `pending` (no source, no `pathlen_m`), so a later stage
  applies it. `wavelength` is `applied`.
- No transmission table and no `zem` are needed.
- Rejected: `--transmission`, `--envelope`, `--mode thinned`, `--scattering`, `--scattering-points`,
  `--store-p`, `--pathlen`, `--telescope-pathlen`. `--mode weighted` is accepted.

## In-telescope path

The extinction inside the telescope covers the path from the telescope sphere plane to the focus.
With path length l (`pathlen_m`), the optical depth is tau(h_T + l) − tau(h_T), with h_T the
telescope altitude. `--pathlen` sets l for all telescopes; `--telescope-pathlen`
(`telescope_pathlen_m`, 1-based ids) overrides it for single telescopes. Default 0, so existing runs
keep their result. Intended l per telescope: the focal length for single-mirror telescopes,
`|p0-s0|+|s0-f0|` for dual-mirror ones (primary, secondary, focus). An id that the run does not
have is an error.

## Limits

- Plane-parallel airmass and straight-line geometry by default. Earth curvature for the airmass is
  optional ([](curvature.md)); curved refraction and curved air delay are not implemented.
- Scattering is plane-parallel and needs sky-frame input and `telescope_r_m` in the run keys; not
  available for photons already in frame G, and rejected with a curved mode.
- The in-telescope path is not derived from the optics (default 0).
- Transmission tables hold one row per integer nm; a wavelength uses the row at its integer part.
- Efficiency-bound tables must cover the wavelength range (rounded outwards).
- Photons within about 3 degrees of the horizon are not bent and are dropped (series limit of the
  ray bending, see [](refraction.md)).
