# Atmosphere stage

`dlufd-propagate` (C++ tool) and `dlufd.propagate` (Python) run the atmosphere stage on a
gstrickd photon stream: sky photons are brought to the telescope plane, then every chunk
gets the atmospheric effects in a fixed order. Each effect is recorded in the efficiency
ledger of the output.

## Stages

For each telescope and event:

| # | Step | Frame | Reference |
|---|---|---|---|
| 1 | Sky photons: straight line to the horizontal plane through the telescope sphere centre; bending (refraction) unless the run says it is applied; air delay; then G | sky → G | [](refraction.md) |
| 2 | Wavelength from 1/lambda^2 for rows without one | G | [](extinction.md) |
| 3 | Extinction p_atm (plane-parallel airmass, zem = emission altitude) | G | [](extinction.md) |
| 4 | Envelope E = product of the envelope factors; thinning or weighting | G | [](extinction.md) |

Adding an effect means adding a step in `Propagator::photons` (`cpp/include/dlufd/propagator.hpp`)
and an entry in `make_ledger`. Single scattering is planned as such a step.

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
factors                 atmosphere,wavelength,qe,mirror
atmosphere.state        applied
atmosphere.source       atm_trans_2000_1_10_0_0_2000.dat
wavelength.state        applied
wavelength.range_nm     [250, 700]
qe.state                envelope
qe.source               qe.dat
qe.wavelength_nm        [200, 800]
qe.value                [0.9, 0.5]
```

- `atmosphere`: applied; source = transmission file.
- `wavelength`: applied; rows without a wavelength get one from the range.
- Each envelope factor: state `envelope` (its value is part of the thinning, not applied in
  full), with its table.

The stage refuses an input ledger that shows `atmosphere` or an envelope factor as applied
(`pending` or absent is fine). Incoming ledgers are passed through before the stage's ledger.
Readers merge the ledger records of a stream in order.

## Output

- Thinned (default): one row per surviving photon, w = 1, frame G, `first_index` counts output
  rows per telescope and event.
- Weighted: all rows with w multiplied by p; rows with w = 0 are removed.
- The column p is written only with `--store-p` (`store_p`).
- `run`: `refraction_applied` = 1 when a profile is given (sky photons are refracted).
- `telescope`: frame set to G.
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
| `--transmission FILE` | sim_telarray transmission table (required) |
| `--site-altitude M` | observation level [m a.s.l.] (required) |
| `--profile FILE` | atmospheric profile (atmprof format): air delay and refraction |
| `--envelope NAME=FILE` | envelope factor, two-column text (nm, value); repeatable |
| `--wavelength-range LO HI` | wavelengths from 1/lambda^2 [nm]; default: run keys |
| `--seed N` | Philox seed (required) |
| `--mode thinned\|weighted` | output mode (default thinned) |
| `--store-p` | keep column p |
| `-i FILE\|-`, `-o FILE\|-` | input and output stream (default stdin, stdout) |

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
)
```

`summary` is a dict with `photons_in`, `photons_out`, `dropped`, `weight_out`.
Files are paths only (no `-`).

## Limits

- Plane-parallel airmass and straight-line geometry (no Earth curvature).
- Path inside the telescope is 0 (`pathlen_tel` of sim_telarray is not set).
- Wavelengths are integer nm, as in sim_telarray.
- Envelope tables must cover the wavelength range (rounded outwards).
- The 60-degree time term of CORSIKA's ray bending is reproduced as it is; see
  [](refraction.md).
