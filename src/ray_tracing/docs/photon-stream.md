# Photon streams

`obdeect ray-tracing trace-stream` traces a photon stream through an optical model
and writes the camera-frame photons. The stream component (`obdeect-stream`) is
always installed with obdeect.

```bash
obdeect ray-tracing trace-stream --model telescope.json --seed 7 -i ground.ps -o camera.ps
obdeect ray-tracing trace-stream --model telescope.json --seed 7 --telescopes 1,2 < in.ps > out.ps
```

Records pass through, except that telescopes not listed in `--telescopes` are
dropped. A summary (events, photons in and out) goes to standard error. Use `-`
(the default) for standard input and output.

## Input frames

| Frame | Photons | Requirements |
| --- | --- | --- |
| G (1) | Ground plane z = 0 through the telescope sphere centre; array axes (x north, y west, z up) | Telescope pointing: `--pointing-zenith-deg` and `--pointing-azimuth-deg` (deg), or the run keys `telescope_pointing_zenith_rad` and `telescope_pointing_azimuth_rad` (rad). The azimuth is the telescope azimuth in the array frame (simtools' `phi`), not the CORSIKA shower azimuth `azimuth_rad` of the event header, which is not used. |
| T (2) | Telescope frame of `OpticalModel.trace` (+z along the optical axis, towards the sky) | Column `z` |

Other frames are rejected. The input columns `x y cx cy t w` are required.
`wavelength` is used when present, and `zem` is copied to the output; a wavelength of 0 is drawn from the
run keys `wavelength_min_nm` and `wavelength_max_nm` (as CORSIKA does). A `cz`
column must agree with `cx` and `cy`. Column `p` (acceptance already applied) is
rejected, because the stage does not reweight by it.

Each photon moves along its direction to an input plane above the incoming
obscurers, using their directional support for the telescope pointing. Its time
advances by the path length over the speed of light. The G rotation is the
ground-to-telescope rotation of simtools, so a shower that comes from the
optical axis maps to direction (0, 0, -1). For example, a telescope at zenith 20°
and array azimuth 0° is `--pointing-zenith-deg 20 --pointing-azimuth-deg 0`;
the CORSIKA shower azimuth (EVTH PHIP) is a different quantity and must not be
passed.

## Output frame (C)

| Column | Unit | Meaning |
| --- | --- | --- |
| `x`, `y` | m | Position on the focal surface, telescope frame; z = 0 is implicit |
| `cx`, `cy`, `cz` | 1 | Direction in the telescope frame; reflected photons have cz > 0 (towards the camera) |
| `t` | ns | Arrival time |
| `w` | photons | Weight (1 per survivor, or the envelope-corrected weight) |
| `wavelength` | nm | Wavelength |
| `zem` | m a.s.l. | Only if the input has it |

Axes follow the telescope frame of the optical model (see
[Optical models](optical-models.md#frames-and-units)). Chunks hold at most 2^16
rows; `first_index` counts the output rows of each telescope and event.

## Random numbers

Each (event, telescope) block draws from a seed derived from `--seed`, the event
number and the telescope id. Bunch ids are `(event_number << 32) | photon_index`,
where `photon_index` is the index of the photon in its telescope and event. The
output does not depend on the chunking or the number of threads.

## Efficiency ledger

The optical model applies `mirror_reflectivity` and, when present,
`telescope_transmission` (or structure shadowing).
Their states are carried in `efficiency_status` records of the photon stream.
The stage reads each factor's state from the ledger, run level or telescope level
(the telescope level overrides the run level):

| State | Action |
| --- | --- |
| absent or `pending` | The factor is applied in full by the trace. |
| `envelope` | The stream carries the envelope f_env, so the output weight is multiplied by 1 / f_env(λ). The trace then gives f / f_env. Photons are traced with `roulette=False` (weights min(1, w / bound) per survivor), so the weights stay unbiased; the table is interpolated at each photon's wavelength. |
| `applied` | An error: the trace would apply the factor again. |

Other factors (for example `qe`) pass on unchanged. After the run record, the
stage writes a ledger record after the incoming run statuses, marking the factors it applies
`applied`, with
`source` set to obdeect and the model hash. Input ledgers that name these
factors are rewritten the same way. For example, a run ledger with
`mirror_reflectivity.state = envelope` (table over 300 to 600 nm) becomes:

```text
factors = mirror_reflectivity,telescope_transmission
mirror_reflectivity.state = applied
mirror_reflectivity.source = obdeect 0.7.2 model <sha256>
telescope_transmission.state = applied
...
```

`telescope_end` records get the sum of the output weights of their telescope.

## Telescope transmission

`--telescope-transmission full` (default) applies the telescope transmission
t(θ), including its on-axis value t0. `--telescope-transmission angular` applies
t(θ) / t0 only, for a camera module that applies t0 itself. Per photon, the
weight is unchanged; the expected sum of output weights of an on-axis photon is
1 / t0 times that of `full`, so the camera multiplies by t0 to recover it.
`angular` needs a model with telescope transmission and a positive on-axis value.
Both modes write the same ledger keys, with the choice recorded in
`telescope_transmission.part`. In `angular` mode the ledger also lists the
factor that is left, with its on-axis value:

```text
factors = mirror_reflectivity,telescope_transmission,telescope_transmission_constant
telescope_transmission.state = applied
telescope_transmission.part = angular
telescope_transmission_constant.state = pending
telescope_transmission_constant.value = 0.9
```

A telescope-level ledger that names `telescope_transmission` gets the same
entries. A ledger that already names `telescope_transmission_constant` is
rejected in `angular` mode.

## Not applied here

Atmosphere, quantum efficiency and camera responses are not applied; they remain
for the later stages.
