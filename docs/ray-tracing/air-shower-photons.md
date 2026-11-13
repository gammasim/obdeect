# Air-shower photons

`trace_bunches` traces photon bunches from an air shower from the input plane to the
camera entrance, with the model's responses. It does not apply the atmosphere or the
quantum efficiency; the caller applies them (see [Transmission convention](#transmission-convention)).
Arguments and example: [Batch tracing](batch-tracing.md#shower-photon-bunches).

## Inputs

Each bunch is `size` photons (fractional sizes allowed) sharing one position,
direction, time, emission height and wavelength.

| Column | Unit | Meaning |
| --- | --- | --- |
| `x`, `y` | input length units | Position on the input plane |
| `cx`, `cy` | 1 | Direction cosines; photons move towards -z, `cz = -sqrt(1 - cx^2 - cy^2)` |
| `time_ns` | ns | Time of the bunch at the input plane |
| `zem` | input length units | Emission height (used by the acceptance table only) |
| `size` | photons | `ceil(size)` sub-photons; the last one has fraction `size - floor(size)` |
| `wavelength_nm` | nm | Wavelength; 0 draws it from a density or a law |
| `bunch_id` | uint64 | Unique per bunch for one seed |

The frame maps input to telescope coordinates: `p = R (u x, u y, u plane_z) + offset`
and `d = R (cx, cy, cz)`, with `u` metres per input length unit and `R` orthonormal.

## Wavelength

- A bunch with `wavelength_nm > 0` keeps that wavelength.
- A bunch with `wavelength_nm = 0` draws one per photon from either:
  - `wavelength_density`: a piecewise-linear density over nm, by inverse transform; or
  - `wavelength_law = (lo_nm, hi_nm)`: density proportional to 1/lambda^2, sampled
    exactly as `lambda = 1 / (1/lo - u (1/lo - 1/hi))`, then truncated to `floor(lambda)`
    as sim_telarray does. Values lie in `[floor(lo), hi)`; the law needs `1 <= lo < hi`,
    so no value is 0 nm.

Giving both raises `ValueError`.

## Thinning

For each sub-photon the keep probability is

`p_keep = fraction x acceptance(lambda, zem) x bound(lambda)`

- `fraction` is 1, or `size - floor(size)` for the last sub-photon of a fractional bunch.
- `acceptance` is an optional table over (wavelength [nm], emission height [m]),
  bilinear, values in [0, 1]. Without a table it is 1.
- `bound` is an upper bound of the responses at that wavelength: the maximum
  reflectivity over incidence angle (primary and secondary), the maximum degradation
  factor, and the on-axis telescope transmission. It is multiplied by `1 + 1e-12`.

A sub-photon is kept if `draw < p_keep`, where `draw` is the first number of the Philox
stream of `(seed, bunch_id, sub_photon)`. Survivors are recorded before tracing.

## Roulette

Each survivor is traced with the responses. A detected photon of weight `w` is kept with
probability `min(1, (w / norm) / bound)`, using the second number of the same draw.

- `roulette=True` (default): kept photons have weight 1.
- `roulette=False`: all detected survivors stay, with weight `min(1, w / norm / bound)`.

`norm` is 1, or the on-axis transmission `t0` with `normalise_transmission=True`. That
setting divides the bound and the roulette weight by `t0`, so output weights are
relative to `t0`; the caller then applies `t0`. If the model has no transmission, or
`t0` is not positive, the divisor is 1. The expected weight of a bunch equals the sum of
weights of a full trace.

## Transmission convention

The telescope transmission is a response of the model. With the default
`normalise_transmission=False`, the output weights carry it in absolute terms. The
simtools shower route uses this setting. Under `Bypass_Optics`, sim_telarray takes the
airmass from the camera-frame direction of each bunch, not from the sky direction. So the
route multiplies each weight by the ratio of the native transmission to the transmission
sim_telarray applies, computed with sim_telarray's formula at the ground airmass and the
bypass airmass. sim_telarray then applies the rest (bypass transmission and quantum
efficiency). The total equals the native chain. See
[Validation](VALIDATION.md#air-shower-photons-phase-3).

## Input plane

The input plane must lie above every incoming obscurer, such as the camera housing.
Photons starting below an obscurer never meet it. The plane is `n . p = value` with
normal `n = R e_z`. Before tracing, `trace_bunches` compares it with each incoming
obscurer's highest point along `n`. If the plane passes below that point, it raises
`ValueError` naming the obscurer and the height along the telescope axis, `value / n_z`.
Raise `plane_z` to fix it.

`frame_offset_m` is the input origin in the telescope frame. For a model with an axes
offset (see [Optical models](optical-models.md#axes-offset)), the native tracer subtracts
that offset. A ground-frame point P is then traced at `R P + frame_offset_m - O`, which is
`R (P - O_ground)` for the ground-frame offset `O_ground` of the sim_telarray
`axes_offsets` record. This assumes `frame_rotation` is the ground-to-telescope rotation
at the compiled zenith angle.

## Time reference

Output `time_ns` is the input `time_ns` plus the path from the input plane through the
optics, times the model's `propagation_group_index`, divided by the speed of light.
No delay is added for telescope height or airmass.

## Determinism

The random numbers of a sub-photon depend only on `(seed, bunch_id, sub_photon)`. The
output does not depend on threads, chunking or bunch order. Bunch ids must be unique per
seed. The recommended id is `(event_number << 32) | bunch_index_in_event`, with a seed per
telescope derived from the run seed and the telescope number.

## Not modelled

- Atmosphere (absorption, airmass) and quantum efficiency. Applying them inside obdeect
  ray tracing is later work.
- Cherenkov emission, shower development and the CORSIKA input format (in simtools).
- Electronics, pixel readout, triggers and event output.
- Scattering or absorption that the model does not describe as a response.
