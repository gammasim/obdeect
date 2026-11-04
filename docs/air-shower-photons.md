# Air-shower photons

`trace_bunches` traces photon bunches from an air shower through an optical model.
It covers the optics only: the photons are traced from the input plane to the
camera with the model's responses. The atmosphere and the quantum efficiency
are not applied here. The caller applies them; the simtools shower route and
sim_telarray do this in two steps (see [Transmission convention](#transmission-convention)).

## Inputs

Each bunch is `size` photons (fractional sizes allowed) that share one position,
direction, time, emission height and wavelength. The columns are:

| Column | Unit | Meaning |
| --- | --- | --- |
| `x`, `y` | input length units | Position on the input plane |
| `cx`, `cy` | 1 | Direction cosines; photons move towards -z of the input frame, `cz = -sqrt(1 - cx^2 - cy^2)` |
| `time_ns` | ns | Time of the bunch at the input plane |
| `zem` | input length units | Emission height (used by the acceptance table only) |
| `size` | photons | Number of photons; `ceil(size)` sub-photons, the last one with fraction `size - floor(size)` |
| `wavelength_nm` | nm | Wavelength; 0 draws the wavelength from a density or a law |
| `bunch_id` | uint64 | Unique per bunch for one seed |

The frame maps the input to the telescope frame: `p = R (u x, u y, u plane_z) + offset`
and `d = R (cx, cy, cz)`, with `u` the metres per input length unit and `R` orthonormal.

## Wavelength

- A bunch with `wavelength_nm > 0` keeps that wavelength.
- A bunch with `wavelength_nm = 0` draws its wavelength per photon from either:
  - `wavelength_density`: a piecewise-linear density over nm, sampled by inverse transform; or
  - `wavelength_law = (lo_nm, hi_nm)`: a density proportional to 1/lambda^2, sampled
    exactly, `lambda = 1 / (1/lo - u (1/lo - 1/hi))`, and truncated to an integer
    number of nm, `floor(lambda)`, as sim_telarray does. Values lie in
    `[floor(lo), hi)`, with `1 <= lo < hi` so that no value is 0 nm.

Giving both raises `ValueError`.

## Thinning

For each sub-photon the keep probability is

`p_keep = fraction x acceptance(lambda, zem) x bound(lambda)`

- `fraction` is 1, or `size - floor(size)` for the last sub-photon of a fractional bunch.
- `acceptance` is an optional table over (wavelength [nm], emission height [m]),
  interpolated bilinearly, with values in [0, 1]. Without a table it is 1.
- `bound` is an upper bound of the model's responses at that wavelength: the
  maximum reflectivity over incidence angle for the primary and secondary mirrors,
  the maximum degradation factor, and the on-axis value of the telescope
  transmission. It is multiplied by `1 + 1e-12`.

A sub-photon is kept if `draw < p_keep`, where `draw` is the first number of the
Philox stream of `(seed, bunch_id, sub_photon)`. The survivors are recorded
before tracing, so the output is sized to the survivors.

## Roulette

Each survivor is traced with the responses. A detected photon with weight `w`
is kept with probability `min(1, (w / norm) / bound)` (`norm = 1`, or the
on-axis transmission with `normalise_transmission=True`), using the second number of
the same Philox draw. With `roulette=True` (default) a kept photon has weight 1.
With `roulette=False` all detected survivors stay, with weight `min(1, w / norm / bound)`.
The expected weight of a bunch equals the sum of the weights of a full trace.

## Transmission convention

The telescope transmission is a response of the model. With the default
`normalise_transmission=False` the output weights carry it in absolute terms.
The shower route in simtools uses this setting. Under `Bypass_Optics` sim_telarray
takes the airmass from the camera-frame direction of each bunch, not from the
sky direction. So the route multiplies each photon weight by the ratio of the
native transmission to the transmission sim_telarray applies, computed with
sim_telarray's formula at the ground airmass and the bypass airmass. sim_telarray
then applies the rest (the bypass transmission and the quantum efficiency).
The total equals the native chain. The correction and its validation are in
[Validation](VALIDATION.md#air-shower-photons-phase-3).

## Time reference

Output `time_ns` is the input `time_ns` plus the path from the input plane through
the optics, times the model's `propagation_group_index`, divided by the speed of light.
obdeect adds no delay for the telescope height or the airmass.

## Input plane

The input plane must lie above every incoming obscurer. `trace_bunches` checks
this before tracing and raises `ValueError` if it does not (see
[Batch tracing](batch-tracing.md#shower-photon-bunches)). The check includes the
tilt from `frame_rotation` and the axes offset of the model.

## Determinism

The random numbers of a sub-photon depend only on `(seed, bunch_id, sub_photon)`.
The output does not depend on the number of threads, the chunking or the order of
the bunches. Bunch ids must be unique per seed.

## Not modelled

- Atmosphere (absorption, airmass) and quantum efficiency. Applying them inside obdeect is later work.
- Cherenkov emission, shower development and the CORSIKA input format. These are in simtools.
- Electronics, pixel readout, triggers, and the event output.
- Scattering or absorption that the model does not describe as a response.
