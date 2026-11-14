# Source model

A source emits photons from a point (an emitter) or as a parallel beam (a star). Photons
are written as photon stream chunks in the frame of the source.

## Model

| Part | Choices |
| --- | --- |
| Geometry | `emitter`: position, pointing, angular profile. `star`: parallel beam (infinite distance) or point source at finite distance, on a disc |
| Angular profile | `isotropic`, `cone` (one angle), `filled_cone`, `lambertian` (cos theta), `gauss` (sigma in deg), `rayleigh`, `table` (angle in deg, value) |
| Spectrum | `fixed` (nm), `uniform`, `inverse_square` (1/lambda^2, Cherenkov), `table` |
| Pulse | `fixed` (ns), `tophat` (full width, delay), `gauss` (sigma, mean), `exponential` (decay, mean), `complex` (flat, Gauss r.m.s., decay, delay), `table` |
| Count | `photons`: mean number emitted into 4 pi (emitter) or on the disc (star) |
| Target | optional, emitters only: `disc` (weighted) or `sphere` (cone) |

Profiles are relative: an emitter with profile `P` emits all `photons` with directions
distributed as `P`. Times are emission times at the source. A star's times are the wavefront
arrival at the disc centre ([ray-tracing convention](../ray-tracing/light-sources.md)).

```python
from obdeect.sources.source import emit

flasher = {
    "frame": "camera",
    "photons": 1e6,
    "emitter": {
        "position_m": [0.0, 0.0, 1.5],
        "pointing": [0.0, 0.0, -1.0],
        "angular": {"kind": "filled_cone", "angle_deg": 11.0},
    },
    "spectrum": {"kind": "fixed", "wavelength_nm": 385.0},
    "pulse": {"kind": "tophat", "width_ns": 2.5},
    "target": {"kind": "disc", "centre_m": [0.0, 0.0, 0.0],
               "normal": [0.0, 0.0, 1.0], "radius_m": 1.2},
}
photons = emit(flasher, seed=3)   # float32 arrays: x y z cx cy cz t w wavelength
```

`iter_chunks` emits large counts (at most 2^16 rows per chunk); `emit_chunk(..., first_index=,
rows=)` emits any slice.

## Frames and units

- Frames: `sky` (absolute), `telescope` (T, z along the optical axis), `camera` (C). Ground is
  not a source frame. The frame is stored in the chunk.
- Lengths in metres, times in ns, wavelengths in nm, angles in degrees (the core takes
  `u = 1 - cos(theta)`).
- Directions are unit vectors of the propagation direction. `cz < 0` is downward, as in
  [ray tracing](../ray-tracing/light-sources.md) and [stream format](../stream/format.md).
- Without a target, `w = 1` for every photon.

## Targets (importance sampling)

Most photons of a calibration lamp miss the telescope. A target generates only photons that
can reach it.

- `disc`: candidates are uniform points on the disc. Weight
  `w = I(theta) / I_max (h / r)^3`, the probability that the emitted photon reaches the disc
  (`h` distance from emitter to disc plane, `r` to the point). Candidate count
  `photons R^2 / (4 h^2)` for an isotropic emitter (general: `I_max / (2 pi integral(I du)) x
  pi R^2 / h^2`). Sum of weights = expected photons on the disc. Efficient when the profile
  covers the disc; narrow cones off the disc give mostly `w = 0`.
- `sphere`: directions sampled in the cone containing the sphere, `w = 1`, count
  `photons (1 - cos theta_max) / 2`. Isotropic emitters only.

## Conventions

- Emitter count: `floor(mean)` plus one with probability equal to the fraction (not Poisson).
  The mean is the same; weights are not thinned.
- Units are SI-like (m, ns, nm, deg). Every photon has weight 1, except disc targets and
  bunches (`bunch_size` photons per row in the tools).
- Gaussian pulse parameters are sigma: `sigma = FWHM / 2.35482`.
- Tool spellings map to dictionary kinds: `Cone:` = `cone`, `FilledCone:` = `filled_cone`,
  `Parallel` = `cone` at 0 deg, `Gauss-rms:` = `gauss`, `Rayleigh` = `rayleigh`,
  `Cherenkov:a:b` = `inverse_square`; `Simple:`, `Gauss:`, `Exp:`, `Complex:` are pulse kinds.
- Tables are passed as arrays, not read from files.
- Refraction and airlight belong to [obdeect propagation](../propagation/index.md). The
  Gaussian spectrum is not available.

Random numbers: [](../stream/random-numbers.md).
