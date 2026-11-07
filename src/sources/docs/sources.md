# Light sources

A source emits photons from a point (an emitter) or as a parallel beam (a star).
Its photons are written as gstrickd photon chunks, in the frame of the source.

## Model

| Part | Choices |
| --- | --- |
| Geometry | `emitter`: position, unit `pointing`, angular profile. `star`: parallel beam (infinite distance) or point source at a finite distance, on a disc |
| Angular profile | `isotropic`, `cone` (one angle), `filled_cone`, `lambertian` (cos theta), `gauss` (sigma in deg), `rayleigh`, `table` (angle in deg, value) |
| Spectrum | `fixed` (nm), `uniform`, `inverse_square` (1/lambda^2, Cherenkov), `table` |
| Pulse | `fixed` (ns), `tophat` (full width, delay), `gauss` (sigma, mean), `exponential` (decay, mean), `complex` (flat, Gauss r.m.s., decay, delay), `table` |
| Count | `photons`: mean number of photons emitted into 4 pi (emitter) or on the disc (star) |
| Target | optional, emitters only: `disc` (weighted) or `sphere` (cone), see below |

Profiles are relative: an emitter with profile `P` emits all `photons` with the
directions distributed as `P`. Times are emission times at the source; a star's
times are the wavefront arrival at the disc centre (obdeect convention).

```python
from dlamb.source import emit

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

Use `iter_chunks` for large counts (at most 2^16 rows per chunk) and
`emit_chunk(..., first_index=, rows=)` to emit any slice of the photons.

## Frames and units

- Frames: `sky` (absolute), `telescope` (T, telescope-local, z along the optical axis),
  `camera` (C). Ground (G) is not a source frame. The frame is stored in the chunk.
- Lengths in metres, times in ns, wavelengths in nm, angles in degrees (in the
  dictionary; the core takes `u = 1 - cos(theta)`).
- Directions are unit vectors of the propagation direction. `cz < 0` is downward,
  as in obdeect and gstrickd. Without a target, `w = 1` for every photon.

## Targets (importance sampling)

Most photons of a calibration lamp miss the telescope. A target generates only
photons that can reach it:

- `disc` (obdeect `Emitter` method): candidates are uniform points on the disc. Each
  candidate has weight `w = I(theta) / I_max * (h / r)^3`, the probability that the
  emitted photon reaches the disc (`h` distance from emitter to the disc plane, `r`
  to the point). The number of candidates is `photons * R^2 / (4 h^2)` for an isotropic
  emitter (general: `I_max / (2 pi integral(I du)) * pi R^2 / h^2`). Sum of weights
  = expected photons on the disc. Efficient when the profile covers the disc.
  Narrow cones off the disc give mostly `w = 0`.
- `sphere` (LightEmission style): directions sampled in the cone that contains the
  sphere, `w = 1`, count `photons * (1 - cos theta_max) / 2`. Isotropic emitters only.

## Random numbers

Each photon draws from Philox (gstrickd `uniform4`) keyed on
`(seed, stage, event, telescope, photon index, stream)`, with stage `0x444C0001` and
streams in `cpp/include/dlamb/streams.hpp`. Results do not depend on chunk size or
order. Two sources with the same seed, event and telescope use the same draws: give
each independent source its own seed.

## Differences to obdeect

- Same distributions (star and emitter formulas ported; moment tests cite obdeect),
  but the random numbers differ: obdeect keys on `(seed, photon id)`.
- The emitter count is `floor(mean)` plus one with probability equal to the fraction
  (LightEmission rounding), not a Poisson draw. The mean is the same; the weights are
  not thinned, so there is no Bernoulli step.
- `sphere` target is new. Photons are gstrickd columns, not obdeect `Ray` objects.

## Differences to LightEmission

- Units are SI-like (m, ns, nm, deg), not cm and bunches of `bsize` photons.
  Every photon has weight 1 except for disc targets.
- `photons` is the total number of emitted photons; LightEmission's fiducial-sphere
  fraction is replaced by an explicit target.
- Gaussian pulse parameters are sigma: convert FWHM with `sigma = FWHM / 2.35482`.
- `Cone:` is `cone`, `FilledCone:` is `filled_cone`, `Parallel` is `cone` with 0 deg,
  `Gauss-rms:` is `gauss`, `Rayleigh` is `rayleigh`, `Cherenkov:a:b` is
  `inverse_square`, and `Simple:`, `Gauss:`, `Exp:`, `Complex:` are the pulse kinds.
- Tables are passed as arrays, not read from files.
- Not ported yet: `BackgroundLightSource` (night-sky background), moving or line
  sources (`length`, `beta`), the `Run` geometry (fiducial sphere, observation
  level), refraction and airlight. The Gaussian spectrum is missing.
- The IACT bunch writes (`telout_`, block 1205) are the adapter's job. LightEmission
  stores `w = -cz`; the mapping to gstrickd is not checked here.
