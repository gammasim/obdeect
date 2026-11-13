# Python usage

```python
from obdeect.ray_tracing import DETECTED, OpticalModel

model = OpticalModel.from_simulation_models(
    "../simulation-models", "LSTN-01", "7.0.0", configuration_seed=1, zenith_deg=20.0
)
photons = model.star_photons(100_000, seed=2, towards_source=(0.0087, 0.0, 0.99996))
arrivals = model.trace(photons, seed=2)
image = arrivals["position_m"][arrivals["status"] == DETECTED]
```

`OpticalModel.from_records(records, read_file, ...)` compiles parameter records
held in memory. `save` and `load` store the model as JSON.

`trace` takes photons as arrays: `position_m` (n, 3), `direction` (n, 3, unit
vectors), `wavelength_nm`, `time_ns`, `weight` and `photon_id` (uint64; a photon's
random numbers depend only on the seed and its id). It returns, per photon:

| Array | Meaning |
| --- | --- |
| `status` | index into `STATUS_NAMES`: detected, missed_primary, missed_secondary, missed_focal_surface, blocked, invalid_input, numerical_failure |
| `terminal_surface` | id of the surface the photon ended on (focal surface, obscurer, mirror) |
| `position_m`, `direction` | at the focal surface, or where the photon was lost |
| `path_length_m`, `time_ns` | geometric path; time with the model's group index |
| `weight` | input weight times the applied responses; zero when lost |
| `primary_point_m`, `primary_surface`, `primary_incidence_deg` | reflection on the primary |
| `secondary_point_m`, `secondary_surface`, `secondary_incidence_deg` | reflection on the secondary (NaN for single-mirror telescopes) |
| `focal_incidence_deg` | angle to the focal-surface normal |

`trace` returns all 14 arrays by default. Request only what you use:

```python
detected = model.trace(
    photons,
    seed=2,
    outputs=("position_m", "direction", "time_ns", "weight", "index"),
    detected_only=True,
)
```

`outputs` lists the arrays to return; an unknown name raises `ValueError`. With
`detected_only=True`, only detected photons are returned, in input order, and `index`
gives each one's row in the input. Values match the full trace for any `threads`.

`classify_losses=False` (default `True`) reports a photon missing the primary as
`missed_primary`, even if an incoming obscurer would have blocked it first. This skips
the obscurer search; detected photons are the same either way.
