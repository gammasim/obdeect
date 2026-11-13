# Batch tracing

## Reuse worker threads

`trace(threads=...)` starts new threads on every call, which is costly for many small
blocks. A `ThreadPool` keeps its threads between calls. `ThreadPool(0)` uses every core.

```python
from obdeect.ray_tracing import ThreadPool, trace_many

pool = ThreadPool(0)
model.trace(photons, seed=2, pool=pool)

# Several models as one job; each request is (model, photons, trace options).
results = trace_many(
    [
        (model_a, photons_a, {"outputs": ("position_m", "weight", "index"), "detected_only": True}),
        (model_b, photons_b, {"apply_responses": False}),
    ],
    seed=2,
    pool=pool,
)
```

Options are those of `trace`: `apply_responses`, `direct`, `outputs`, `detected_only`
and `classify_losses`. Results match one `trace` call per request, for any pool size.
A call must not reuse its pool from inside another call on that pool.

## Shower photon bunches

`trace_bunches` thins and traces air-shower photon bunches. Sampling and output rules
are in [Air-shower photons](air-shower-photons.md).

```python
out = model.trace_bunches(
    {
        "x": x, "y": y, "cx": cx, "cy": cy,        # input plane, input length units
        "time_ns": t, "zem": zem, "size": size,    # size: photons per bunch
        "wavelength_nm": wavelength,               # 0: drawn from a density or law
        "bunch_id": bunch_id,                      # unique per seed
    },
    seed=run_seed,                # draws keyed on (seed, event, telescope, bunch)
    event=event_number,
    telescope=telescope_id,
    frame_rotation=rotation,      # 3 x 3, input frame to telescope frame
    frame_offset_m=offset_m,      # telescope frame, metres
    length_unit_m=0.01,           # metres per input length unit (cm here)
    plane_z=0.0,                  # height of the input plane, input units
    acceptance=(wavelength_nm_axis, zem_m_axis, values),  # optional, in [0, 1]
    wavelength_density=(nm, density),  # or wavelength_law=(lo_nm, hi_nm)
    normalise_transmission=False,
    pool=pool,                    # or threads=...
)
```

Give at most one of `wavelength_density` and `wavelength_law`. Invalid input, or an input
plane below an incoming obscurer, raises `ValueError` before tracing.
