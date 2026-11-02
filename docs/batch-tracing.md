# Batch tracing

## Reuse worker threads

`trace(threads=...)` starts its threads on every call, which is costly when an
event is traced as many small blocks. A `ThreadPool` keeps its threads between
calls: `ThreadPool(n)` (0 uses every core) runs one job at a time, and
`model.trace(photons, seed=2, pool=pool)` uses it instead of `threads`.
`trace_many` traces the photons of several optical models, as one job, in a
single call:

```python
from obdeect import ThreadPool, trace_many

pool = ThreadPool(0)
results = trace_many(
    [
        (model_a, photons_a, {"outputs": ("position_m", "weight", "index"), "detected_only": True}),
        (model_b, photons_b, {"apply_responses": False}),
    ],
    seed=2,
    pool=pool,
)
```

Each request is `(model, photons, options)` with the options `apply_responses`,
`direct`, `outputs`, `detected_only` and `classify_losses` of `trace`; the result is the list of
the dicts that one `trace` call per request returns, with the same values for
any pool size. A call must not use the same pool from inside another call that
runs on it.

## Shower photon bunches

Air-shower photons come as bunches: `size` photons (may be fractional) with one
position, direction, time, emission height and wavelength. `trace_bunches`
thins and traces them without making every photon a full trace:

```python
out = model.trace_bunches(
    {
        "x": x,
        "y": y,
        "cx": cx,
        "cy": cy,
        "time_ns": t,
        "zem": zem,
        "size": size,
        "wavelength_nm": wavelength,
        "bunch_id": bunch_id,
    },
    seed=run_seed_of_this_telescope,
    frame_rotation=rotation,  # 3 x 3, input frame to telescope frame
    frame_offset_m=offset_m,  # telescope frame; the model's axes offset is subtracted
    length_unit_m=0.01,  # metres per input length unit (cm here)
    plane_z=0.0,  # height of the input plane, input units; above all incoming obscurers
    acceptance=(wavelength_nm_axis, zem_m_axis, values),  # optional, in [0, 1]
    wavelength_density=(nm, density),  # optional; or wavelength_law=(lo_nm, hi_nm)
    normalise_transmission=False,  # optional, see below
    pool=pool,  # or threads=...
)
```

The inputs are float32 arrays (`bunch_id` is uint64). Photons move towards -z of
the input frame: `cz = -sqrt(1 - cx^2 - cy^2)`. Position is
`R (u x, u y, u plane_z) + offset` and direction `R (cx, cy, cz)`, with
`u = length_unit_m`; `R` must be orthonormal. A bunch of `wavelength_nm = 0`
takes a wavelength per photon from `wavelength_density` (linear between the
nodes, over nm).

Alternatively, `wavelength_law=(lo_nm, hi_nm)` draws the wavelength of those
bunches with a density proportional to 1/lambda^2 on [lo_nm, hi_nm], by exact
inverse transform. This is the sampling sim_telarray uses over the CORSIKA
header range. The wavelength is then truncated to an integer number of nm,
`floor(lambda)`, as in sim_telarray. The values lie in [floor(lo_nm), hi_nm),
so for a non-integer `lo_nm` the lowest value can be below `lo_nm`. The law
needs finite 0 < lo_nm < hi_nm. Give at most one of `wavelength_density` and
`wavelength_law`; giving both raises `ValueError`. Bunches with a fixed
wavelength are not affected.

`frame_offset_m` is the input origin in the telescope frame, in metres. A
model with an axes offset (see [Optical models](optical-models.md)) stores it as
`axes_offset_m`, and `trace_bunches` subtracts that vector from `frame_offset_m`.
A ground-frame point P is therefore traced at `R P + frame_offset_m - O`, which
is `R (P - O_ground)` for the ground-frame offset `O_ground` of the
sim_telarray `axes_offsets` record. This assumes that `frame_rotation` is the
ground-to-telescope rotation at the zenith angle the model was compiled for.
A model without an axes offset is traced as before.

`normalise_transmission=True` (default `False`) divides the response bound and
the roulette weight by the on-axis value `t0` of the telescope transmission (its
`on_axis`). If the model has no telescope transmission, or `on_axis` is not
positive, the divisor is 1. The kept-photon probability becomes
`fraction * acceptance * bound / t0`. The roulette ratio `weight / bound` does
not change, so the expected weight of a bunch photon is that of a full trace
divided by `t0`. sim_telarray already applies `t0` in its bunch loop, so the
telescope transmission then enters as `t(theta) / t0`. The output weights are
relative to `t0`; the caller applies `t0` itself, as sim_telarray does.

For each photon of a bunch (`ceil(size)` of them; the last one counts with the
fraction `size - floor(size)`) obdeect

1. looks up the external `acceptance` (atmosphere, quantum efficiency, ...; a
   table you supply, with wavelength in nm and emission height in metres; none
   means 1) and the upper bound of the model's responses at that wavelength
   (`response_bound`: reflectivities, degradation maps, telescope
   transmission);
2. keeps the photon with probability `fraction * acceptance * bound`, so that
   most lost photons cost no ray tracing;
3. traces the kept photon with the responses, and with `roulette=True`
   (default) keeps it again with probability `weight / bound`. Survivors have
   weight 1 and their expected number equals the sum of weights of a full
   trace. With `roulette=False` all of them stay and carry `min(1, weight / bound)`
   (the ratio is capped at 1).

The result is a dict of columns in input order, bunch by bunch and photon by
photon: `position_m`, `direction`, `time_ns`, `weight`, `wavelength_nm`,
`bunch_index` (row of the bunch in the input) and `sub_photon` (number within
the bunch). It does not depend on `threads`, on `pool`, or on how the bunches
are split over calls. Invalid input (a frame that is not orthonormal, an
acceptance outside [0, 1], a bunch of wavelength 0 without a density or law, a
bad size or direction, both `wavelength_density` and `wavelength_law`, an
invalid law) raises `ValueError` before any tracing.

The input plane `z = plane_z` (input units) must lie above every obscurer on the
incoming leg, such as the camera housing. Photons that start below such an
obscurer never meet it. Before tracing, `trace_bunches` computes, for each incoming
obscurer, the lowest height of the plane along the telescope axis over that
obscurer's footprint, and raises `ValueError` if it is below the obscurer's top.
For an untilted plane this is the plane height at the input origin. The tilt of
the plane is taken from `frame_rotation`, so the whole plane is bounded, not only
its origin. The message names the obscurer. Raise `plane_z` to fix this.

The random numbers of a photon depend on `(seed, bunch_id, sub_photon)`.
**Bunch ids must be unique per seed.** The recommended id is
`(event_number << 32) | bunch_index_in_event`, and a different `seed` per
telescope, derived by the caller from the run seed and the telescope number
(otherwise two telescopes would draw the same numbers for equal ids).
`classify_losses=False` is as for `trace`.
