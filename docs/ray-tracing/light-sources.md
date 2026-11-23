# Light sources

Choose a source by the light path you study. The package samples stars and point
sources, calibration emitters, and external air-shower photons.

| Source | Typical use | Interface |
| --- | --- | --- |
| Star at infinity | Point-spread function, effective area, off-axis response | `star_photons` |
| Point source at finite distance | Focus and distance-dependent imaging | `star_photons(distance_m=...)` |
| Flasher | Direct illumination of the camera | `emitter_photons`, then `trace(direct=True)` |
| Illuminator | External calibration light passing through the mirrors | `emitter_photons`, then `trace` |
| Air-shower photon bunches | Cherenkov light from an external shower simulation | `trace_bunches` |

## Stars and finite-distance point sources

`star_photons` samples the entrance disk. At infinity, rays are parallel; at a finite
distance, their directions follow the point source. `towards_source` points to the source
in the telescope frame. `distance_m` is measured from the dish vertex.

```python
# On-axis star at infinity (the default).
photons = model.star_photons(100_000, seed=2, wavelength_nm=400.0)

# On-axis point source 10 km from the dish vertex.
photons = model.star_photons(100_000, seed=2, distance_m=10_000.0)
arrivals = model.trace(photons, seed=2)
```

Use `towards_source` for an off-axis source. On the command line, use
`obdeect ray-tracing trace --field-deg X Y` and optionally `--distance-m`.
This sampler fixes the number of incident photons; it does not model an emitter's
luminosity or pulse shape.

## Flashers and illuminators

Both use a point emitter with a position, wavelength, angular intensity profile and
optional pulse shape. `emitted` is the mean number of photons emitted into all directions.
Sampling targets a disk and uses photon weights for the illumination; it does not generate
every emitted photon.

A **flasher** illuminates the camera directly from inside the telescope. Set the target
disk at the camera and trace with `direct=True`: obscurers are tested, but mirrors and
telescope transmission are bypassed.

An **illuminator** shines onto the telescope from outside. Target the entrance disk (the
default) and use the normal mirror path. For bright pulses, `trace_emitter` or
`trace_emitter_blocks` generates and traces photons in blocks.

Example: a flasher with an 11-degree cone and a 2.5 ns flat pulse:

```python
camera = model.focal_surface_z_m
photons = model.emitter_photons(
    1e6,  # photons emitted into 4 pi
    seed=3,
    position_m=(0.2, 0.2, camera - 16.5),  # telescope frame
    intensity=([0.0, 11.0, 11.0001], [1.0, 1.0, 0.0]),  # 11 deg filled cone
    pulse=([-1.25, 1.25], [1.0, 1.0]),  # 2.5 ns flat pulse
    wavelength_nm=385.0,
    disk=((0.0, 0.0, camera), (0.0, 0.0, 1.0), 1.2),  # the camera
)
arrivals = model.trace(photons, seed=3, direct=True)
```

The candidate count is Poisson distributed and weights lie in [0, 1]. Keeping each
candidate with the probability of its weight gives the photons crossing the target disk;
the sum of weights estimates their mean number.

## Photons from below

Photons may move up as well as down. `towards_source` may have any nonzero z component: a
source below the telescope (`towards_source` with negative z) gives upward photons, which
start on the input plane below the optics.

`star_photons` starts such photons on a plane a margin below the lowest optic. Bunches from
`trace_bunches` with cz > 0 start on the lower input plane, which lies a margin of 1 m below
every obscurer (any leg) and every mirror (primary and secondary) along the telescope axis.

The trace branch is chosen by the direction in the telescope frame: a photon whose traced
direction has z >= 0 takes the first-hit path below, whatever the sign of cz in the input frame.
For a telescope at a large zenith angle, a ground-frame photon with cz < 0 can move up in the
telescope frame; it arrives from behind the dish and is blocked by its back. Photons moving
down in the telescope frame are traced exactly as before.

`trace` and `trace_bunches` report the first hit of an upward photon:

- an obscurer, or the back face of a mirror: `blocked`, with the surface recorded;
- the reflecting face of the primary: the sequential path continues (the photon may be detected);
- the reflecting face of the secondary: `non_sequential` (phase 1 does not trace further);
- the focal surface from below (on its receiving side): `direct_focal`;
- nothing: `missed_primary`.

Phase 1 counts `non_sequential` and `direct_focal` photons as losses with weight zero; they
are the light that the planned phase 2 (non-sequential tracing) will recover. A photon that
moves horizontally (zero component along the axis of its input plane) has no plane to start on
and is rejected with an error.

## Sphere rows

In a sphere run (see [Photon stream](photon-stream.md#input-frames)) each ground-frame row is a
start point inside or on the telescope sphere, and it is traced from that point: its position is
transformed into the telescope frame as a point (rotation, then the axes offset), with its own z.
Its time is the row's time at that point. Nothing is moved to an input plane, so the plane checks
do not apply to these rows.

The branch follows the direction in the telescope frame, as for bunches from below: a row with
traced z < 0 takes the sequential path (the primary's reflecting face, then the secondary and the
focal surface), and a row with traced z >= 0 takes the first-hit path. A horizontal ground row
(cz = 0) on a tilted telescope has traced z = R d_z and so takes whichever branch its traced
direction selects; with traced z = 0 it takes the first-hit path.

The sequential path looks ahead from the start point only. A photon that starts between the
primary and the focal surface and moves towards the primary is reflected by it; one that moves
towards the camera reaches the focal surface from below (`direct_focal`). Obscurers behind the
start point are not met, which is correct for a start point inside the telescope.

For the same photons, the sphere row of each downward photon is the entry point of its line into
the sphere, and the tracer gives the same detected photons as the plane run (checked in
[Validation](VALIDATION.md#sphere-runs-ray-tracing-phase-5d)).

## Air-shower photons and custom inputs

`trace_bunches` traces photon bunches from an external air-shower simulation; it does not
generate the shower. See [Shower photon bunches](batch-tracing.md#shower-photon-bunches).

For other sources, pass photon arrays to `trace`. The required fields and output arrays
are in [Python usage](python-guide.md).
