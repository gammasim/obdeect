# Optical models

## Frames and units

Lengths are metres; angles are degrees unless named `_rad`. In the telescope frame,
+z points along the optical axis towards the sky, and the dish vertex is at
z = -`mirror_offset`. The x and y axes follow the mirror-list coordinates; when the
telescope points above the horizon, x points down.

`star_photons` samples a disk centred on the dish vertex, perpendicular to the optical
axis, with radius 1.2 times the dish radius. `sampling_area_m2` is its area. The distance
of a point source is measured from the disk centre.

## File format

An optical model (`obdeect.optical-model.v2`) is one JSON document:

- `format`, `sha256` (of the document without it), `provenance` (parameter versions,
  file hashes, seed, zenith angle, structure flag), `notes`;
- `sampling`: the star-light disk;
- `optics`: `primary`, optional `secondary` (each `facets` or an `asphere` with
  `segments`, a `reflectivity` table, `scatter`, optional `degradation` map; facets may store the
  dish vertex `vertex_z_m` = -`mirror_offset` in metres, written by the compiler for segmented
  primaries and used for the in-telescope path, see [](../propagation/propagation.md)),
  `focal_surface` (plane or asphere), `obscurers` (cylinder, polyhedron, plate, frustum;
  acting on the legs `incoming`, `after_primary`, `after_secondary`),
  `telescope_transmission`, `propagation_group_index`, and optional `axes_offset_m`.

The C++ loader accepts either a structured JSON-compatible Python object or JSON text. It accepts
only known keys, checks every value and reports the path of the first error. Surface and obscurer
ids are unique.

## Axes offset

Some telescopes have the optical origin displaced from the rotation axes. The sim_telarray
`axes_offsets` record gives the horizontal offset between the azimuth and altitude axes,
then the altitude-to-optical-axis offset, in its units (cm in sim_telarray).

The compiler stores the offset in the telescope frame, in metres, as
`optics.axes_offset_m` = [x, y, z]. The key is written only when the offset is non-zero,
and needs the compilation zenith angle. For components `(O, T)` and zenith `theta`, the
offset is `(O cos(theta) - T, 0, O sin(theta))`. The loader requires three numbers.

Only `trace_bunches` applies the offset (see [Air-shower photons](air-shower-photons.md#input-plane)).
`trace`, `star_photons` and the emitter paths do not use it.

For mirrors, obscurers and the telescope-structure table, see
[Optical elements](optical-elements.md).

## Reflecting sides and photons from below

Photons from above (downward, the usual case) meet the primary on its reflecting face. Photons
from below (upward, `cz > 0` in the input frame) meet the back of the primary, and the other
faces of the optics from their own sides. Phase 1 of the ray tracer decides the side of each
first hit from the surface normal (`trace_photon`, see [Light sources](light-sources.md#photons-from-below)):

| Surface | Normal returned by the tracer | Reflecting side | Photon from below |
| --- | --- | --- | --- |
| Primary facets | towards the reflecting side | the normal side: `n . d < 0` | back face (blocked) |
| Primary asphere | `+z` | above (`n . d < 0`) | back face (blocked) |
| Secondary asphere | `+z` (a graph surface) | below, facing the primary (`n . d > 0`) | reflecting face (non-sequential) |
| Focal surface | `+z` | receives light from below (`n . d > 0`) | receiving side (direct focal) |

The secondary is an asphere whose sag is reversed: its reflecting face looks down, against
the normal. A facet is hit on its reflecting side only if its normal points to it, so the back
of a facet is not a mirror hit for downward photons; the two-sided search is used for photons
from below only.

Status values (`STATUS_NAMES`, in this order): `detected`, `missed_primary`, `missed_secondary`,
`missed_focal_surface`, `blocked`, `invalid_input`, `numerical_failure`, and, for photons from
below, `non_sequential` (the reflecting face of the secondary, out of sequence) and `direct_focal`
(the focal surface, on its receiving side, before any mirror). Every status but `detected` has
weight zero. Full non-sequential tracing (bounces off all surfaces until detection or escape) is
planned as phase 2; until then these statuses are counted losses.
