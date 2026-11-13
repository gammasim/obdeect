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
  `segments`, a `reflectivity` table, `scatter`, optional `degradation` map),
  `focal_surface` (plane or asphere), `obscurers` (cylinder, polyhedron, plate, frustum;
  acting on the legs `incoming`, `after_primary`, `after_secondary`),
  `telescope_transmission`, `propagation_group_index`, and optional `axes_offset_m`.

The C++ loader accepts only known keys, checks every value and names the JSON path of
the first error. Surface and obscurer ids are unique.

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
