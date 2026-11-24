# Optical models

## Frames and units

Lengths are metres; angles are degrees unless named `_rad`. In the telescope frame,
+z points along the optical axis towards the sky, and the dish vertex is at
z = -`mirror_offset`. The x and y axes follow the mirror-list coordinates; when the
telescope points above the horizon, x points down.

`star_photons` samples a disk centred on the dish vertex, perpendicular to the optical
axis, with radius 1.2 times the dish radius. `sampling_area_m2` is its area. The distance
of a point source is measured from the disk centre. Star photons start on a plane `start_z_m`,
1 m above the highest optical element: the larger of the vertex and facet heights used before,
and the bound of the optics (asphere rims and obscurer walls included, as `trace_bunches` checks
its input planes). Light is unchanged by the start height; the photon times are referenced to
the start plane.

## File format

An optical model (`obdeect.optical-model.v2`) is one JSON document:

- `format`, `sha256` (of the document without it), `provenance` (parameter versions,
  file hashes, seed, zenith angle, structure flag; its keys `configuration_seed`, `zenith_deg`,
  `structure`, `parameters` and `files_sha256` belong to the compiler), `notes`;
- `sampling`: the star-light disk (`disk_radius_m` > 0, `disk_z_m`, `start_z_m`);
- `optics`: `primary`, optional `secondary` (each `facets` or an `asphere` with
  `segments`, a `reflectivity` table, `scatter`, optional `degradation` map; facets may store the
  dish vertex `vertex_z_m` = -`mirror_offset` in metres, written by the compiler for segmented
  primaries and used for the in-telescope path, see [](../propagation/propagation.md)),
  `focal_surface` (plane or asphere), `obscurers` (cylinder, polyhedron, plate, frustum;
  acting on the legs `incoming`, `after_primary`, `after_secondary`),
  `telescope_transmission`, `propagation_group_index`, and optional `axes_offset_m`.

The C++ loader accepts either a structured JSON-compatible Python object or JSON text. It accepts
only known keys, checks every value and reports the path of the first error. Surface and obscurer
ids are unique. The `sampling`, `provenance` and `notes` blocks are validated too.

Text is checked against `sha256`: `load_model(text)` compares the entry with the canonical hash
(see `obdeect/ray_tracing/model_hash.hpp`), and so does the Python loader (`OpticalModel`, through
`obdeect.ray_tracing.hashing`). A structured object has no canonical text in C++, so the Python
caller checks it. The text parser and the canonical hash accept the same documents: the four JSON
white-space characters, `\uXXXX` escapes (surrogate pairs decoded, lone surrogates rejected) and a
nesting depth of 64.

## Compiling from simulation-models

`obdeect.ray_tracing.importers.simulation_models` is the only code that knows simulation-models
parameter names and sim_telarray conventions; the core model does not depend on it.

- The air group index comes from the telescope's site (its atmosphere at the observation level).
  A telescope without a site (for example `MSTx-FlashCam`) needs `propagation_group_index`
  (`obdeect ray-tracing compile --propagation-group-index`); `compile_model` has no default.
- Values with no optical effect are accepted exactly: a zero shadow offset selects the automatic
  position of the secondary shadow, a zero outer power of the telescope transmission means 1, a
  zero mirror-list height or focal length is the automatic prescription, and the panel id (the
  seventh mirror-list column) labels a mirror without acting on it. Every other value the tracer
  cannot use is rejected: a method-0 transmission amplitude, a dual-mirror outer power other than
  0 or 1, a negative shadow offset, a single zenith term of the focus offset, an overshooting
  reflectivity, a zero dish length, and an RPOL table with a 2-D scheme other than 1, a zero scale
  or a z limit.
- Parameters without a no-op default (`mirror_class`, `camera_body_shape`, `dish_shape_length`)
  are required.

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

The tracer follows every photon through one loop: after each reflection it
searches the surfaces again until the photon is detected, absorbed or escapes. Photons from
above (downward, the usual case) meet the primary on its reflecting face first; photons from
below (upward, `cz > 0` in the input frame) meet the back of the primary and, if the light passes
the primary's hole, the reflecting face of the secondary. The side of a hit follows from the
surface normal the tracer returns:

| Surface | Normal returned by the tracer | Reflecting side | Hit from the other side |
| --- | --- | --- | --- |
| Primary facets | towards the reflecting side | the normal side: `n . d < 0` | blocked (back face) |
| Primary asphere | `+z` | above (`n . d < 0`) | blocked (back face) |
| Secondary asphere | `+z` (a graph surface) | below, facing the primary (`n . d > 0`) | blocked (back face) |
| Focal surface | `+z` | the receiving side (below for a single mirror, above for a dual mirror; see below) | passed (see below) |

The secondary is an asphere whose sag is reversed: its reflecting face looks down, against the
normal. The reflected direction does not depend on the side.

### Focal surface

The focal surface detects light on its receiving side: the side that the light of the sequence
arrives from. For a single mirror the light leaves the primary upwards, so the receiving side is
the lower one; for a dual mirror the light leaves the secondary downwards, and the focal surface
lies below the secondary, so the receiving side is the upper one.

- On the leg after the last mirror (the final leg), a hit from either side is detected: the
  sequential path, unchanged.
- On any other leg (before the last mirror, including light from below before any reflection),
  an asphere focal surface is detected inside its extent (`outer_radius_m`) from its receiving
  side. A plane focal surface has no extent: a hit from its receiving side is `direct_focal`, a
  counted loss. A hit from the other side is passed: it is not a surface the light meets, and the
  camera body shadows it (obscurers on the legs; for a dual mirror the camera back plate on the
  after-primary leg).

The model has no fiducial radius. A plane focal surface is unbounded; the camera body obscurers
(`camera_body_diameter`, `camera_depth`) give it a lateral extent only as shadows, not as a
detection area. The asphere's outer radius is the only lateral extent of a focal surface.

### Loop rules

- Obscurers: those of the current leg, the first one nearer than the surface hit blocks the photon
  (ties go to the surfaces: primary, secondary, focal surface, then obscurers). Light from below
  meets the obscurers of every leg before its first reflection, as it crosses them all on its way
  up. The secondary is tested on the incoming leg only where it blocks the sky
  (`secondary_blocks_incoming`) or the photon comes from below.
- Legs: `incoming` until the first reflection, then `after_primary` or `after_secondary` after the
  last reflection on that mirror.
- A photon with no reflection yet that meets nothing is `missed_primary`; after a reflection on the
  primary of a dual mirror, `missed_secondary`; otherwise `missed_focal_surface`.
- Reflections are limited to 16 (`kMaxReflections`); a photon that would reflect again ends as
  `bounce_limit` on that mirror.
- A reflection keeps the mirror's reflectivity (at the incidence angle), degradation and scatter.
  The first reflection on each mirror draws its random numbers as before (`primary_scatter`,
  `secondary_scatter`); later reflections use `reflection_stream(k)`, k their reflection number.
  The output `primary_*` and `secondary_*` columns record the first hit on each mirror;
  `reflections` counts all reflections.

Status values (`STATUS_NAMES`, in this order): `detected`, `missed_primary`, `missed_secondary`,
`missed_focal_surface`, `blocked`, `invalid_input`, `numerical_failure`, `non_sequential`,
`direct_focal`, `bounce_limit`. `non_sequential` is no longer produced: the reflecting face of the
secondary is traced on (its value is kept so that the numbering is stable). `direct_focal` is
produced only for a plane focal surface (no extent) hit from its receiving side before the last
mirror. Every status but `detected` has weight zero.
