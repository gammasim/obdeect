# Scope and limitations

The trace ends at the focal surface. A camera window can optionally be included in the
optical model; effects not represented there remain outside the trace. Their impact
depends on the instrument and its model parameters.

## Misalignment of aspheric (dual-mirror) mirrors

The primary and secondary of the Schwarzschild-Couder telescopes (SST, SCT) are even
polynomial surfaces (`mirror_class` 2) built from segments. In the real telescope each
segment, and the secondary as a whole, can be tilted or shifted from the design.

The alignment errors (`mirror_align_random_distance`, `_horizontal`, `_vertical`) and
`random_focal_length` act only on facets of segmented dishes (`mirror_class` 0). Dual-mirror
surfaces are traced at their design position; `primary_segmentation` and
`secondary_segmentation` define gaps between segments. `mirror_reflection_random_angle` tilts
the surface normal at each reflection (a convention, not a ray angle: see
[mirror scatter conventions](optics.md#mirror-scatter-conventions)). For the resulting PSF
comparisons, see [Validation](VALIDATION.md).

## Curvature and refraction of the camera window (LST and others)

A camera window is a transparent plate, flat or curved, in front of the pixels. The optical
model can trace refraction at both surfaces, Fresnel loss, bulk transmittance and added flight
time. The rim and flange absorb photons. See [Optical elements](optical-elements.md).

**Size of the effects** for a flat plate of thickness t and refractive index n (PMMA n = 1.49
as an example; take each camera window's thickness from the camera documentation):

| t | lateral shift at 10 / 16 / 25 deg incidence | focus shift t(1 - 1/n) | extra delay (n_g - 1) t / c |
| --- | --- | --- | --- |
| 3 mm | 0.17 / 0.28 / 0.46 mm | 1.0 mm | 5 ps |
| 5 mm | 0.29 / 0.47 / 0.77 mm | 1.6 mm | 8 ps |
| 10 mm | 0.58 / 0.95 / 1.5 mm | 3.3 mm | 17 ps |

The lateral shift is t sin θ (1 - cos θ / sqrt(n² - sin² θ)); the delay uses a group index
n_g ≈ 1.5. The mean incidence angle on the LST focal plane is 16 deg
([validation](VALIDATION.md)). For comparison: the LST D80 is about 3 cm, its pixels about
5 cm wide, and a focus error δ blurs the image by about δ D / f = δ / 1.2 for the LST. A shell
of uniform thickness t and radius of curvature R has optical power of order (n - 1) t / R²: a
focal length of kilometres for millimetres and metres, negligible against 28 m.

The focus shift is the only effect in the table that is not negligible. The interplay with
`focus_offset` is being studied: it was tuned to the measured PSF without a traced window, and
a curved shell changes focus differently across the field. The 5 mm row gives representative
LST-scale values; the actual result depends on the window shape and ray incidence.

## Other effects not traced

- **Telescope structure of LST, SST and SCT.** Only the MST has a `telescope_structure` table
  in the models examined. The other telescopes use `telescope_transmission`, a parametrisation
  of shadowing.
- **Pixels, light guides, photodetectors:** downstream camera simulation.
- **Refraction in air** between an illuminator and the telescope: a few arcsec over hundreds
  of metres.
- **Diffraction and polarisation:** irrelevant at the size of IACT mirrors and for unpolarised
  Cherenkov light.

## Non-sequential tracing

Every photon is traced until it is detected, absorbed, or escapes, after up to 16 reflections
(`bounce_limit` beyond). Remaining limits:

- A focal surface is detected only inside its extent: an asphere's outer radius. A plane focal
  surface has no extent, so light from its receiving side before the last mirror is counted as
  `direct_focal` (a loss) rather than detected. The camera body obscurers give the lateral shadow.
- Light from below meets the obscurers of every leg before its first reflection; after a
  reflection the obscurers of the current leg apply. Structure on the after-primary leg is therefore
  not met by a photon that has reflected off the secondary and comes back down.
- The bunch roulette (`trace_bunches`) uses `response_bound` as its scale. A photon is kept with
  probability min(1, weight / bound). A photon whose weight exceeds the bound (a direct detection
  without reflection, or a dual-mirror detection after the primary alone) is kept with probability
  one and output with weight weight / bound, so the expected output weight stays the traced weight;
  without roulette the output weight is weight / bound (it can exceed one for such photons). The
  estimator is unbiased (native `bunch_roulette_unbiased`).
- Multi-bounce photons are traced with the same reflectivity, degradation and scatter as the
  first reflection, at each reflection's own incidence angle.
- The first hit of each mirror is recorded in the output (`primary_*`, `secondary_*`); the others
  are counted in `reflections` only.

## Sphere runs

- The enclosure check bounds each mirror and obscurer by the box of its support functions, so a
  sphere within a few centimetres of the optics can be rejected although it encloses them. The
  check never accepts a sphere that does not enclose the optics.
- The focal surface is a plane of unbounded extent in the model. Its only check is that the plane
  reaches the sphere (its height above the sphere centre is at most R). An asphere focal surface
  is bounded like the mirrors.
- A row emitted inside the sphere is traced from its emission point; the tracer does not check that
  the point lies outside the optics. Emission inside a mirror or an obscurer is not detected as
  such.
- The output `path` adds the traced path inside the telescope in metres (for a plane run, the move
  to the input plane is included, which can make the path of a photon emitted below the input
  plane negative: such input is rejected). The path of a reflected photon is the sum of its legs.
- Each telescope's radius comes from `telescope_r_m` indexed by telescope id; a run without it, or
  with a shorter array, is rejected.
