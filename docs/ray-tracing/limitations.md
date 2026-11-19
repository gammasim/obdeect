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
