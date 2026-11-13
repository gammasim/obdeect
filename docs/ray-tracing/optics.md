# Optics

obdeect ray tracing traces photons from an entrance surface to the camera entrance.
The optical model describes mirrors, obscurers, the focal surface and optical
responses.

```{toctree}
:maxdepth: 1

optical-elements
telescopes
shadowing
limitations
```

## Tracing a photon

1. The source supplies position, direction, wavelength, time, weight and id.
   Star light (`star_photons`) is sampled on a disk (radius 1.2
   times the dish radius, perpendicular to the optical axis); a source at finite
   distance emits from its position.
2. Incoming leg: nearest of primary mirror, secondary back side (if it blocks
   incoming light) and obscurers. Lost photons record where and on what.
3. Primary: the nearest valid mirror hit produces a specular reflection.
   Scattering by `mirror_reflection_random_angle` (two Gaussian components): for
   facets the reflected direction is rotated, for polynomial mirrors the normal
   is tilted by half the angle.
4. Secondary (dual mirror): same, with the after-primary obscurers.
5. Focal surface: plane or polynomial; the photon's position, direction, time
   and incidence angles are recorded. obdeect ray tracing stops at the camera entrance;
   pixels, light guides and the camera window are not part of the optical model.

## Mirror scatter conventions

`mirror_reflection_random_angle` is applied with the convention of the
simulation it comes from, not as a measured RMS of the reflected-ray angle. Per
reflection one Gaussian width is drawn (the second component with the given
fraction), then two independent deviations a1 and a2 with that width:

- **Segmented mirrors** (`method` `outgoing_angles`): the specularly reflected
  ray is rotated by a1 about one tangent axis of the facet and by a2 about the
  other. A rotation about an axis in the tangent plane moves the ray by
  sigma times the sine of the angle between ray and axis, so the mean squared
  deviation is sigma^2 (1 + cos^2 i) for the incidence angle i: sigma per axis
  at normal incidence, and the in-plane component shrinks with the cosine of i.
  The axes are orthonormalised against the normal at the hit point (on a curved
  facet the vertex tangent is not perpendicular to it).
- **Polynomial mirrors** (`method` `surface_slopes`): the slopes dz/dx and
  dz/dy of the surface at the hit point are each changed by half of a1 and a2,
  so the reflected ray deviates by about the full angle where the surface is
  flat. A slope is not an angle: a slope change ds tilts the normal by
  ds cos^2(theta) in the radial and ds cos(theta) in the tangential direction,
  with theta the tilt of the surface. On strongly tilted parts of a secondary
  the effective width is therefore smaller than the nominal one, and
  anisotropic. Compare PSFs, not parameter values, between this convention and
  one defined as a ray angle.

Both are spatially uncorrelated between photons, equal over the whole mirror
and independent of wavelength. They describe the core of the PSF; a wide
scattering halo is not modelled.

Weights multiply the responses met (unless geometry only). Each random number
comes from the Philox-4x32-10 generator shared by all components, keyed by the seed,
event and telescope (photon streams; 0 otherwise), photon id (below 2^40) and stream, so
results do not depend on threads, blocks or photon order. Times use the group
index stored in the model.
