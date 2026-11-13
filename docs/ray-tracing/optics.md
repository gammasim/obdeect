# Optics

Photons travel from the source to the camera entrance. The optical model holds the
mirrors, obscurers, the focal surface and the optical responses.

```{toctree}
:maxdepth: 1

optical-elements
telescopes
shadowing
limitations
```

## Tracing a photon

1. The source supplies position, direction, wavelength, time, weight and id.
   Star light (`star_photons`) is sampled on a disk (radius 1.2 times the dish radius,
   perpendicular to the optical axis); a source at finite distance emits from its position.
2. Incoming leg: nearest of primary mirror, secondary back side (if it blocks incoming
   light) and obscurers. Lost photons record where and on what.
3. Primary: the nearest valid mirror hit gives a specular reflection, scattered by
   `mirror_reflection_random_angle` (two Gaussian components). For facets the reflected
   direction is rotated; for polynomial mirrors the normal is tilted by half the angle.
4. Secondary (dual mirror): same, with the after-primary obscurers.
5. Focal surface: plane or polynomial. Position, direction, time and incidence angles are
   recorded. Tracing stops at the camera entrance; pixels, light guides and the camera
   window are not part of the optical model.

## Mirror scatter conventions

`mirror_reflection_random_angle` follows the convention of the simulation it comes from.
It is not a measured RMS of the reflected-ray angle. Per reflection, one Gaussian width is
drawn (the second component with the given fraction), then two independent deviations
a1 and a2 with that width:

- **Segmented mirrors** (`method` `outgoing_angles`): the reflected ray is rotated by a1
  about one tangent axis of the facet and by a2 about the other. A rotation about a tangent
  axis moves the ray by sigma times the sine of the angle between ray and axis, so the mean
  squared deviation is sigma^2 (1 + cos^2 i) for incidence angle i: sigma per axis at normal
  incidence, and the in-plane component shrinks with the cosine of i. The axes are
  orthonormalised against the normal at the hit point (on a curved facet the vertex tangent
  is not perpendicular to it).
- **Polynomial mirrors** (`method` `surface_slopes`): the slopes dz/dx and dz/dy at the hit
  point each change by half of a1 and a2, so the reflected ray deviates by about the full
  angle where the surface is flat. A slope is not an angle: a slope change ds tilts the
  normal by ds cos^2(theta) radially and ds cos(theta) tangentially, with theta the tilt of
  the surface. On strongly tilted parts of a secondary the effective width is smaller than
  the nominal one, and anisotropic. Compare PSFs, not parameter values, between this
  convention and one defined as a ray angle.

Both are uncorrelated between photons, equal over the whole mirror and independent of
wavelength. They describe the PSF core; a wide scattering halo is not modelled.

Weights multiply the responses met (unless geometry only). Each random number comes from
the Philox-4x32-10 generator, keyed by the seed, event and telescope (photon streams; 0
otherwise), photon id (below 2^40) and stream. Results therefore do not depend on threads,
blocks or photon order. Times use the group index stored in the model.
