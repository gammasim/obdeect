# Atmospheric propagation

Between its emission and the telescope, a photon crosses up to several tens of kilometres of air.
This component moves photon streams along that path. It computes where the photon arrives and
when, how much light the atmosphere absorbs, and optionally how much light scatters towards the
telescope. The result is a stream of photons in the ground frame that the optics step can
use.

All effects are applied in a fixed order and recorded in the
[efficiency status](../stream/efficiency-status.md) of the stream, so a later module knows what is
done. To use it, start with the [atmosphere stage](propagation.md); the other pages describe the
physics.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`play;pst-color-primary` Atmosphere stage
:link: propagation
:link-type: doc
:class-card: sd-shadow-sm

The command `obdeect propagation run`: steps, options, output and Python interface.
:::

:::{grid-item-card} {fas}`cloud;pst-color-primary` Extinction
:link: extinction
:link-type: doc
:class-card: sd-shadow-sm

Transmission, wavelength assignment and thinning.
:::

:::{grid-item-card} {fas}`wave-square;pst-color-primary` Refraction and timing
:link: refraction
:link-type: doc
:class-card: sd-shadow-sm

Ray bending and air delay from an atmospheric profile.
:::

:::{grid-item-card} {fas}`burst;pst-color-primary` Scattering
:link: scattering
:link-type: doc
:class-card: sd-shadow-sm

Single scattering of light-source photons towards the telescope.
:::

:::{grid-item-card} {fas}`table;pst-color-primary` Legacy tables
:link: legacy-tables
:link-type: doc
:class-card: sd-shadow-sm

Convert older extinction tables to the current format.
:::

::::

```{toctree}
:hidden:
:maxdepth: 1

propagation
extinction
refraction
scattering
legacy-tables
```
