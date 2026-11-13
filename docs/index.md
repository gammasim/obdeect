# obdeect

```{image} _static/images/obdeect_21ratio_small.jpg
:alt: obdeect, a simulation chain for imaging atmospheric Cherenkov telescopes
:align: center
:class: homepage-logo
```

A simulation chain for imaging atmospheric Cherenkov telescopes. The modules exchange
photons as photon streams: light sources or CORSIKA air showers, the atmosphere, the
telescope optics, and the camera simulation in sim_telarray. Each module is a component
with a C++20 core, a Python API and commands of the form `obdeect <component> <command>`.

:::{container} homepage-actions

```{button-link} https://github.com/gammasim/obdeect
:color: primary
:shadow:

{fab}`github` GitHub Repository
```

```{button-link} https://github.com/gammasim/obdeect/blob/main/LICENSE
:color: secondary
:outline:

{fas}`scale-balanced` License BSD-3
```

:::

---

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`stream;pst-color-primary` Photon stream
:link: stream/index
:link-type: doc
:class-card: sd-shadow-sm

The format between the modules, and conversion from CORSIKA and to sim_telarray.
:::

:::{grid-item-card} {fas}`cloud;pst-color-primary` Atmospheric propagation
:link: propagation/index
:link-type: doc
:class-card: sd-shadow-sm

Extinction, refraction, timing, thinning and scattering in the atmosphere.
:::

:::{grid-item-card} {fas}`lightbulb;pst-color-primary` Light sources
:link: sources/index
:link-type: doc
:class-card: sd-shadow-sm

Flashers, LEDs, lasers, stars, night-sky background and muons.
:::

:::{grid-item-card} {fas}`eye;pst-color-primary` Ray tracing
:link: ray-tracing/index
:link-type: doc
:class-card: sd-shadow-sm

Mirrors, shadowing and optical responses, from the telescope to the camera.
:::

::::

```{toctree}
:hidden:
:maxdepth: 2

installation
stream/index
propagation/index
sources/index
ray-tracing/index
changelog
```
