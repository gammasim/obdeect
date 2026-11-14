# obdeect

```{image} _static/images/obdeect_21ratio_small.jpg
:alt: obdeect, a simulation chain for imaging atmospheric Cherenkov telescopes
:align: center
:class: homepage-logo
```

obdeect simulates how light travels from its source to the camera of an imaging
atmospheric Cherenkov telescope. The light comes from CORSIKA air showers or from artificial
sources such as flashers, LEDs, lasers and stars. It crosses the atmosphere, is traced through
the telescope optics, and ends in the camera simulation of sim_telarray.

The steps are separate programs that pass **photon streams** to each other, so you can
combine them as in a Unix pipeline and replace any one of them.

:::{container} homepage-actions

```{button-link} user-guide/index
:color: primary
:shadow:

{fas}`bolt` Start with the user guide
```

```{button-link} https://github.com/gammasim/obdeect
:color: secondary
:outline:

{fab}`github` GitHub repository
```

```{button-link} https://github.com/gammasim/obdeect/blob/main/LICENSE
:color: secondary
:outline:

{fas}`scale-balanced` License BSD-3
```

:::

```{raw} html
<svg class="chain-figure" viewBox="0 0 760 150" role="img"
     aria-label="Light sources and CORSIKA feed the atmosphere, then the optics, then sim_telarray; photon streams connect all steps."
     font-family="inherit" font-size="14" fill="currentColor">
  <g fill="none" stroke="var(--pst-color-primary)" stroke-width="2">
    <rect x="4" y="12" width="150" height="40" rx="8"/>
    <rect x="4" y="82" width="150" height="40" rx="8"/>
    <rect x="224" y="47" width="140" height="56" rx="8"/>
    <rect x="434" y="47" width="140" height="56" rx="8"/>
    <rect x="644" y="47" width="112" height="56" rx="8" stroke-dasharray="5 4"/>
  </g>
  <g fill="none" stroke="currentColor" stroke-width="1.5" opacity="0.7">
    <path d="M154 32 C190 32 190 75 224 75"/>
    <path d="M154 102 C190 102 190 75 224 75"/>
    <path d="M364 75 H434"/>
    <path d="M574 75 H644"/>
  </g>
  <g text-anchor="middle">
    <text x="79" y="37">Air showers</text>
    <text x="79" y="107">Light sources</text>
    <text x="294" y="72" font-weight="600">Atmosphere</text>
    <text x="294" y="92" font-size="12" opacity="0.8">propagation</text>
    <text x="504" y="72" font-weight="600">Telescope</text>
    <text x="504" y="92" font-size="12" opacity="0.8">ray tracing</text>
    <text x="700" y="72">sim_telarray</text>
    <text x="700" y="92" font-size="12" opacity="0.8">camera</text>
    <text x="394" y="143" font-size="12" opacity="0.8">every arrow is a photon stream</text>
  </g>
</svg>
```

## Where to go

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`bolt;pst-color-primary` User guide
:link: user-guide/index
:link-type: doc
:class-card: sd-shadow-sm

Install obdeect, run a first pipeline and connect it to sim_telarray.
:::

:::{grid-item-card} {fas}`stream;pst-color-primary` Photon stream
:link: stream/index
:link-type: doc
:class-card: sd-shadow-sm

The format that connects the steps, with converters from CORSIKA and to sim_telarray and GrOptics.
:::

:::{grid-item-card} {fas}`cloud;pst-color-primary` Atmospheric propagation
:link: propagation/index
:link-type: doc
:class-card: sd-shadow-sm

Extinction, refraction, timing and scattering between the emission point and the telescope.
:::

:::{grid-item-card} {fas}`lightbulb;pst-color-primary` Light sources
:link: sources/index
:link-type: doc
:class-card: sd-shadow-sm

Flashers, LEDs, lasers, point sources, night-sky background and muons.
:::

:::{grid-item-card} {fas}`eye;pst-color-primary` Ray tracing
:link: ray-tracing/index
:link-type: doc
:class-card: sd-shadow-sm

Mirrors, shadowing and optical responses, from the telescope entrance to the camera.
:::

:::{grid-item-card} {fas}`clock-rotate-left;pst-color-primary` Changelog
:link: changelog
:link-type: doc
:class-card: sd-shadow-sm

What changed in each release.
:::

::::

```{toctree}
:hidden:
:maxdepth: 2
:caption: User guide

user-guide/index
user-guide/installation
user-guide/quickstart
user-guide/pipelines
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: Components

stream/index
propagation/index
sources/index
ray-tracing/index
```

```{toctree}
:hidden:
:caption: Project

changelog
```
