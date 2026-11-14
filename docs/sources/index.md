# Light sources

Besides air showers, Cherenkov telescopes are calibrated and tested with light that people
control: flashers on the dish, LEDs in the camera, lasers, stars, point sources and the night-sky
background. This component simulates such light and writes it as a photon stream, ready for the
atmosphere, the optics and sim_telarray. It replaces the sim_telarray light-emission package.

A source is not a hard-wired device. It is described by a position, a direction and an angular
profile, a spectrum and a pulse shape, so one model covers all programs. The photons start at the
sky, telescope-sphere, telescope or camera level; to avoid simulating photons that never reach the
telescope, a source can generate only photons aimed at a target.

Start with the [source model](sources.md), or pick a program:

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {fas}`lightbulb;pst-color-primary` Flasher
:link: programs/flasher
:link-type: doc
:class-card: sd-shadow-sm

Flashers on the dish axis or at the camera corners.
:::

:::{grid-item-card} {fas}`microchip;pst-color-primary` LED
:link: programs/led
:link-type: doc
:class-card: sd-shadow-sm

LEDs in front of camera pixels, to test triggers.
:::

:::{grid-item-card} {fas}`circle-dot;pst-color-primary` Point source
:link: programs/point-source
:link-type: doc
:class-card: sd-shadow-sm

A point light source above the telescope.
:::

:::{grid-item-card} {fas}`arrows-left-right;pst-color-primary` Beam scan
:link: programs/beam-scan
:link-type: doc
:class-card: sd-shadow-sm

A parallel beam scanned across the camera area.
:::

:::{grid-item-card} {fas}`bolt;pst-color-primary` Laser beam
:link: programs/laser-beam
:link-type: doc
:class-card: sd-shadow-sm

A laser crossing the field of view.
:::

:::{grid-item-card} {fas}`moon;pst-color-primary` Night sky
:link: programs/night-sky
:link-type: doc
:class-card: sd-shadow-sm

Night-sky background on a telescope sphere.
:::

:::{grid-item-card} {fas}`wind;pst-color-primary` Muon track
:link: programs/muon-track
:link-type: doc
:class-card: sd-shadow-sm

Cherenkov light of one muon.
:::

:::{grid-item-card} {fas}`diagram-project;pst-color-primary` Use with sim_telarray
:link: ../user-guide/pipelines
:link-type: doc
:class-card: sd-shadow-sm

Pipelines from each source to sim_telarray.
:::

:::{grid-item-card} {fas}`book;pst-color-primary` Specification
:link: spec/index
:link-type: doc
:class-card: sd-shadow-sm

Exact definitions, for implementers.
:::

::::

```{toctree}
:hidden:
:maxdepth: 1

sources
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Programs

programs/flasher
programs/led
programs/point-source
programs/beam-scan
programs/laser-beam
programs/night-sky
programs/muon-track
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Reference

spec/index
```
