# Light sources

Besides air showers, Cherenkov telescopes are calibrated and tested with light that people
control: flashers on the dish, LEDs in the camera, lasers, stars, point sources and the night-sky
background. This component simulates such light and writes it as a photon stream (`.pst`), ready
for the atmosphere, the optics and sim_telarray. It replaces the sim_telarray light-emission package.

A source is not a hard-wired device. It is a geometry (point, line, plane wave or sky background)
with angular, spectral and pulse profiles. A run file lists the sources of each event. To avoid
simulating photons that never reach the telescope, a source can generate only photons aimed at a
target.

```sh
obdeect sources run flasher --help
obdeect sources check run.toml
obdeect sources run run.toml -o out.pst
```

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`lightbulb;pst-color-primary` Source model
:link: sources
:link-type: doc
:class-card: sd-shadow-sm

Geometries, profiles, targets, frames and units, with examples.
:::

:::{grid-item-card} {fas}`terminal;pst-color-primary` Command line
:link: cli
:link-type: doc
:class-card: sd-shadow-sm

`obdeect sources run` and `check`.
:::

:::{grid-item-card} {fas}`sliders;pst-color-primary` Presets
:link: presets
:link-type: doc
:class-card: sd-shadow-sm

Ready-made setups: flashers, LEDs, point sources, laser beams, lidar, night sky,
muons and simple air showers.
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

Exact definitions, for implementers: run file, run records, sampling.
:::

::::

```{toctree}
:hidden:
:maxdepth: 1

sources
cli
presets
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Reference

spec/index
```
