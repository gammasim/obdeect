# Ray tracing

The telescope turns parallel light into an image. This component traces photons through a
compiled **optical model**: the mirrors, the telescope structure that shadows them, the focal
surface and the optical responses such as mirror reflectivity. A photon ends at the camera
entrance, or is lost on the way, and the trace records where.

You can trace stars and point sources, calibration light and the photons of a photon stream.
The optical model comes from the simulation-models database, so the same code serves any telescope
that it describes. The stream step connects the atmosphere to sim_telarray's camera simulation.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`bolt;pst-color-primary` Command line
:link: usage
:link-type: doc
:class-card: sd-shadow-sm

Compile an optical model, trace a star and plot the result.
:::

:::{grid-item-card} {fas}`code;pst-color-primary` Python
:link: python-guide
:link-type: doc
:class-card: sd-shadow-sm

Trace photon arrays and analyse the arrivals.
:::

:::{grid-item-card} {fas}`lightbulb;pst-color-primary` Light sources
:link: light-sources
:link-type: doc
:class-card: sd-shadow-sm

Stars, point sources, flashers and air-shower photons as input.
:::

:::{grid-item-card} {fas}`stream;pst-color-primary` Photon streams
:link: photon-stream
:link-type: doc
:class-card: sd-shadow-sm

Trace a photon stream from the atmosphere to the camera.
:::

:::{grid-item-card} {fas}`eye;pst-color-primary` Optics
:link: optics
:link-type: doc
:class-card: sd-shadow-sm

Mirrors, photon paths, optical responses and shadowing.
:::

:::{grid-item-card} {fas}`check-circle;pst-color-primary` Validation
:link: VALIDATION
:link-type: doc
:class-card: sd-shadow-sm

Accuracy checks, comparisons with sim_telarray and performance.
:::

:::{grid-item-card} {fas}`book;pst-color-primary` Python API
:link: api
:link-type: doc
:class-card: sd-shadow-sm

Optical-model methods, compilation and analysis helpers.
:::

::::

```{toctree}
:hidden:
:maxdepth: 2
:caption: User guide

usage
python-guide
light-sources
photon-stream
batch-tracing
air-shower-photons
optical-models
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: Reference

optics
VALIDATION
api
```

```{toctree}
:caption: Notes
:hidden:
:maxdepth: 1

analysis-v0.7.2
```
