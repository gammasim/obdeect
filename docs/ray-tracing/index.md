# Ray tracing

Optical ray tracing for gamma-ray astronomy, with a C++20 core, a Python
interface and a command line. Trace photons through mirrors and telescope
structure to the camera entrance, using a compiled optical model.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`bolt;pst-color-primary` User guide
:link: user-guide
:link-type: doc
:class-card: sd-shadow-sm

Install the ray-tracing package and run your first trace from the command line or Python.
:::

:::{grid-item-card} {fas}`eye;pst-color-primary` Optics
:link: optics
:link-type: doc
:class-card: sd-shadow-sm

Understand mirrors, photon paths, optical responses and shadowing.
:::

:::{grid-item-card} {fas}`check-circle;pst-color-primary` Validation
:link: VALIDATION
:link-type: doc
:class-card: sd-shadow-sm

Explore accuracy checks, comparisons with other codes, and measured performance.
:::

:::{grid-item-card} {fas}`book;pst-color-primary` Python API
:link: api
:link-type: doc
:class-card: sd-shadow-sm

Look up optical-model methods, compilation and analysis helpers.
:::

::::

```{toctree}
:hidden:
:maxdepth: 2

user-guide
optics
VALIDATION
photon-stream
api
```

```{toctree}
:caption: Notes
:hidden:
:maxdepth: 1

analysis-v0.7.2
```
