# Ray tracing

Optical ray tracing for gamma-ray astronomy. It traces photons through mirrors and
telescope structure to the camera entrance, using a compiled optical model. The
core is C++20, with a Python module and a command line.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`bolt;pst-color-primary` User guide
:link: user-guide
:link-type: doc
:class-card: sd-shadow-sm

Install, compile an optical model and trace your first photons.
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

Accuracy checks, comparisons with sim_telarray and measured performance.
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
