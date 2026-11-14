# Photon stream

A photon stream is the format that connects the modules of obdeect. It is a binary stream of
records: one run, its events, the telescopes of each event, and the photons that reach each
telescope in chunks of columns (position, direction, time, weight and wavelength). The stream
also records which efficiency factors each module has applied, so that the next module knows
what is left to do. Streams are written and read in sequence, so modules can run side by side in a pipe.

The component provides the format with readers and writers in C++, Python and C, and the
converters at both ends of the chain. `obdeect stream from-corsika` reads CORSIKA IACT output,
and `obdeect stream to-simtel` writes input for sim_telarray, the camera simulation. For the
GrISU chain, `obdeect stream to-grisu` writes the photon text read by GrOptics, which feeds
GrISUDet. Two tools, `cat` and `tee`, join and copy streams.

Every other component builds on this one: all use the same frames, units, random-number rules
and efficiency records, which are defined here once.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`terminal;pst-color-primary` Command line
:link: command-line
:link-type: doc
:class-card: sd-shadow-sm

Join, copy and convert streams; connect to sim_telarray.
:::

:::{grid-item-card} {fas}`table-list;pst-color-primary` Format
:link: format
:link-type: doc
:class-card: sd-shadow-sm

Records, metadata, photon columns and frames.
:::

:::{grid-item-card} {fas}`scale-balanced;pst-color-primary` Efficiency status
:link: efficiency-status
:link-type: doc
:class-card: sd-shadow-sm

How modules record which efficiencies they have applied.
:::

:::{grid-item-card} {fas}`dice;pst-color-primary` Random numbers
:link: random-numbers
:link-type: doc
:class-card: sd-shadow-sm

One generator, reproducible results, independent modules.
:::

:::{grid-item-card} {fas}`code;pst-color-primary` Programming interfaces
:link: programming
:link-type: doc
:class-card: sd-shadow-sm

Read and write streams from C++, Python and C.
:::

:::{grid-item-card} {fas}`file-export;pst-color-primary` GrOptics and GrISU
:link: grisu
:link-type: doc
:class-card: sd-shadow-sm

Export photons as GrISU photon text.
:::

::::

```{toctree}
:hidden:
:maxdepth: 1

command-line
format
efficiency-status
random-numbers
programming
grisu
```
