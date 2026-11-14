# Specification

This section defines every light source and program of obdeect sources in physics and geometry
terms. It is written for implementers: with these pages and the obdeect stream headers, you can
reimplement obdeect sources. If you only want to use a source, read the
[source model](../sources.md) and the pages under Programs instead.

Two terms are used precisely. "Agree statistically" means agreement within four standard
deviations. "Exact" means identical counts, record order, metadata and deterministic columns.

```{toctree}
:maxdepth: 1

conventions
distributions
source-core
stream-writer
sky-sources
background
track
cli-conventions
tool-flasher
tool-led
tool-point-source
tool-beam-scan
tool-laser-beam
tool-night-sky
tool-muon-track
python-api
```

## Pages

| Page | Content |
| --- | --- |
| conventions | units, frames, random-number contract, float32 rule |
| distributions | density tables, pulses, spectra, angular profiles |
| source-core | `Source`: emitter and star geometry, targets, counts, validation |
| stream-writer | obdeect stream records, metadata keys, columns of every writer |
| sky-sources | sky-frame run/plan layer shared by the sky tools |
| background | night-sky source |
| track | muon-like Cherenkov track source |
| cli-conventions | option parsing shared by the tools |
| tool-* | one page per executable |
| python-api | dictionary API and plot helper |

## Build contract

- Header-only, standard-library-only C++20 core in `cpp/include/obdeect/sources/`; obdeect stream supplies
  the stream writer/reader, `PhotonChunk`, `uniform4` (Philox), stage ids and frame codes.
- Executables `obdeect sources flasher`, `-led`, `-point-source`, `-beam-scan`, `-laser-beam`,
  `-night-sky`, `-muon-track` from `cpp/tools/<name>.cpp`; a Python extension `obdeect.sources._core` (nanobind) and the package `obdeect.sources`.
