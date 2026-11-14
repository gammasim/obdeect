# Specification

This section specifies, in physics and geometry terms, every light source and tool of obdeect sources. An
implementer who has only these pages, the tests under `cpp/tests` and `python/tests`, and the
obdeect stream headers must be able to reimplement obdeect sources. "Agree statistically" means agreement within
the four-sigma checks of the tests; "exact" means bit-identical counts, record order, metadata
and deterministic columns.

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
  `-night-sky`, `-muon-track` from `cpp/tools/<name>.cpp`; tests `cpp/tests/test_<name>.cpp`
  (one CTest each); Python extension `obdeect.sources._core` (nanobind) and package `python/obdeect sources`.
