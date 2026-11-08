# Specification

This section specifies, in physics and geometry terms, every light source and tool of dlamb. An
implementer who has only these pages, the tests under `cpp/tests` and `python/tests`, and the
gstrickd headers must be able to reimplement dlamb. "Agree statistically" means agreement within
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
tests
```

## Pages

| Page | Content |
| --- | --- |
| conventions | units, frames, random-number contract, float32 rule |
| distributions | density tables, pulses, spectra, angular profiles |
| source-core | `Source`: emitter and star geometry, targets, counts, validation |
| stream-writer | gstrickd records, metadata keys, columns of every writer |
| sky-sources | sky-frame run/plan layer shared by the sky tools |
| background | night-sky source |
| track | muon-like Cherenkov track source |
| cli-conventions | option parsing shared by the tools |
| tool-* | one page per executable |
| python-api | dictionary API and plot helper |
| tests | list of tests, the C++ and Python API they require, acceptance criteria |

## Build contract

- Header-only, standard-library-only C++20 core in `cpp/include/dlamb/`; gstrickd supplies
  the stream writer/reader, `PhotonChunk`, `uniform4` (Philox), stage ids and frame codes.
- Executables `dlamb-flasher`, `-led`, `-point-source`, `-beam-scan`, `-laser-beam`,
  `-night-sky`, `-muon-track` from `cpp/tools/<name>.cpp`; tests `cpp/tests/test_<name>.cpp`
  (one CTest each); Python extension `dlamb._core` (nanobind) and package `python/dlamb`.
- Header and namespace names used by the tests are listed on the [tests](tests.md) page.
