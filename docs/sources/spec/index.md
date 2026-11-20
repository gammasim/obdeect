# Specification

This section defines the light sources and the run in physics and geometry terms. It is written
for implementers: with these pages and the obdeect stream headers, you can reimplement obdeect
sources. If you only want to use a source, read the [source model](../sources.md), the
[command line](../cli.md) and the [presets](../presets.md) instead.

Two terms are used precisely. "Agree statistically" means agreement within four standard
deviations. "Exact" means identical counts, record order, metadata and deterministic columns.

```{toctree}
:maxdepth: 1

conventions
distributions
source-core
run
run-file
python-api
```

## Pages

| Page | Content |
| --- | --- |
| conventions | units, frames, random-number contract, float32 rule |
| distributions | density tables, pulses, spectra, angular profiles |
| source-core | geometries, targets, validation, counts, sampling |
| run | run records and metadata, seeds, bunches, auto target, cache |
| run-file | every key of the run file with type and default |
| python-api | configuration, run and emission functions |

## Build contract

- Header-only, standard-library-only C++20 core in `src/sources/cpp/include/obdeect/sources/`;
  obdeect stream supplies the photon-stream writer, `PhotonChunk`, `uniform4` (Philox), stage ids
  and frame codes.
- Command `obdeect sources run|check` and the Python package `obdeect.sources`; the compiled
  extension `obdeect.sources._core` (nanobind) is the only C++ binding.
