# gstrickd Agent Guide

## Scope

`gstrickd` defines the photon stream exchanged between the modules of the simulation
chain (light sources, atmosphere, optics, camera): run, event and telescope records
with float32 photon columns, frames, MC truth and an efficiency status. It also holds
the only EventIO code of the chain: `corsika2ps` (CORSIKA IACT output to stream) and
`ps2simtel` (stream to sim_telarray input), plus `ps-cat` and `ps-tee`. Physics
(atmosphere, optics, light sources, camera) is out of scope.

## Layout

- `cpp/include/gstrickd/`: C++ core headers.
- `cpp/bindings/core.cpp`: nanobind extension; the only C++ file that may
  depend on Python.
- `cpp/tests/`: native tests (`test_*.cpp`) and shared helpers (`check.hpp`).
- `python/gstrickd/`: Python boundary (API, CLI, optional plotting).
- `python/tests/`: pytest suite.
- `docs/`: Sphinx sources. `docs/changes/` is towncrier input.

## Design rules

- Never ever modify the main README.md. Suggest instead changes, don't apply them.
- Keep the core C++20 standard-library-only. Python, plotting and file formats of
  other programs belong at the boundary, never in a per-photon kernel.
- Per-photon kernels are pure functions of (chunk, parameters): no allocation, I/O,
  Python calls or virtual dispatch per photon. Float32 SoA photon columns,
  counter-based random numbers (Philox) keyed on (seed, stage id, event, telescope,
  photon index). The scalar version is the reference for any SIMD or GPU version.
- Write kernels as templates on the scalar type (later differentiability).
- Make units, frames, random seeds and applied efficiencies explicit. Reject invalid
  input; never silently discard fields.
- Do not hardwire CTAO names, telescope geometry, model versions or physics
  parameters. Configuration comes from validated input written by simtools.
- Plotting is optional (`.[plot]`); matplotlib is never a core dependency.
- Use Python >=3.14 (as in `pyproject.toml`).
- The EventIO reader and writer are this project's own code. Never depend on
  sim_telarray's EventIO code; use sim_telarray only as a test reference.
- The stream format is versioned. Readers reject unknown major versions.

## Changes and validation

- Make the narrowest complete change and preserve unrelated work.
- Add a Towncrier fragment under `docs/changes/` for changes that should appear in
  the release notes.
- Add or update a focused test for every behavioural change. Analytic kernel tests
  come before cross-tool or performance tests.
- Results must not depend on photon ordering, chunk size or thread count.
- Do not weaken tolerances or replace reference data to hide a discrepancy.

## Build, test, and lint

```bash
python -m pip install -e '.[dev]'
cmake --preset debug
cmake --build --preset debug
ctest --test-dir build/debug --output-on-failure
python -m pytest python/tests
ruff format --check python
ruff check python
find cpp -type f \( -name '*.cpp' -o -name '*.hpp' \) -print0 | \
xargs -0 clang-format --dry-run --Werror
pre-commit run --all-files
```

Docs build with warnings as errors: `sphinx-build -W --keep-going -b html docs docs/_build/html`.
