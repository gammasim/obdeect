# obdeect Agent Guide

## Scope

`obdeect` is one project with four components, released together under one version
(tags `vX.Y.Z`). Each component is a portable C++20 core with a thin Python interface:

- `src/stream`: photon-stream format, EventIO to and from CORSIKA and sim_telarray.
- `src/propagation`: photon propagation through the atmosphere.
- `src/sources`: artificial light sources.
- `src/ray_tracing`: optical ray tracer for telescope optics and structure.

Air-shower generation, electronics, triggers and event writing are out of scope.
Read the relevant code and tests first. Use the component `docs/` for architecture
and science changes.

## Layout

Per component `src/<c>/`:

- `cpp/include/obdeect/<c>/`: C++ core headers.
- `cpp/bindings/core.cpp`: nanobind extension; the only C++ file that may depend on Python.
- `cpp/tests/`: native tests (`test_*.cpp`), helpers (`check.hpp`), fixtures in `cpp/tests/data/`.
- `python/obdeect/<c>/`: Python boundary (API, CLI, optional plotting).
- `python/tests/`: pytest suite.
- `docs/`: Sphinx sources.

Python packages: `obdeect.stream`, `.propagation`, `.sources`, `.ray_tracing`.
Wheels: `obdeect-stream`, `obdeect-propagation`, `obdeect-sources`, `obdeect-ray-tracing`.
`pip install obdeect` installs all. CLI: `obdeect <component> <command>`.
Root `docs/changes/` is the Towncrier input for all components.

## Design rules

- Never modify the main `README.md`. Suggest changes instead, do not apply them.
- Terminology: call a compiled telescope geometry and optical-response artifact an "optical model".
  In documentation, say "efficiency status", never "ledger"; say "photons in the sky frame", not
  "sky photons"; call the `--envelope` tables "efficiency bounds". Photon-stream files use `.pst`,
  never `.ps` (PostScript).
- Keep it simple, efficient, deterministic and portable. Prefer small, explicit data
  structures and testable kernels over frameworks, global state, hidden defaults or
  abstraction in the hot path.
- Keep the core C++20 standard-library-only. Python, plotting, file formats of other
  programs and simulator dependencies belong at the boundary, never in a kernel.
- Per-photon kernels are pure functions of (chunk, parameters): no allocation, I/O,
  Python calls or virtual dispatch per photon. Use float32 SoA photon columns (or
  immutable compiled models and contiguous buffers where stated by the component).
  Write kernels as templates on the scalar type. The scalar version is the reference
  for any SIMD or GPU version.
- Random numbers: counter-based (Philox), keyed on (seed, stage id, event, telescope, photon index).
- Make units, frames, random seeds, statuses, efficiencies and loss accounting explicit.
  Reject invalid or unsupported input; never silently discard fields.
- Do not hardwire CTAO names, telescope geometry, model versions, file paths or physics
  parameters into the core, the public API or defaults. Configuration comes from
  validated input data (simtools writes it). CTAO support is an explicit importer or
  versioned fixture only, never required and never an implicit fallback.
- Plotting is optional (`.[plot]`); matplotlib is never a core dependency.
- Use Python >=3.14 (as in the root `pyproject.toml`). Do not fall back to earlier versions.
- The stream format is versioned; readers reject unknown major versions.
- The EventIO reader and writer are this project's own code. Use sim_telarray only as a test reference.

## Component rules

- `stream`: Results must not depend on photon ordering, chunk size or thread count.
  Record each efficiency in the stream.
- `propagation`: Atmospheric effects are an ordered list of independent effects, each
  recorded in the efficiency status of the stream. Extinction follows `docs/extinction.md`;
  its tests hold the reference values.
- `sources`: Sources are defined by validated input (position, direction, angular,
  spectral and pulse profiles), never by hardwired device names. Sources emit at sky,
  telescope-sphere (G), telescope (T) or camera (C) level.
- `ray_tracing`: Use immutable compiled optical models and contiguous SoA buffers. After
  output buffers are sized, the per-photon path must not allocate, do I/O, call Python or
  dispatch virtually. Terminology: "optical model", not telescope-specific names.
  When renaming a public concept, update implementation symbols, imports, exception
  types, serialized keys, fixture variables and test call sites together.
  Use `docs/VALIDATION.md` and `docs/optical-models.md` for architectural or scientific changes.

## Changes and validation

- Make the narrowest complete change and preserve unrelated work.
- Add a Towncrier fragment under root `docs/changes/` for changes that belong in the release notes.
- Add or update a focused test for every behavioural change. Analytic kernel and primitive
  tests come before optical-model, cross-tool or performance tests.
- Results must not depend on photon ordering, block or chunk size, or thread count.
- Do not weaken tolerances or replace reference data to hide a discrepancy. Record
  provenance for fixtures and cross-tool comparisons. Keep examples and fixtures labelled
  and configuration-driven.

## Build, test and lint

Run from the repository root. One combined build covers all components:

```bash
python -m pip install -e '.[dev]'     # editable install; rebuilds on C++ changes
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
python -m pytest src/<component>/python/tests   # e.g. src/ray_tracing
ruff format --check src
ruff check src
find src -type f \( -name '*.cpp' -o -name '*.hpp' \) -print0 | \
xargs -0 clang-format --dry-run --Werror
pre-commit run --all-files
```

Ruff uses a 100-character line limit (root `pyproject.toml`). Run `ruff format` before
`ruff check`. A clean handoff needs both `ruff format --check src` and `ruff check src`.
After renames, grep the Python files for removed symbols and run the focused tests first.

Run `clang-tidy -p build <changed-file>` for changed C++ files when available.
`pre-commit run --all-files` runs ruff, codespell, markdownlint, yamllint, actionlint,
pyproject-fmt, shellcheck and the towncrier check. Run the full checks for shared kernels,
interfaces, configuration or tooling changes.

CI also builds a `-DOBDEECT_USE_FMA=ON` Release and an AddressSanitizer/UBSan Debug build.
Run them locally for floating-point kernel or memory-handling changes:

```bash
cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-sanitize --parallel
ctest --test-dir build-sanitize --output-on-failure
```

Docs build with warnings as errors, per component:
`sphinx-build -W --keep-going -b html src/<component>/docs src/<component>/docs/_build/html`.
Regenerate ray-tracing figures only from new validation results, via `src/ray_tracing/docs/figures/make_figures.py`.
