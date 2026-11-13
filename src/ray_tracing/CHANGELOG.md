# Changelog

This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and uses semantic versioning after the first stable release.

<!-- towncrier release notes start -->

## [v0.7.2] - 2026-11-06

### Added

- Expanded the agent guide (`AGENTS.md`) with the repository layout, CI-equivalent sanitizer and FMA build commands, pre-commit coverage, and corrected paths and references. (#27)
- Add a --version to the obdeect executable. (#32)

### Fixed

- Accept two components for axes offsets (as expected by the simulation models). (#31)
- Include the Python headers before the standard library in the bindings, fixing source builds with Python 3.15 that failed on `_POSIX_C_SOURCE` redefinition warnings. Wheels are now also built for Python 3.15. (#33)

### Maintenance

- Reduce external python dependencies. (#28)
- Pull-request documentation builds only build the dev pages; releases are rebuilt on `main` and tags. Changelog fragments of type `maintenance` are now accepted. (#33)


## [v0.7.1] - 2026-11-04

Bug fix release.

### Fixed

- Support both components of the `axes_offsets` telescope parameter. The compiled optical model
  now includes the transverse offset, allowing models such as SSTS design v2.0.0 to run with
  obdeect. (#31)

## [v0.7.0] - 2026-11-04

Trace air-shower Cherenkov photons and pass them on to `sim_telarray` using the
`bypass_optics` mode.

### Added

- Air-shower photon bunches: `trace_bunches` takes `wavelength_law=(lo_nm, hi_nm)`, an exact
  1/λ² wavelength law with sim_telarray-style truncation to integer nm. It also takes
  `normalise_transmission`. Thinning sizes the output from the photons that survive.
- Telescope axes offsets are compiled into the optical model and applied when tracing. Input
  planes below an incoming obscurer are rejected, and `max_incoming_obscurer_z_m()` reports the
  obscurer height.
- New documentation: an air-shower photons guide, batch-tracing updates, and validation results
  in `docs/VALIDATION.md`.

### Changed

- Matplotlib is no longer a core dependency. Install `obdeect[plot]` for plotting; `obdeect plot`
  without it fails with a message saying so.
- Reduced external Python dependencies. `scipy` is no longer a dev dependency, and `astropy` is
  only needed for regenerating figures (`obdeect[figures]`).
- CI runs the Python lint through pre-commit, and the test job installs `.[plot]`.

## [v0.6.0]

### Changed

- The PyPI distribution is named `obdeect`. Install with `pip install obdeect`.
  The repository is hosted in the `gammasim` organisation.

## [v0.5.0]

Maintenance release: CI workflows and changelog. No change to the code or the Python API.

### Changed

- The documentation workflow uses `actions/deploy-pages` v5 (Node.js 24).
- The publishing workflow uses `actions/download-artifact` v8. A download whose hash does not
  match the expected hash now fails the workflow.
- The changelog is rebuilt from the git history for all releases from 0.1.0 on.

## [v0.4.0]

Performance release for the optical ray tracer, with additions to the Python API and the
command line. The optical-model format is unchanged (`obdeect.optical-model.v2`).

### Added

- `OpticalModel.trace` keyword arguments `outputs`, `detected_only`, `classify_losses` and `pool`.
- `ThreadPool` (persistent worker pool) and `trace_many` (several telescopes in one call).
- `trace_bunches`, tracing of CORSIKA-style photon bunches.
- `OpticalModel.obscurer_hits` and the command `obdeect obscurers`: the obscurers on the path
  of each photon, counted as `crossed`, `first` and `sole`.
- Optional fused multiply-add build, `-DOBDEECT_USE_FMA=ON` (x86-64, `-mfma`; off by default).
- CI job that builds and tests with `OBDEECT_USE_FMA=ON`.
- Documentation site (Sphinx, MyST Markdown) with the user guide, shadowing, telescope and optical
  element descriptions, light sources, batch tracing, physics limitations and the Python API.
  Figures show the telescope structure and the deviations from sim_telarray.
- Versioned documentation: development and tagged versions are published together to GitHub Pages,
  with a version switcher.

### Changed

- Obscurers are indexed per optical leg: a bounding-box hierarchy for every leg, and a uniform
  x-y grid for the incoming leg and for the legs after the primary, walked cell by cell along the
  ray.
- Aspheric mirrors stop at the single root of the ray when the slope bound guarantees one root.
  A root in a segment gap is a miss, without scanning the aperture. Hits are unchanged.
- Segment membership on aspheric mirrors uses a cell grid with cached verdicts.
- Asphere sag and slope are evaluated as polynomials in r² (rounding-level change).
- `telescope_transmission` is checked by `response_bound`, which raises `ModelError` for an
  unsupported transmission instead of accepting it.
- Replay comparison has tiers, used to check results across builds.
- Source comments and docstrings describe obdeect's own conventions instead of sim_telarray's.
  No change in behaviour.

### Fixed

- Secondary baffle: a fourth value not below R1 is read as the outer radius, as in sim_telarray,
  instead of as the wall thickness. No current model uses this value.
- Bunch acceptance rejects axes that are not finite or not increasing; wavelength density nodes
  must start above 0 nm.
- A failure while starting the thread pool no longer leaves joinable threads behind; an
  exception in a task is rethrown to the caller and the remaining chunks are skipped.

## [v0.3.0]

Phase 2: calibration light.

### Added

- Point emitters: `emitter_photons` makes weighted photons towards a target disk, with Poisson
  statistics.
- `trace_emitter` and `trace_emitter_blocks`, which trace emitter light in blocks of 250,000
  photons, so memory stays bounded.
- Direct tracing of photons that start inside the telescope (flashers), with `trace(direct=True)`.
- Documentation of the hand-off of photon bunches to sim_telarray for dual-mirror telescopes.

### Changed

- Segment footprint frames are computed once per optical model. Results are bit-identical; SST
  tracing is about 1.7 times faster.

## [v0.2.1]

### Added

- Linux arm64 wheels.

## [v0.2.0]

Phase 1: optical models for in-process use by simtools, with sim_telarray equivalence.

### Added

- Python API and `obdeect` command line with `compile`, `trace` and `plot`.
- Compilation of CTAO optical models from simulation-models: mirror lists, telescope structure
  (one table with explicit shapes), telescope transmission, camera transport, patch productions.
- Photon diagnostics and CTAO validation scripts.
- Deterministic optical transport.

### Changed

- Detector and ray-tracing seeds are separate.
- Legacy optical-model seeds are read compatibly.

### Fixed

- Optical input validation and deterministic boundary selection.
- GCC header portability and the macOS wheel deployment target.

## [v0.1.0]

First version: a C++20 ray tracer with a Python package.

### Added

- C++20 artificial 400-nm source, spherical primary and focal-screen tracer.
- Simple MST-inspired camera and four-mast pre-primary shadow model.
- CSV trace-path output and optional Matplotlib diagnostic plot.
- CMake, Python packaging, tests, linting CI, BSD-3-Clause license and citation metadata.
