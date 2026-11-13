# Changelog

All components are released together with one version (tags `vX.Y.Z`); the sections below keep the history of the former packages.

<!-- towncrier release notes start -->

## Ray tracing (released as obdeect up to v0.7.2)

### [v0.7.2] - 2026-11-06

#### Added

- Expanded the agent guide (`AGENTS.md`) with the repository layout, CI-equivalent sanitizer and FMA build commands, pre-commit coverage, and corrected paths and references. (#27)
- Add a --version to the obdeect executable. (#32)

#### Fixed

- Accept two components for axes offsets (as expected by the simulation models). (#31)
- Include the Python headers before the standard library in the bindings, fixing source builds with Python 3.15 that failed on `_POSIX_C_SOURCE` redefinition warnings. Wheels are now also built for Python 3.15. (#33)

#### Maintenance

- Reduce external python dependencies. (#28)
- Pull-request documentation builds only build the dev pages; releases are rebuilt on `main` and tags. Changelog fragments of type `maintenance` are now accepted. (#33)


### [v0.7.1] - 2026-11-04

Bug fix release.

#### Fixed

- Support both components of the `axes_offsets` telescope parameter. The compiled optical model
  now includes the transverse offset, allowing models such as SSTS design v2.0.0 to run with
  obdeect. (#31)

### [v0.7.0] - 2026-11-04

Trace air-shower Cherenkov photons and pass them on to `sim_telarray` using the
`bypass_optics` mode.

#### Added

- Air-shower photon bunches: `trace_bunches` takes `wavelength_law=(lo_nm, hi_nm)`, an exact
  1/λ² wavelength law with sim_telarray-style truncation to integer nm. It also takes
  `normalise_transmission`. Thinning sizes the output from the photons that survive.
- Telescope axes offsets are compiled into the optical model and applied when tracing. Input
  planes below an incoming obscurer are rejected, and `max_incoming_obscurer_z_m()` reports the
  obscurer height.
- New documentation: an air-shower photons guide, batch-tracing updates, and validation results
  in `docs/VALIDATION.md`.

#### Changed

- Matplotlib is no longer a core dependency. Install `obdeect[plot]` for plotting; `obdeect plot`
  without it fails with a message saying so.
- Reduced external Python dependencies. `scipy` is no longer a dev dependency, and `astropy` is
  only needed for regenerating figures (`obdeect[figures]`).
- CI runs the Python lint through pre-commit, and the test job installs `.[plot]`.

### [v0.6.0]

#### Changed

- The PyPI distribution is named `obdeect`. Install with `pip install obdeect`.
  The repository is hosted in the `gammasim` organisation.

### [v0.5.0]

Maintenance release: CI workflows and changelog. No change to the code or the Python API.

#### Changed

- The documentation workflow uses `actions/deploy-pages` v5 (Node.js 24).
- The publishing workflow uses `actions/download-artifact` v8. A download whose hash does not
  match the expected hash now fails the workflow.
- The changelog is rebuilt from the git history for all releases from 0.1.0 on.

### [v0.4.0]

Performance release for the optical ray tracer, with additions to the Python API and the
command line. The optical-model format is unchanged (`obdeect.optical-model.v2`).

#### Added

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

#### Changed

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

#### Fixed

- Secondary baffle: a fourth value not below R1 is read as the outer radius, as in sim_telarray,
  instead of as the wall thickness. No current model uses this value.
- Bunch acceptance rejects axes that are not finite or not increasing; wavelength density nodes
  must start above 0 nm.
- A failure while starting the thread pool no longer leaves joinable threads behind; an
  exception in a task is rethrown to the caller and the remaining chunks are skipped.

### [v0.3.0]

Phase 2: calibration light.

#### Added

- Point emitters: `emitter_photons` makes weighted photons towards a target disk, with Poisson
  statistics.
- `trace_emitter` and `trace_emitter_blocks`, which trace emitter light in blocks of 250,000
  photons, so memory stays bounded.
- Direct tracing of photons that start inside the telescope (flashers), with `trace(direct=True)`.
- Documentation of the hand-off of photon bunches to sim_telarray for dual-mirror telescopes.

#### Changed

- Segment footprint frames are computed once per optical model. Results are bit-identical; SST
  tracing is about 1.7 times faster.

### [v0.2.1]

#### Added

- Linux arm64 wheels.

### [v0.2.0]

Phase 1: optical models for in-process use by simtools, with sim_telarray equivalence.

#### Added

- Python API and `obdeect` command line with `compile`, `trace` and `plot`.
- Compilation of CTAO optical models from simulation-models: mirror lists, telescope structure
  (one table with explicit shapes), telescope transmission, camera transport, patch productions.
- Photon diagnostics and CTAO validation scripts.
- Deterministic optical transport.

#### Changed

- Detector and ray-tracing seeds are separate.
- Legacy optical-model seeds are read compatibly.

#### Fixed

- Optical input validation and deterministic boundary selection.
- GCC header portability and the macOS wheel deployment target.

### [v0.1.0]

First version: a C++20 ray tracer with a Python package.

#### Added

- C++20 artificial 400-nm source, spherical primary and focal-screen tracer.
- Simple MST-inspired camera and four-mast pre-primary shadow model.
- CSV trace-path output and optional Matplotlib diagnostic plot.
- CMake, Python packaging, tests, linting CI, BSD-3-Clause license and citation metadata.

## Photon stream (released as gstrickd)

### [v0.1.0] - 2026-11-08

First functional release. gstrickd defines the photon stream that connects the modules of a
Cherenkov-telescope simulation chain (light sources, atmosphere, optics, camera), and converts
CORSIKA output into this stream and the stream into sim_telarray and GrISU/GrOptics input.

#### Added

- Photon-stream format 1.0 (`docs/format.md`): run, event and telescope records; photons as
  float32 columns (position, direction, time, emission height, weight, wavelength) in chunks of
  at most 65536 rows; frames sky, ground, telescope and camera; MC truth of runs and events; an
  efficiency status that records for each factor (atmosphere, QE, mirror reflectivity, ...)
  whether it is applied, pending or applied as an envelope. Readers reject unknown major
  versions.
- Reader and writer in C++20 (header-only, standard library only, `gstrickd/stream.hpp`), in
  Python (`gstrickd.Reader`, `gstrickd.Writer`, photon columns as NumPy arrays) and in C
  (`gstrickd/gstrickd.h`, `libgstrickd_c`) for programs in other languages;
  `gstrickd.stage.transform` for writing a processing stage in Python; conformance streams in
  `cpp/tests/data/conformance/` for testing independent implementations.
- Counter-based random numbers (Philox) keyed on seed, stage, event, telescope and photon, and
  the stage ids of all chain modules, so that results do not depend on thread count or chunk
  size.
- `corsika2ps`: CORSIKA IACT EventIO files (long, compact and split photon blocks) to a photon
  stream. The EventIO reader and writer are part of gstrickd.
- `ps2simtel`: photon stream to sim_telarray input for ground-, camera- and telescope-frame
  photons, with a camera map per telescope and the required `Bypass_Optics` value reported.
  It checks the efficiency status, so that no factor is applied twice or not at all
  (`--allow-pending`; `--scaled-factor qe` when the QE is thinned before ray tracing), writes
  the run and event headers for streams from light sources, and an atmospheric profile with
  `--atmprof`. More than 1 GB/s through a pipe.
- `ps2grisu`: ground-frame photon streams as GrISU/GrOptics photon text input, with telescope
  matching from a `TLLOC` file.
- `ps-cat` (summary and concatenation of streams) and `ps-tee` (one stream to several
  consumers). Every tool reports its version and the stream-format version with `--version`.
- Wheels for Linux (x86-64, ARM) and macOS with the C++ headers, the C library and a CMake
  package (`find_package(gstrickd 0.1)`; `python -m gstrickd --cmake-dir`, `--include-dir`,
  `--lib-dir`). Supported platforms are Linux and macOS.


### [v0.0.1] - 2026-11-07

#### Added

- Initial repository skeleton (build, packaging, CI, documentation); no functionality yet.

## Light propagation (released as dlufd)

### [v0.1.0] - 2026-11-08

First functional release of dlufd's atmospheric photon-propagation tools.

#### Added

- `dlufd-propagate` and `dlufd.propagate` process gstrickd photon streams, moving sky-frame photons to the telescope plane and applying atmospheric effects. (#4)
- Atmospheric extinction from transmission tables, wavelength assignment for photons without a wavelength, and early thinning with an efficiency envelope. Extinction can include the photon path inside each telescope, with a path length set per telescope. C++ and Python APIs are available. (#1, #7)
- Refraction and air-travel-time calculations from an `atmprof` refractive-index profile, including group-index timing. (#2)
- Optional single scattering of light-source photons using Rayleigh or Henyey-Greenstein scattering. Scattered photons are directed towards telescope spheres and use an independent random stream. Scattering is off by default. (#3, #5)
- Geometry-only propagation with `--no-extinction` (`extinction=False` in Python): assign missing wavelengths, refract photons and calculate air delay without applying extinction or thinning. Photon weights are preserved so a later stage can apply the atmosphere. (#9)
- `dlufd-convert-table` and `dlufd.convert_table` convert supported legacy extinction tables to dlufd's transmission-table format. (#14)

#### Maintenance

- Photon-stream handling uses gstrickd's shared format and atmosphere-stage identifier. CMake finds the installed package before fetching it. (#8)
- Linux and macOS are supported; Windows is not supported. (#12)


### [v0.0.1] - 2026-11-07

#### Added

- Initial repository skeleton (build, packaging, CI, documentation); no functionality yet.

## Light sources (released as dlamb)

### [v0.1.0] - 2026-11-08

The first functional release of dlamb provides a C++20 photon-generation core,
a Python API, and command-line tools for producing gstrickd photon streams from
artificial and natural light sources. Sources can be emitted in the sky, telescope,
or camera frame, with configurable geometry and angular, spectral, and pulse
profiles. This release establishes the source models and stream-writing interface;
atmospheric propagation, telescope optics, and camera simulation remain outside
dlamb's scope.

#### Added

- C++ core and Python API (`dlamb.source`) for generating photon chunks from point emitters and stars, with angular, spectral, and pulse profiles and disc or sphere targets. (#1)
- Command-line generators for flashers (`dlamb-flasher`, with dish and secondary types), LEDs (`dlamb-led`), point sources, beam scans, laser beams, night-sky background, and muon tracks. The tools write gstrickd photon streams. (#2, #3, #4)
- Direct camera-plane output for dish flashers via `dlamb-flasher --type dish --frame camera`, including straight-line propagation and the corresponding path-length time offset; documentation for using dlamb with sim_telarray. (#8)

#### Maintenance

- Build against an installed gstrickd package so stream constants and format checks use the gstrickd definitions. (#7)
- Documented source geometry, units, options, and a physics-based specification of the source models, distributions, stream records, random-number keys, and test expectations. (#9, #12)
- Reimplemented the light sources against the specification and tests. (#13)
- Supported build platforms are Linux and macOS. (#11)


### [v0.0.1] - 2026-11-07

#### Added

- Initial repository skeleton (build, packaging, CI, documentation); no functionality yet.
