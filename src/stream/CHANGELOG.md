# Changelog

This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and uses semantic versioning after the first stable release.

<!-- towncrier release notes start -->

## [v0.1.0] - 2026-11-08

First functional release. gstrickd defines the photon stream that connects the modules of a
Cherenkov-telescope simulation chain (light sources, atmosphere, optics, camera), and converts
CORSIKA output into this stream and the stream into sim_telarray and GrISU/GrOptics input.

### Added

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


## [v0.0.1] - 2026-11-07

### Added

- Initial repository skeleton (build, packaging, CI, documentation); no functionality yet.
