# Changelog

This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and uses semantic versioning after the first stable release.

<!-- towncrier release notes start -->

## [v0.1.0] - 2026-11-08

The first functional release of dlamb provides a C++20 photon-generation core,
a Python API, and command-line tools for producing gstrickd photon streams from
artificial and natural light sources. Sources can be emitted in the sky, telescope,
or camera frame, with configurable geometry and angular, spectral, and pulse
profiles. This release establishes the source models and stream-writing interface;
atmospheric propagation, telescope optics, and camera simulation remain outside
dlamb's scope.

### Added

- C++ core and Python API (`dlamb.source`) for generating photon chunks from point emitters and stars, with angular, spectral, and pulse profiles and disc or sphere targets. (#1)
- Command-line generators for flashers (`dlamb-flasher`, with dish and secondary types), LEDs (`dlamb-led`), point sources, beam scans, laser beams, night-sky background, and muon tracks. The tools write gstrickd photon streams. (#2, #3, #4)
- Direct camera-plane output for dish flashers via `dlamb-flasher --type dish --frame camera`, including straight-line propagation and the corresponding path-length time offset; documentation for using dlamb with sim_telarray. (#8)

### Maintenance

- Build against an installed gstrickd package so stream constants and format checks use the gstrickd definitions. (#7)
- Documented source geometry, units, options, and a physics-based specification of the source models, distributions, stream records, random-number keys, and test expectations. (#9, #12)
- Reimplemented the light sources against the specification and tests. (#13)
- Supported build platforms are Linux and macOS. (#11)


## [v0.0.1] - 2026-11-07

### Added

- Initial repository skeleton (build, packaging, CI, documentation); no functionality yet.
