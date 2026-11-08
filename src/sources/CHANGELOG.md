# Changelog

This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and uses semantic versioning after the first stable release.

<!-- towncrier release notes start -->

## [v0.1.0] - 2026-11-08

### Added

- Light sources: stars and point emitters with angular, spectral and pulse profiles and disc or sphere targets, emitted as gstrickd photon chunks (C++ core and Python API `dlamb.source`). (#1)
- Add dlamb tools writing gstrickd photon streams: `dlamb-flasher` (dish and secondary flashers) and `dlamb-led`. (#2)
- Sky-level light sources `dlamb-point-source`, `dlamb-beam-scan`, `dlamb-laser-beam`, `dlamb-night-sky` write gstrickd photon streams (sky frame); night-sky background source in `dlamb/background.hpp`. (#3)
- Add the light-source tool `dlamb-muon-track`, writing gstrickd photon streams. Light-source tools have descriptive names (`dlamb-flasher`, `dlamb-point-source`, ...). (#4)
- One `dlamb-flasher` (`--type dish|secondary`) replaces the two flasher tools; `dlamb-pixel-led` is now `dlamb-led`; the figure-eight, rotating-LED and geometry-report demonstration tools are removed. (#6)
- `dlamb-flasher --type dish --frame camera` writes the flasher photons directly on the camera plane (frame C, straight-line propagation, time includes path / c); new page "Using dlamb with sim_telarray". (#8)

### Maintenance

- Build against an installed gstrickd package (`find_package(gstrickd)`, or the package of the Python interpreter); dlamb's stage ids and format checks use gstrickd's constants, and the Python copies of frame, column and chunk constants are tested against the compiled gstrickd. (#7)
- Documentation, comments and help texts describe dlamb's own geometry, units and options; each tool page names the earlier program it replaces in one sentence. (#9)
- Windows is no longer supported: builds are tested and supported on Linux and macOS only. (#11)
- Add a specification of all light sources under docs/spec: options, geometry, distributions, stream records, random-number keys and test expectations, written in physics terms for a clean-room reimplementation. (#12)
- Light sources reimplemented from the physics specification (`docs/spec`) and the tests, without reference to the earlier code. (#13)


## [v0.0.1] - 2026-11-07

### Added

- Initial repository skeleton (build, packaging, CI, documentation); no functionality yet.
