# Changelog

This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and uses semantic versioning after the first stable release.

<!-- towncrier release notes start -->

## [v0.1.0] - 2026-11-08

First functional release of dlufd's atmospheric photon-propagation tools.

### Added

- `dlufd-propagate` and `dlufd.propagate` process gstrickd photon streams, moving sky-frame photons to the telescope plane and applying atmospheric effects. (#4)
- Atmospheric extinction from transmission tables, wavelength assignment for photons without a wavelength, and early thinning with an efficiency envelope. Extinction can include the photon path inside each telescope, with a path length set per telescope. C++ and Python APIs are available. (#1, #7)
- Refraction and air-travel-time calculations from an `atmprof` refractive-index profile, including group-index timing. (#2)
- Optional single scattering of light-source photons using Rayleigh or Henyey-Greenstein scattering. Scattered photons are directed towards telescope spheres and use an independent random stream. Scattering is off by default. (#3, #5)
- Geometry-only propagation with `--no-extinction` (`extinction=False` in Python): assign missing wavelengths, refract photons and calculate air delay without applying extinction or thinning. Photon weights are preserved so a later stage can apply the atmosphere. (#9)
- `dlufd-convert-table` and `dlufd.convert_table` convert supported legacy extinction tables to dlufd's transmission-table format. (#14)

### Maintenance

- Photon-stream handling uses gstrickd's shared format and atmosphere-stage identifier. CMake finds the installed package before fetching it. (#8)
- Linux and macOS are supported; Windows is not supported. (#12)


## [v0.0.1] - 2026-11-07

### Added

- Initial repository skeleton (build, packaging, CI, documentation); no functionality yet.
