# Atmospheric propagation

Extinction, refraction, timing and optional single scattering for Cherenkov-telescope photon
streams. C++20 core (standard library only) with a Python API.

Status: alpha. APIs may change before version 1.0.

- [Extinction](extinction.md): transmission, wavelengths, thinning.
- [Refraction](refraction.md): ray bending and air delay.
- [Propagation](propagation.md): the atmosphere stage on photon streams (`obdeect propagation run`).
- [Scattering](scattering.md): single scattering of light-source photons.

```{toctree}
:maxdepth: 1

extinction
refraction
propagation
scattering
```
