# dlufd

Atmospheric propagation for Cherenkov-telescope simulations: extinction, refraction and timing from the refractive-index profile, optional single scattering, and early photon thinning, as a standard-library-only C++20 core with a Python API, working on gstrickd photon streams.

Status: pre-alpha. Refraction and timing are implemented (see [refraction](refraction.md)); extinction and optional single scattering (see [scattering](scattering.md)) are implemented; wavelength assignment and thinning are part of extinction.

```{toctree}
:maxdepth: 1

changelog
extinction
refraction
scattering
```
