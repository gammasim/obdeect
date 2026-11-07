# dlufd

Atmospheric propagation for Cherenkov-telescope simulations: extinction, refraction and timing from the refractive-index profile, optional single scattering, and early photon thinning, as a standard-library-only C++20 core with a Python API, working on gstrickd photon streams.

Status: pre-alpha. Extinction, refraction and timing are implemented; the atmosphere stage on photon streams is in [propagation](propagation.md). Single scattering is not.

```{toctree}
:maxdepth: 1

changelog
extinction
refraction
propagation
```
