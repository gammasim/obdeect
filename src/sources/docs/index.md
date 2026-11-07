# dlamb

Artificial light sources for Cherenkov-telescope simulations: flashers, illuminators, lasers, LEDs and stars, emitting photons at sky, telescope-sphere, telescope or camera level, as a standard-library-only C++20 core with a Python API writing gstrickd photon streams.

Status: pre-alpha. Light sources (profiles, stars, emitters, targets), the gstrickd stream writer and the four camera-level programs (ff-1m, ff-gct, fpls, pixled) are implemented.

```{toctree}
:maxdepth: 1

sources
changelog
```

```{toctree}
:maxdepth: 1
:caption: Programs

programs/ff-1m
programs/ff-gct
programs/fpls
programs/pixled
```
