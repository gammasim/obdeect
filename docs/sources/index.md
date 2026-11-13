# Light sources

Artificial light sources for Cherenkov-telescope simulations: flashers, illuminators, lasers, LEDs and stars, emitting photons at sky, telescope-sphere, telescope or camera level, as a standard-library-only C++20 core with a Python API writing photon streams.

Status: alpha; the APIs may change before version 1.0. Light sources (profiles, stars, emitters, targets), the obdeect stream writer and all light-source programs listed under Programs are implemented.

Installation of obdeect sources and the build against the stream sources: [installation](../installation.md).

```{toctree}
:maxdepth: 1

sources
sim_telarray
spec/index
```

```{toctree}
:maxdepth: 1
:caption: Programs

programs/flasher
programs/led
programs/point-source
programs/beam-scan
programs/laser-beam
programs/night-sky
programs/muon-track
```
