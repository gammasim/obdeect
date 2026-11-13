# obdeect sources

Artificial light sources for Cherenkov-telescope simulations: flashers, illuminators, lasers, LEDs and stars, emitting photons at sky, telescope-sphere, telescope or camera level, as a standard-library-only C++20 core with a Python API writing photon streams.

Status: alpha; the APIs may change before version 1.0. Light sources (profiles, stars, emitters, targets), the obdeect stream writer and all light-source programs listed under Programs are implemented.

**Building against obdeect stream.** obdeect sources builds against the stream component (obdeect stream headers and stage ids). CMake uses, in this order, the stream
target of the combined obdeect build, the stream sources in `../stream` of the repository, or an
installed obdeect stream package (`obdeect_stream_DIR`, then `python -m obdeect.stream --cmake-dir`). A source
distribution does not contain the stream sources, so it needs the installed package.

```{toctree}
:maxdepth: 1

sources
sim_telarray
spec/index
changelog
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
