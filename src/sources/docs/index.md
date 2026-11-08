# dlamb

Artificial light sources for Cherenkov-telescope simulations: flashers, illuminators, lasers, LEDs and stars, emitting photons at sky, telescope-sphere, telescope or camera level, as a standard-library-only C++20 core with a Python API writing gstrickd photon streams.

Status: pre-alpha. Light sources (profiles, stars, emitters, targets), the gstrickd stream writer and all light-source programs listed under Programs are implemented.

**Building against gstrickd.** CMake looks for the installed gstrickd package
(`find_package(gstrickd 0.1 CONFIG)`) in this order: `gstrickd_DIR` or `CMAKE_PREFIX_PATH`,
then the package of the Python interpreter (`python -m gstrickd --cmake-dir`). If neither is
found, it fetches a pinned gstrickd commit from GitHub. The configure output names the source
used. Use `-DDLAMB_REQUIRE_INSTALLED_GSTRICKD=ON` to fail instead of fetching, as for an offline
build; a source distribution does not contain gstrickd.

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
