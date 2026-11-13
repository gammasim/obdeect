# Installation

Install all components, or only the ones you need, in a virtual environment:

```sh
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install obdeect              # all components at one version
python -m pip install obdeect-sources      # one component (also installs obdeect-stream)
obdeect --version                          # installed components and stream format version
```

| Component | Package | Python | Command |
|---|---|---|---|
| photon stream | `obdeect-stream` | `obdeect.stream` | `obdeect stream ...` |
| atmospheric propagation | `obdeect-propagation` | `obdeect.propagation` | `obdeect propagation ...` |
| light sources | `obdeect-sources` | `obdeect.sources` | `obdeect sources ...` |
| ray tracing | `obdeect-ray-tracing` | `obdeect.ray_tracing` | `obdeect ray-tracing ...` |

The packages require Python 3.14 or later. If no wheel is available for your platform,
installation builds the extensions and requires a C++20 compiler and CMake 3.20 or later.
For plots, install the `plot` extra, e.g. `python -m pip install "obdeect-ray-tracing[plot]"`.

## From source

From a checkout, build and test all components:

```sh
cmake --preset release && cmake --build --preset release && ctest --preset release
python -m pip install -e ./src/stream -e ./src/propagation -e ./src/sources -e "./src/ray_tracing[dev]"
for c in stream propagation sources ray_tracing; do python -m pytest src/$c/python/tests; done
```

## C++ consumers

The stream component installs its headers and a CMake package. Other components build
against, in this order, the stream target of the combined build, the stream sources in
`src/stream`, the copy bundled in a component's source distribution (`stream-src`), or an
installed package (`obdeect_stream_DIR`, then `python -m obdeect.stream --cmake-dir`). A
source distribution of a component carries the stream sources of its own version, so it
builds in an isolated environment; `obdeect-stream` of the same version remains the runtime
dependency.

## Documentation

```sh
python -m pip install "./src/ray_tracing[docs]"
sphinx-build -W --keep-going -b html docs docs/_build/html
```
