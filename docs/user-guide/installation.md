# Installation

obdeect needs Python 3.14 or later. Install everything, or only the components you need,
in a virtual environment:

```sh
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install obdeect              # all components, one version
python -m pip install obdeect-sources      # one component (also installs obdeect-stream)
obdeect --version                          # installed components and photon-stream format version
```

| Component | Package | Python module | Command |
|---|---|---|---|
| photon stream | `obdeect-stream` | `obdeect.stream` | `obdeect stream ...` |
| atmospheric propagation | `obdeect-propagation` | `obdeect.propagation` | `obdeect propagation ...` |
| light sources | `obdeect-sources` | `obdeect.sources` | `obdeect sources ...` |
| ray tracing | `obdeect-ray-tracing` | `obdeect.ray_tracing` | `obdeect ray-tracing ...` |

All components share one version. If no wheel exists for your platform, pip builds the
extensions and needs a C++20 compiler and CMake 3.20 or later.

Plots are optional. Add the `plot` extra to install Matplotlib:

```sh
python -m pip install "obdeect-ray-tracing[plot]"
```

## From source

From a checkout, build and test all components:

```sh
cmake --preset release && cmake --build --preset release && ctest --preset release
python -m pip install -e ./src/stream -e ./src/propagation -e ./src/sources -e "./src/ray_tracing[dev]"
for c in stream propagation sources ray_tracing; do python -m pytest src/$c/python/tests; done
```

Optional CMake switches:

- `-DOBDEECT_RAY_TRACING_USE_FMA=ON` (default OFF) evaluates the aspheric mirror polynomials
  with `std::fma`. Use it only on hardware with fused multiply-add; elsewhere `std::fma` is
  emulated and slow. Results change at the rounding level.
- Results also differ in the last bits between platforms, because `std::pow`, `std::hypot`
  and the trigonometric functions differ. See
  [Validation](../ray-tracing/VALIDATION.md#rounding-level-changes).

## Use the C++ libraries in your own code

The stream component installs its headers and a CMake package; see
[Programming interfaces](../stream/programming.md#c). The other components build against the
stream target of the combined build, the stream sources in `src/stream`, the copy bundled in
their source distribution, or an installed package (`obdeect_stream_DIR`, then
`python -m obdeect.stream --cmake-dir`), in this order. A source distribution carries the
stream sources of its own version, so it builds in an isolated environment. `obdeect-stream`
of the same version remains the runtime dependency.

## Build this documentation

```sh
python -m pip install "./src/ray_tracing[docs]"
sphinx-build -W --keep-going -b html docs docs/_build/html
```
