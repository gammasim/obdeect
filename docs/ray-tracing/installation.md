# Installation

Install the common requirements first: see [installation](../installation.md).

The ray-tracing package provides the `obdeect ray-tracing` command and the
`obdeect.ray_tracing` Python module:

```sh
python -m pip install "obdeect-ray-tracing[plot]"   # [plot] adds Matplotlib for plots
obdeect ray-tracing --help
```

## Build options

- `-DOBDEECT_RAY_TRACING_USE_FMA=ON` (CMake, default OFF) evaluates the aspheric polynomials
  with `std::fma`. Use it only on hardware with fused multiply-add: without it
  `std::fma` is emulated and very slow. Results differ from builds without it
  at the rounding level.
- Results also remain platform-dependent in the last bits, because the math
  library (`std::pow`, `std::hypot`, trigonometric functions) differs between
  platforms; see [Validation](VALIDATION.md#rounding-level-changes).

## First trace

Continue with the [user guide](user-guide.md) and the [command-line usage](usage.md).
