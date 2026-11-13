# Installation

Install the [common requirements](../installation.md) first.

```sh
python -m pip install "obdeect-ray-tracing[plot]"   # [plot] adds Matplotlib
obdeect ray-tracing --help
```

This provides the `obdeect ray-tracing` command and the `obdeect.ray_tracing` module.

## Build options

- `-DOBDEECT_RAY_TRACING_USE_FMA=ON` (CMake, default OFF) evaluates the aspheric
  polynomials with `std::fma`. Use it only on hardware with fused multiply-add;
  otherwise `std::fma` is emulated and very slow. Results differ at the rounding level.
- Results also differ in the last bits between platforms, because `std::pow`,
  `std::hypot` and trigonometric functions differ. See
  [Validation](VALIDATION.md#rounding-level-changes).

Next: the [user guide](user-guide.md).
