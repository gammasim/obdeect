# Command line

Compile an optical model once, then reuse it to trace and plot. Use a telescope and
version available in your simulation-models checkout.

```sh
# Compile an optical model.
obdeect ray-tracing compile --simulation-models ../simulation-models --model LSTN-01 \
  --version 7.0.0 --seed 1 --zenith-deg 20 --output lst.json

# Trace a star 0.5 deg off axis; print the image summary, keep every photon.
obdeect ray-tracing trace --model lst.json --seed 2 --photons 100000 \
  --field-deg 0.5 0 --output lst-arrivals.npz

# Draw the dish seen from the sky and the focal-plane image.
obdeect ray-tracing plot --model lst.json --arrivals lst-arrivals.npz --output lst.png

# Which obscurers (camera, baffle, masts) shade a star 2 deg off axis.
obdeect ray-tracing obscurers --model lst.json --seed 2 --photons 100000 --field-deg 2 0
```

- `compile` reads full and patch productions. A patch is laid over the productions
  in its `model_version_history`.
- `compile --seed` fixes the random mirror errors (focal lengths, alignment).
  `--zenith-deg` is the zenith angle for zenith-dependent focus offsets and alignment errors.
- `trace --seed` fixes source sampling and mirror scattering. `--geometry-only` keeps
  the photon weights unchanged (no reflectivity, telescope transmission or degradation).
- Results do not depend on `--threads`, photon order or block size.

Run `obdeect ray-tracing <command> --help` for all options. For structural shadow
diagnostics, see [Shadowing](shadowing.md).
