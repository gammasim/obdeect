# Command line

Use simtools to compile an optical model once, then reuse it to trace and plot.

```sh
# Compile an optical model with the selected site and telescope.
simtools-compile-obdeect-optical-model --model_version 7.0.0 --site north_site \
  --telescope LSTN-01 --instrument_seed 1 --output_file lst.json

# Trace a star 0.5 deg off axis; print the image summary, keep every photon.
obdeect ray-tracing trace --model lst.json --seed 2 --photons 100000 \
  --field-deg 0.5 0 --output lst-arrivals.npz

# Draw the dish seen from the sky and the focal-plane image.
obdeect ray-tracing plot --model lst.json --arrivals lst-arrivals.npz --output lst.png

# Which obscurers (camera, baffle, masts) shade a star 2 deg off axis.
obdeect ray-tracing obscurers --model lst.json --seed 2 --photons 100000 --field-deg 2 0
```

- `trace --seed` fixes source sampling and mirror scattering. `--geometry-only` keeps
  the photon weights unchanged (no reflectivity, telescope transmission or degradation).
- Results do not depend on `--threads`, photon order or block size.

Run `obdeect ray-tracing <command> --help` for all options. For structural shadow
diagnostics, see [Shadowing](shadowing.md).
