# obdeect-ray-tracing

Optical ray tracer for imaging atmospheric Cherenkov telescopes. The C++20 core uses only the
standard library. simtools compiles optical models from simulation models; this package loads
those models, traces photons, and plots results. Status: pre-alpha.

```sh
pip install obdeect-ray-tracing  # this component only (also installs obdeect-stream)
pip install obdeect              # all obdeect components
```

Requires Python 3.14 or higher.

## Commands

```sh
simtools-compile-obdeect-optical-model --model_version 7.0.0 --site north_site \
  --telescope MSTN-01 --instrument_seed 1 --output_file mst.json
obdeect ray-tracing trace --model mst.json --seed 2 --photons 10000 --field-deg 0.5 0
obdeect ray-tracing trace-stream --model mst.json --seed 7 -i ground.pst -o camera.pst
obdeect ray-tracing obscurers ...
obdeect ray-tracing plot ...
```

Full documentation and the repository README: https://github.com/gammasim/obdeect
