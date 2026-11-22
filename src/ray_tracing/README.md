# obdeect-ray-tracing

Optical ray tracer for imaging atmospheric Cherenkov telescopes. The C++20 core uses only the
standard library. The Python API compiles optical models from structured simulation-model
records or a simulation-models checkout. The native loader validates structured inputs directly
and reads saved JSON models. Status: pre-alpha.

```sh
pip install obdeect-ray-tracing  # this component only (also installs obdeect-stream)
pip install obdeect              # all obdeect components
```

Requires Python 3.14 or higher.

## Commands

```sh
obdeect ray-tracing compile --simulation-models ../simulation-models --model MSTN-01 \
  --version 7.0.0 --seed 1 --output mst.json
obdeect ray-tracing trace --model mst.json --seed 2 --photons 10000 --field-deg 0.5 0
obdeect ray-tracing trace-stream --model mst.json --seed 7 -i ground.pst -o camera.pst
obdeect ray-tracing obscurers ...
obdeect ray-tracing plot ...
```

Full documentation and the repository README: https://github.com/gammasim/obdeect
