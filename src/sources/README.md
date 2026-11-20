# obdeect-sources

Artificial light sources for Cherenkov-telescope simulations: flashers, LEDs, point sources,
laser beams, satellite lidar, night-sky background, muon tracks and simple air showers. Sources
emit photons in the sky, telescope or camera frame. The C++20 core uses only the standard library.
It has a Python API and writes photon streams.

```sh
pip install obdeect-sources     # this component only (also installs obdeect-stream)
pip install obdeect             # all obdeect components
```

Requires Python 3.14 or higher.

## Commands

```sh
obdeect sources run flasher --help                       # presets: options by --help
obdeect sources run run.toml -o out.pst                  # run file
obdeect sources check run.toml                           # validate a run file
```

The `lidar` and `laser-beam` presets emit a sampled point pulse; use
`obdeect propagation run --mode weighted --scattering TABLE --scattering-points N`
for scattered light. Use `obdeect sources run lidar --help` for geometry and pulse options.

Full documentation and the repository README: https://github.com/gammasim/obdeect
