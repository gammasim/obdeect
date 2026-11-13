# obdeect-sources

Artificial light sources for Cherenkov-telescope simulations: flashers, illuminators, lasers,
LEDs and stars. Sources emit photons at sky, telescope-sphere, telescope or camera level. The
C++20 core uses only the standard library. It has a Python API and writes photon streams.

```sh
pip install obdeect-sources     # this component only (also installs obdeect-stream)
pip install obdeect             # all obdeect components
```

Requires Python 3.14 or higher.

## Commands

```sh
obdeect sources flasher --type dish --photons 1e6     # also --type secondary, --frame camera
obdeect sources led --camera-file cam.dat --photons 500
obdeect sources point-source
obdeect sources beam-scan
obdeect sources laser-beam
obdeect sources night-sky
obdeect sources muon-track
```

Full documentation and the repository README: https://github.com/gammasim/obdeect
