# obdeect-propagation

Atmospheric propagation for Cherenkov-telescope simulations: extinction, refraction and timing
from the refractive-index profile, optional single scattering, and early photon thinning. The
C++20 core uses only the standard library. It has a Python API and works on photon streams.

```sh
pip install obdeect-propagation  # this component only (also installs obdeect-stream)
pip install obdeect              # all obdeect components
```

Requires Python 3.14 or higher.

## Commands

```sh
obdeect propagation run --transmission atm.dat --site-altitude 2156 --seed 1 -i in.ps -o out.ps
obdeect propagation convert-table ...       # convert an atmospheric transmission table
```

Full documentation and the repository README: https://github.com/gammasim/obdeect
