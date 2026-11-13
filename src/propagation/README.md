# dlufd

Atmospheric propagation of Cherenkov and artificial-source photons.

Atmospheric propagation for Cherenkov-telescope simulations: extinction, refraction and timing from the refractive-index profile, optional single scattering, and early photon thinning, as a standard-library-only C++20 core with a Python API, working on gstrickd photon streams.

`dlufd` is one module of a modular Cherenkov-telescope simulation chain
(CORSIKA7 or light sources, atmosphere, optics, camera). Modules exchange photons as
[gstrickd](https://github.com/gernotmaier/gstrickd) photon streams;
[simtools](https://github.com/gammasim/simtools) configures and runs the chain.

## User Installation

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install dlufd
```

(requires python 3.14 or higher)

## Developer Installation

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
cmake --preset debug
cmake --build --preset debug
```

Test your installation with:

```sh
ctest --test-dir build/debug --output-on-failure
python -m pytest python/tests
```

## License and Citation

BSD 3-Clause License, see [LICENSE](LICENSE). Please cite `dlufd` using the
metadata in [CITATION.cff](CITATION.cff).
