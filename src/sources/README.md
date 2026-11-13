# obdeect sources

Artificial light sources for Cherenkov-telescope simulations.

Artificial light sources for Cherenkov-telescope simulations: flashers, illuminators, lasers, LEDs and stars, emitting photons at sky, telescope-sphere, telescope or camera level, as a standard-library-only C++20 core with a Python API writing photon streams.

`obdeect sources` is one module of a modular Cherenkov-telescope simulation chain
(CORSIKA7 or light sources, atmosphere, optics, camera). Modules exchange photons as
[obdeect stream](https://github.com/gammasim/obdeect) photon streams;
[simtools](https://github.com/gammasim/simtools) configures and runs the chain.

## User Installation

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install obdeect-sources
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

BSD 3-Clause License, see [LICENSE](LICENSE). Please cite `obdeect sources` using the
metadata in [CITATION.cff](CITATION.cff).
