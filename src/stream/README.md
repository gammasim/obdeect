# obdeect stream

Photon-stream format and I/O for Cherenkov-telescope simulation chains.

A versioned, columnar photon-stream format with a standard-library-only C++20 reader/writer, a Python API, and adapters from CORSIKA IACT output and to sim_telarray input, connecting the modules of a Cherenkov-telescope simulation chain.

`obdeect stream` is one module of a modular Cherenkov-telescope simulation chain
(CORSIKA7 or light sources, atmosphere, optics, camera). Modules exchange photons as
[obdeect stream](https://github.com/gammasim/obdeect) photon streams;
[simtools](https://github.com/gammasim/simtools) configures and runs the chain.

## User Installation

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install obdeect-stream
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

BSD 3-Clause License, see [LICENSE](LICENSE). Please cite `obdeect stream` using the
metadata in [CITATION.cff](CITATION.cff).
