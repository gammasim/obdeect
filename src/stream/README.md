# obdeect-stream

Photon-stream format and I/O for Cherenkov-telescope simulation chains. A versioned, columnar
photon-stream format with a C++20 reader and writer (standard library only), a Python API, and
adapters from CORSIKA IACT output and to sim_telarray input. Modules of the chain exchange
photons as photon streams.

```sh
pip install obdeect-stream      # this component only
pip install obdeect             # all obdeect components
```

Requires Python 3.14 or higher.

## Commands

```sh
obdeect stream cat [FILE ...]               # concatenate streams; --summary for record counts
obdeect stream tee FILE                     # copy a stream passing through
obdeect stream from-corsika < run.iact      # CORSIKA IACT EventIO to photon stream
obdeect stream to-simtel --camera-map cam.txt run.pst > sim.simtel   # photon stream to sim_telarray input
obdeect stream to-grisu --telescope-config TLLOC.dat run.pst > run.grisu  # GrISU photon text
```

Full documentation and the repository README: https://github.com/gammasim/obdeect
