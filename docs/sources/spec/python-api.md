# Python API

Package `obdeect.sources`. Errors are `ValueError`; the message starts with the key path
(see [run file](run-file.md)).

## `obdeect.sources.config`

- `load(path) -> dict`: read a `.toml` or `.json` run file (not normalized).
- `normalize(config) -> dict`: validate and fill defaults. The input is not changed. Pure Python.

```python
from obdeect.sources import config

run = config.normalize(config.load("run.toml"))
run["run"]["seed"]  # 1 (the default, an int)
```

## `obdeect.sources.run`

- `write(config, output, ...) -> summary`: write the photon stream (`.pst`) to `output`. The
  summary has the counts of [`RunSummary`](run.md#types). The normalized run description is
  stored as the run metadata key `run_description`.
- `check(config) -> dict`: validate and compute the expected values without writing photons.

```python
from obdeect.sources import config, run

description = config.normalize(config.load("run.toml"))
print(run.check(description))
summary = run.write(description, "out.pst")
```

## `obdeect.sources.source`

Single emissions: an emission is one entry of `events[i].emissions` (with `kind`). Functions
accept sparse emissions and normalize them first. An `auto` target is rejected: it needs the
telescopes of a run.

- `normalize_emission(emission) -> dict`.
- `expected_photons(emission) -> float`: mean candidate rows.
- `expected_weight(emission) -> float`: expected sum of weights.
- `photon_count(emission, *, seed, event=0, telescope=0) -> int`.
- `emit(emission, *, seed, event=0, telescope=0) -> dict`: all rows as float32 arrays, keys
  `x y z cx cy cz t w wavelength`. Holds everything in memory.
- `iter_chunks(emission, *, seed, event=0, telescope=0)`: chunks of at most 65536 rows.
- `emit_chunk(emission, *, seed, event=0, telescope=0, first_index=0, rows=65536)`.
- Constants: `CHUNK_ROWS = 65536`, `COLUMNS`.

```python
from obdeect.sources.source import emit

point = {
    "kind": "point",
    "photons": 1e6,
    "position_m": [0.0, 0.0, 300.0],
    "angular": {"kind": "isotropic"},
    "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
    "target": {"kind": "sphere", "centre_m": [0.0, 0.0, 0.0], "radius_m": 12.5},
}
photons = emit(point, seed=3)
photons["w"].sum()  # expected photons on the sphere, approximately
```

## `obdeect.sources.plot` (optional, needs matplotlib)

- `profile_figure(emission, *, seed=1, photons=200000)`: histograms of angle, wavelength and
  emission time of the first photons (weighted by `w`). An emission with no photons raises
  `ValueError`.

## `obdeect.sources._core`

Compiled extension (nanobind) used by the functions above. Not a public interface.
