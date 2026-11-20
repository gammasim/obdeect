# Quick start

Three steps: create photons, look at them, move them through the atmosphere. Install obdeect
first ([Installation](installation.md)).

## 1. Create photons

A point light source 10 km above the telescope, 50 events:

```console
$ obdeect sources run point-source --observation-level 2150 --spectrum fixed:wavelength_nm=400 \
    --position 0,0,10000 --photons 1e8 --events 50 --telescope-radius 12.5 -o point-source.pst
events 50, rows 1952, sum of weights 1952, expected 1953.13
```

`--photons` is the mean number of photons into 4 pi per event. Only the photons that reach the
telescope are written (about 39 per event here), each with its weight. The last line is the
summary on standard error.

The result is a [photon stream](../stream/format.md): a binary file of runs, events,
telescopes and photons. The extension `.pst` is only a convention. Every obdeect command also
reads and writes streams on standard input and output, so commands combine with pipes.

## 2. Look inside

```console
$ obdeect stream cat --summary point-source.pst
```

The summary lists the records per type and the photons and weight per telescope.

## 3. Propagate through the atmosphere

`obdeect propagation run` moves the photons to the telescope and applies timing and extinction.
It needs an atmospheric transmission table for the observation level (2150 m, the observation
level of the example above):

```console
$ obdeect propagation run --transmission atm_trans.dat --site-altitude 2150 --seed 1 \
    -i point-source.pst -o ground.pst
```

## Next

- [Pipelines](pipelines.md): combine components, and feed sim_telarray.
- [Light sources](../sources/index.md): what you can simulate besides a point source.
- [Ray tracing](../ray-tracing/index.md): trace the photons through a telescope.
