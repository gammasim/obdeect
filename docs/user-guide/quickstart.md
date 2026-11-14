# Quick start

Three steps: create photons, look at them, move them through the atmosphere. Install obdeect
first ([Installation](installation.md)).

## 1. Create photons

A point light source 10 km above the telescope, 50 events:

```console
$ obdeect sources point-source -N 50 -o point-source.pst
```

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
It needs an atmospheric transmission table for the observation level (2150 m, the default of
`point-source`):

```console
$ obdeect propagation run --transmission atm_trans.dat --site-altitude 2150 --seed 1 \
    -i point-source.pst -o ground.pst
```

## Next

- [Pipelines](pipelines.md): combine components, and feed sim_telarray.
- [Light sources](../sources/index.md): what you can simulate besides a point source.
- [Ray tracing](../ray-tracing/index.md): trace the photons through a telescope.
