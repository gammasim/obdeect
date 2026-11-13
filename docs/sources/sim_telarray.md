# Using obdeect sources with sim_telarray

`obdeect stream to-simtel` ([obdeect stream](../stream/usage.md)) converts a photon stream
into sim_telarray input. obdeect sources replaces the sim_telarray light-emission package.
The stream frame decides which steps lie between the two.

| Frame | Written by | Between obdeect sources and `obdeect stream to-simtel` |
| --- | --- | --- |
| camera (C) | `obdeect sources led`, `obdeect sources flasher --type dish --frame camera` | nothing; `--camera-map` lists the telescope |
| telescope (T) | `obdeect sources flasher --type secondary` | nothing; needs the telescope-frame input of `obdeect stream to-simtel` |
| sky | `obdeect sources flasher --type dish` (default), `point-source`, `beam-scan`, `laser-beam`, `night-sky`, `muon-track` | [`obdeect propagation run`](../propagation/index.md) |

`obdeect stream to-simtel` prints the sim_telarray option that each telescope needs.

## LED in front of camera pixels (camera frame)

```console
$ obdeect sources led --camera-file cam.dat --use-trg-all --photons 500 \
    | obdeect stream to-simtel --camera-map cam.txt | sim_telarray -c cfg -C Bypass_Optics=3 -
```

## Flasher on the dish axis (camera frame)

`--frame camera` sends photons in a straight line from the flasher to the camera plane z = 0.
No atmosphere. Arrival time = emission time + path / c.

```console
$ obdeect sources flasher --type dish --frame camera --distance 1600 --photons 2.5e6 \
    | obdeect stream to-simtel --camera-map cam.txt | sim_telarray -c cfg -C Bypass_Optics=3 -
```

Without `--frame camera`, photons stay at the flasher (sky frame) and need
`obdeect propagation run`, as for the other sky-frame tools.

## Flasher at the camera corners, secondary mirror (telescope frame)

```console
$ obdeect sources flasher --type secondary --fire all --delays 0,0,0,0 \
    | obdeect stream to-simtel | sim_telarray -c cfg -
```

This needs the telescope-frame input of `obdeect stream to-simtel`, added to obdeect stream
separately.

## Point source, beam, laser, night sky, muon (sky frame)

`--no-extinction` propagates the geometry only; extinction is left to sim_telarray.

```console
$ obdeect sources point-source -z 10000000 -n 1e7 \
    | obdeect propagation run --no-extinction ... | obdeect stream to-simtel | sim_telarray -c cfg -
```

The same pipeline works for `beam-scan`, `laser-beam`, `night-sky` and `muon-track`.
