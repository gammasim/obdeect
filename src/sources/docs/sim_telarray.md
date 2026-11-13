# Using dlamb with sim_telarray

dlamb writes gstrickd photon streams. `obdeect stream to-simtel` (gstrickd) converts a stream into
sim_telarray input, so dlamb replaces the light-emission package of sim_telarray. The frame of
the stream decides which steps lie between dlamb and `obdeect stream to-simtel`.

| Frame | Written by | Between dlamb and `obdeect stream to-simtel` |
| --- | --- | --- |
| camera (C) | `obdeect sources led`, `obdeect sources flasher --type dish --frame camera` | nothing; `--camera-map` lists the telescope |
| telescope (T) | `obdeect sources flasher --type secondary` | nothing; needs the telescope-frame input of `obdeect stream to-simtel` |
| sky | `obdeect sources flasher --type dish` (default), `obdeect sources point-source`, `obdeect sources beam-scan`, `obdeect sources laser-beam`, `obdeect sources night-sky`, `obdeect sources muon-track` | `obdeect propagation run` |

`obdeect stream to-simtel` prints the sim_telarray option that each telescope needs.

## LED in front of camera pixels (camera frame)

```console
$ obdeect sources led --camera-file cam.dat --use-trg-all --photons 500 \
    | obdeect stream to-simtel --camera-map cam.txt | sim_telarray -c cfg -C Bypass_Optics=3 -
```

## Flasher on the dish axis (camera frame)

`--frame camera` sends the photons in a straight line from the flasher to the camera plane
z = 0. There is no atmosphere. The arrival time is the emission time plus path / c.

```console
$ obdeect sources flasher --type dish --frame camera --distance 1600 --photons 2.5e6 \
    | obdeect stream to-simtel --camera-map cam.txt | sim_telarray -c cfg -C Bypass_Optics=3 -
```

Without `--frame camera` the photons stay at the flasher (sky frame) and need
`obdeect propagation run` as for the other sky-frame tools.

## Flasher at the camera corners, secondary mirror (telescope frame)

```console
$ obdeect sources flasher --type secondary --fire all --delays 0,0,0,0 \
    | obdeect stream to-simtel | sim_telarray -c cfg -
```

This needs the telescope-frame input of `obdeect stream to-simtel`, which is added to gstrickd separately.

## Point source, beam, laser, night sky, muon (sky frame)

`--no-extinction` propagates the geometry only; the extinction is left to sim_telarray.

```console
$ obdeect sources point-source -z 10000000 -n 1e7 \
    | obdeect propagation run --no-extinction ... | obdeect stream to-simtel | sim_telarray -c cfg -
```

The same pipeline works for `obdeect sources beam-scan`, `obdeect sources laser-beam`, `obdeect sources night-sky` and
`obdeect sources muon-track`.
