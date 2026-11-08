# Using dlamb with sim_telarray

dlamb writes gstrickd photon streams. `ps2simtel` (gstrickd) converts a stream into
sim_telarray input, so dlamb replaces the light-emission package of sim_telarray. The frame of
the stream decides which steps lie between dlamb and `ps2simtel`.

| Frame | Written by | Between dlamb and `ps2simtel` |
| --- | --- | --- |
| camera (C) | `dlamb-led`, `dlamb-flasher --type dish --frame camera` | nothing; `--camera-map` lists the telescope |
| telescope (T) | `dlamb-flasher --type secondary` | nothing; needs the telescope-frame input of `ps2simtel` |
| sky | `dlamb-flasher --type dish` (default), `dlamb-point-source`, `dlamb-beam-scan`, `dlamb-laser-beam`, `dlamb-night-sky`, `dlamb-muon-track` | `dlufd-propagate` |

`ps2simtel` prints the sim_telarray option that each telescope needs.

## LED in front of camera pixels (camera frame)

```console
$ dlamb-led --camera-file cam.dat --use-trg-all --photons 500 \
    | ps2simtel --camera-map cam.txt | sim_telarray -c cfg -C Bypass_Optics=3 -
```

## Flasher on the dish axis (camera frame)

`--frame camera` sends the photons in a straight line from the flasher to the camera plane
z = 0. There is no atmosphere. The arrival time is the emission time plus path / c.

```console
$ dlamb-flasher --type dish --frame camera --distance 1600 --photons 2.5e6 \
    | ps2simtel --camera-map cam.txt | sim_telarray -c cfg -C Bypass_Optics=3 -
```

Without `--frame camera` the photons stay at the flasher (sky frame) and need
`dlufd-propagate` as for the other sky-frame tools.

## Flasher at the camera corners, secondary mirror (telescope frame)

```console
$ dlamb-flasher --type secondary --fire all --delays 0,0,0,0 \
    | ps2simtel | sim_telarray -c cfg -
```

This needs the telescope-frame input of `ps2simtel`, which is added to gstrickd separately.

## Point source, beam, laser, night sky, muon (sky frame)

`--no-extinction` propagates the geometry only; the extinction is left to sim_telarray.

```console
$ dlamb-point-source -z 10000000 -n 1e7 \
    | dlufd-propagate --no-extinction ... | ps2simtel | sim_telarray -c cfg -
```

The same pipeline works for `dlamb-beam-scan`, `dlamb-laser-beam`, `dlamb-night-sky` and
`dlamb-muon-track`.
