# Pipelines

Each obdeect command reads a [photon stream](../stream/format.md) and writes one, so you build a
simulation by piping commands together. `obdeect stream to-simtel` closes the pipe: it turns
a stream into sim_telarray input.

```text
source ──▶ atmosphere ──▶ optics ──▶ to-simtel ──▶ sim_telarray
```

Which steps you need depends on the **frame** in which the source writes its photons:

| Frame | Written by | Steps before `obdeect stream to-simtel` |
| --- | --- | --- |
| sky | CORSIKA (`obdeect stream from-corsika`), `obdeect sources point-source`, `beam-scan`, `laser-beam`, `night-sky`, `muon-track`, and `flasher --type dish` by default | atmosphere (`obdeect propagation run`), then optics (`obdeect ray-tracing trace-stream`) unless sim_telarray traces the optics |
| telescope (T) | `obdeect sources flasher --type secondary` | none |
| camera (C) | `obdeect sources led`, `obdeect sources flasher --type dish --frame camera` | none |

`obdeect stream to-simtel` prints the sim_telarray option (`Bypass_Optics`) that each telescope
needs. Replace `<config>` below with your sim_telarray configuration.

## Air showers through the full chain

```console
$ zstd -dc run.corsika.zst \
    | obdeect stream from-corsika \
    | obdeect propagation run --transmission atm.dat --site-altitude 2156 --seed 1 \
    | obdeect ray-tracing trace-stream --model telescope.json --seed 7 \
        --pointing-zenith-deg 20 --pointing-azimuth-deg 0 \
    | obdeect stream to-simtel --camera-map cam.txt > sim.simtel
```

- `telescope.json` is a compiled [optical model](../ray-tracing/optical-models.md).
- `cam.txt` is the **camera map**: it lists the telescopes whose photons are in the camera
  frame. See [Command line](../stream/command-line.md#camera-map-and-telescope-frames).

## Flasher on the dish

```console
$ obdeect sources flasher --type dish --photons 1e6 \
    | obdeect propagation run --transmission atm.dat --site-altitude 2156 --seed 1 \
    | obdeect ray-tracing trace-stream --model telescope.json --seed 7 \
    | obdeect stream to-simtel --camera-map cam.txt > sim.simtel
```

With `--frame camera`, the flasher sends photons in a straight line to the camera plane
(no atmosphere, arrival time = emission time + path / c) and sim_telarray takes over from there:

```console
$ obdeect sources flasher --type dish --frame camera --distance 1600 --photons 2.5e6 \
    | obdeect stream to-simtel --camera-map cam.txt \
    | sim_telarray -c <config> -C Bypass_Optics=3 -
```

## Flasher at the camera corners (dual-mirror telescopes)

The four flashers fire onto the secondary mirror. The photons are in the telescope frame, and
`obdeect stream to-simtel` needs a camera map that marks the telescope as dual-mirror
([details](../stream/command-line.md#telescope-frame-photons-dual-mirror-telescopes)):

```console
$ obdeect sources flasher --type secondary --fire all --delays 0,0,0,0 \
    | obdeect stream to-simtel --camera-map cam.txt \
    | sim_telarray -c <config> -C Bypass_Optics=1 -
```

## LED in front of camera pixels

The LED source needs the pixel positions and trigger groups, so it reads the sim_telarray
**camera definition** (`--camera-file cam.dat`). `obdeect stream to-simtel` does not read that
file. It needs only the **camera map** (`--camera-map cam.txt`), a three-column list of the
telescopes with camera-frame photons:

```console
$ obdeect sources led --camera-file cam.dat --use-trg-all --photons 500 \
    | obdeect stream to-simtel --camera-map cam.txt \
    | sim_telarray -c <config> -C Bypass_Optics=3 -
```

## Point source, beam, laser, night sky and muon

These sources write sky-frame photons. If sim_telarray applies the extinction itself, run the
atmosphere step with `--no-extinction`, which moves the photons to the telescope and nothing
more:

```console
$ obdeect sources point-source -z 10000000 -n 1e7 \
    | obdeect propagation run --no-extinction --site-altitude 2150 --seed 1 \
    | obdeect stream to-simtel \
    | sim_telarray -c <config> -
```

The same pipeline works for `beam-scan`, `laser-beam`, `night-sky` and `muon-track`.
