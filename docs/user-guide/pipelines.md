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
| sky | CORSIKA (`obdeect stream from-corsika`), `obdeect sources run` with the presets `point-source`, `beam-scan`, `laser-beam`, `lidar`, `night-sky`, `muon`, `shower` | atmosphere (`obdeect propagation run`), then optics (`obdeect ray-tracing trace-stream`) unless sim_telarray traces the optics |
| telescope (T) | `obdeect sources run flasher --type secondary` | none |
| camera (C) | `obdeect sources run flasher --type dish`, `obdeect sources run led` | none |

`obdeect stream to-simtel` prints the sim_telarray option (`Bypass_Optics`) that each telescope
needs: for camera-frame photons 1 (single mirror) or 2 (dual mirror), for telescope-frame photons
1. The examples below use 1 for single-mirror telescopes. Replace `<config>` below with your
sim_telarray configuration.

## Air showers through the full chain

```console
$ zstd -dc run.corsika.zst \
    | obdeect stream from-corsika \
    | obdeect propagation run --transmission atm.dat --site-altitude 2156 --seed 1 \
        --pathlen "$(obdeect propagation pathlen telescope.json)" \
    | obdeect ray-tracing trace-stream --model telescope.json --seed 7 \
        --pointing-zenith-deg 20 --pointing-azimuth-deg 0 \
    | obdeect stream to-simtel --camera-map cam.txt > sim.simtel
```

- `telescope.json` is a compiled [optical model](../ray-tracing/optical-models.md).
- `cam.txt` is the **camera map**: it lists the telescopes whose photons are in the camera
  frame. See [Command line](../stream/command-line.md#camera-map-and-telescope-frames).

## Dish flasher in front of the camera

The dish flasher writes camera-frame photons (`--position` in the camera frame, pointing at the
camera plane):

```console
$ obdeect sources run flasher --type dish --position 0,0,1.6 --camera-radius 0.9 \
    --observation-level 2150 --spectrum fixed:wavelength_nm=400 --photons 1e6 \
    | obdeect stream to-simtel --camera-map cam.txt \
    | sim_telarray -c <config> -C Bypass_Optics=1 -
```

## Flasher at the camera corners (dual-mirror telescopes)

The four flashers fire onto the secondary mirror. The photons are in the telescope frame, and
`obdeect stream to-simtel` needs a camera map that marks the telescope as dual-mirror
([details](../stream/command-line.md#telescope-frame-photons-dual-mirror-telescopes)).
`cam.txt` then contains the line `1 2 0.0`:

```console
$ obdeect sources run flasher --type secondary --flasher-offset 0.5 --flasher-depth 0.3 \
    --flasher-inclination 2 --mirror-camera-distance 16 --fiducial-radius 12.5 \
    --observation-level 2150 --spectrum fixed:wavelength_nm=400 --photons 1e6 \
    --delays 0 0 0 0 --fire all \
    | obdeect stream to-simtel --camera-map cam.txt \
    | sim_telarray -c <config> -C Bypass_Optics=1 -
```

## LED in front of camera pixels

The LED source needs the pixel positions and trigger groups, so it reads the sim_telarray
**camera definition** (`--camera-file cam.dat`). `obdeect stream to-simtel` does not read that
file. It needs only the **camera map** (`--camera-map cam.txt`), a list of the telescopes with
camera-frame photons:

```console
$ obdeect sources run led --camera-file cam.dat --use-trg-all --photons 500 \
    --camera-radius 0.9 --observation-level 2150 --spectrum fixed:wavelength_nm=400 \
    | obdeect stream to-simtel --camera-map cam.txt \
    | sim_telarray -c <config> -C Bypass_Optics=1 -
```

## Point source, beam, laser, lidar, night sky, muon and shower

These sources write sky-frame photons. If sim_telarray applies the extinction itself, run the
atmosphere step with `--no-extinction`, which moves the photons to the telescope and nothing
more:

```console
$ obdeect sources run point-source --observation-level 2150 --spectrum fixed:wavelength_nm=400 \
    --position 0,0,10000 --photons 1e8 --telescope-radius 12.5 \
    | obdeect propagation run --no-extinction --site-altitude 2150 --seed 1 \
    | obdeect stream to-simtel \
    | sim_telarray -c <config> -
```

The same pipeline works for the other sky-frame presets (`beam-scan`, `laser-beam`, `lidar`,
`night-sky`, `muon`, `shower`). Their options are in [presets](../sources/presets.md).
