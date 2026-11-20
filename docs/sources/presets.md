# Presets

A preset builds a run description from a few options, so common setups need no run file. Use it
in place of the run file: `obdeect sources run <preset> ...`. Every preset needs
`--observation-level` (m above sea level). `obdeect sources run <preset> --help` lists all
options.

The examples write to `/dev/null`; replace it with `-o out.pst` to keep the photons. Numbers
are small for a quick test. Profile options take `KIND:NAME=VALUE,...`, for example
`--spectrum fixed:wavelength_nm=400`.

## flasher

Camera flasher, with `--type`:

- `dish` (default): one point emission at `--position` in the camera frame, pointing down into
  the camera plane. Frame `camera`, target `auto` with the sphere of `--camera-radius`.
- `secondary`: four point emissions in the telescope frame, one per flasher. Flasher `i` sits at
  `(+-offset, +-offset, mirror-camera-distance + depth)` and points at `-90 deg + inclination`.
  Frame `telescope`, target `auto` with `--fiducial-radius`. `--fire` picks the flashers of an
  event and `--delays` sets their time offsets.

Each value of `--photons` is one event. Required: `--spectrum`, `--photons`; for `dish`
`--position` and `--camera-radius`; for `secondary` the `--flasher-*`, `--mirror-camera-distance`
and `--fiducial-radius` options.

```sh
obdeect sources run flasher --type dish --position 0,0,1.6 --camera-radius 0.9 \
    --observation-level 2150 --spectrum fixed:wavelength_nm=400 --photons 1e4 -o /dev/null
```

## led

LEDs in front of selected camera pixels. Frame `camera`. Each lit pixel is one point emission,
parallel along -z into the camera plane, with target `auto` (sphere of `--camera-radius`).
The pixels come from `--camera-file`, a sim_telarray camera definition with `Pixel` and
`Trigger` lines. Choose the pixels with one of `--use LIST`, `--use-trg-rand N` or
`--use-trg-all`.

Required: `--camera-file`, one pixel choice, `--photons`, `--camera-radius`, `--spectrum`.

```sh
obdeect sources run led --camera-file cam.dat --use-trg-all --photons 500 \
    --camera-radius 0.9 --observation-level 2150 --spectrum fixed:wavelength_nm=400 -o /dev/null
```

## point-source

A point light source at `--position` (sky frame, z is the height above the observation level),
emitting into 4 pi along `--direction` (default down). One event per `--photons` value. Frame
`sky`, target `auto`, so only photons that reach the telescope are written.

Required: `--position`, `--photons`, `--spectrum`, `--telescope-radius` (one telescope at the
origin).

```sh
obdeect sources run point-source --observation-level 2150 --spectrum fixed:wavelength_nm=400 \
    --position 0,0,10000 --photons 1e8 --events 50 --telescope-radius 12.5 -o /dev/null
```

## beam-scan

A parallel downward beam (cone of 0 deg) at each position of a scan over the camera area. The
scan starts at `--scan-start`, runs along x for `--scan-extent-x`, then along y for
`--scan-extent-y`, in steps of `--scan-step`, with the beam at height `--height`. Each position
is one point emission of one event. Frame `sky`, target `none`: all photons are written.

Required: `--photons-per-position`, `--spectrum`, `--camera-radius`, `--scan-start X Y`,
`--scan-extent-x`, `--scan-extent-y`, `--scan-step`, `--height`.

```sh
obdeect sources run beam-scan --observation-level 2150 --spectrum fixed:wavelength_nm=400 \
    --photons-per-position 1e4 --camera-radius 1.5 --scan-start -1 -1 \
    --scan-extent-x 2 --scan-extent-y 2 --scan-step 0.5 --height 1000 -o /dev/null
```

## laser-beam

A pulsed laser pointing down, at `--laser-position` (sky frame, z above the observation level).
One point emission per event: the pulse is `E lambda / (h c)` primary photons, with pulse energy
`--pulse-energy` (J) and `--wavelength-nm` (exact SI constants; the spectrum is fixed at that
wavelength). `--divergence` is the full beam angle. Frame `sky`, target `none`. Scattering is not
generated here; propagation computes it (see the chain below).

Required: `--laser-position X,Y,Z`, `--beam-zenith` (> 90), `--beam-azimuth`, `--divergence`,
`--pulse-energy`, `--wavelength-nm`, one of `--telescope-radius` or `--telescope X,Y,Z,R`.

```sh
obdeect sources run laser-beam --observation-level 2150 --laser-position 200,0,3000 \
    --beam-zenith 170 --beam-azimuth 180 --divergence 0.05 --pulse-energy 1e-3 \
    --wavelength-nm 355 --telescope-radius 12.5 --rows 1000 -o /dev/null
```

## lidar

A pulsed satellite lidar pointing down. Same model as `laser-beam`, with the emission at the
satellite: `--satellite-altitude M` above the observation level at `--satellite-offset X,Y`
(default 0,0), or `--satellite-position X,Y,Z`. Each pulse is one event (`--events`).

Required: `--satellite-altitude` or `--satellite-position`, `--beam-zenith` (> 90),
`--beam-azimuth`, `--divergence`, `--pulse-energy`, `--wavelength-nm`, one of `--telescope-radius`
or `--telescope X,Y,Z,R`. `--start-height M` moves the emission along the beam to that height
(for example 100000 for a satellite); `--rows` sets the sampled rows per pulse.

```sh
obdeect sources run lidar --observation-level 2150 --satellite-altitude 500000 \
    --satellite-offset 1000,0 --beam-zenith 180 --beam-azimuth 0 --divergence 0.01 \
    --pulse-energy 1e-3 --wavelength-nm 532 --telescope-radius 12.5 --start-height 100000 \
    --rows 1000 -o /dev/null
```

Chain of `laser-beam` and `lidar` (the scattered light, with weighted mode because the primary
weights exceed 2048):

```sh
obdeect sources run lidar ... -o primaries.pst
obdeect propagation run --transmission T --site-altitude H --seed N --mode weighted \
    --scattering S --scattering-points 64 -i primaries.pst -o scattered.pst
```

Model and limits, shared by both presets:

- t = 0 is the laser firing. With `--start-height`, `time_offset_ns` is the path from the laser to
  the start divided by c (vacuum), so the times stay relative to the firing.
- The beam must point down (`--beam-zenith` > 90) and tilt at most `--max-tilt` (default 60 deg)
  from nadir. The scattering model is plane-parallel.
- Only single scattering. A beam that passes a telescope sphere (closest approach within the
  sphere plus the beam radius at that point, below the start) is an error, unless
  `--allow-direct-beam` is given.
- `--beam-azimuth` is astronomical: clockwise from north (sky frame x north, y west).

## night-sky

Night-sky background: one `sky_background` emission per event on the sphere of one telescope,
with the axis pointing down (the telescope looks up). Photons arrive from within
`--viewcone-deg` of the axis. `--photons` is per event. Frame `sky`, target `none`.
`--telescope-zenith` and `--telescope-azimuth` set the pointing (default: straight up).

Required: `--spectrum`, `--photons`, `--viewcone-deg`, `--telescope-position`,
`--telescope-radius`.

```sh
obdeect sources run night-sky --observation-level 2150 \
    --spectrum uniform:low_nm=300,high_nm=600 --photons 1e5 --viewcone-deg 60 \
    --telescope-position 0,0,0 --telescope-radius 12.5 -o /dev/null
```

## muon

Cherenkov light of a muon-like track on the telescope sphere. The track runs from
`--start-height` down to the impact point (`--impact-x`, `--impact-y`) in the direction of
`--zenith` and `--azimuth`. Its light is one line emission, with the cone at the Cherenkov
angle from the refractive index, and target `auto` (the telescope sphere). Frame `sky`. The
time is zero at the impact point. Each core position of `--core-scatter-count` is one event
entry; `--events` repeats each one.

The index is required: one of `--n-minus-one ETA`, `--index-exponential ETA0 SCALE_HEIGHT_M` or
`--index-profile FILE`. The Cherenkov band is `--wavelength-min` to `--wavelength-max` (nm).
`--primary` (default muon) sets mass and charge; `--momentum` or `--beta` sets the speed.

Required: `--observation-level`, `--wavelength-min`, `--wavelength-max`, an index option,
`--telescope-radius`, `--start-height`.

```sh
obdeect sources run muon --observation-level 2150 --wavelength-min 300 --wavelength-max 600 \
    --n-minus-one 2.9e-4 --telescope-radius 12.5 --start-height 10000 --zenith 20 \
    --azimuth 0 --events 10 -o /dev/null
```

## shower

A simple model of a straight air shower: Cherenkov and fluorescence light. The axis runs from
`--start-height` to the core (`--core-x`, `--core-y`) in the direction of `--zenith` and
`--azimuth`. It is split into steps of `--step` metres. Per step, one Cherenkov line (cone at the
Cherenkov angle, target `auto`) and one isotropic fluorescence line (target `auto`). The number
of particles follows a Gaisser-Hillas profile (`--nmax`, `--xmax`, `--x0`, `--gh-lambda`) in
slant depth; the atmosphere is `--atmosphere-depth0` with `--scale-height`.

The numbers are inputs, not calibrated values. The model has no lateral distribution, no
electron spread and no multiple scattering. Rows grow with `--nmax` times `--step`, so keep
both small for a test.

Required: `--observation-level`, `--wavelength-min`, `--wavelength-max`, an index option,
`--telescope-radius`, `--start-height`, `--zenith`, `--azimuth`, `--step`, `--nmax`, `--xmax`,
`--x0`, `--gh-lambda`, `--atmosphere-depth0`, `--scale-height`, `--fluorescence-yield`,
`--fluorescence-spectrum`.

```sh
obdeect sources run shower --observation-level 2150 --wavelength-min 300 --wavelength-max 600 \
    --n-minus-one 2.9e-4 --telescope-radius 12.5 --start-height 10000 --zenith 20 --azimuth 0 \
    --step 100 --nmax 100 --xmax 600 --x0 -50 --gh-lambda 70 --atmosphere-depth0 800 \
    --scale-height 8000 --fluorescence-yield 4 --fluorescence-spectrum fixed:wavelength_nm=337 \
    -o /dev/null
```
