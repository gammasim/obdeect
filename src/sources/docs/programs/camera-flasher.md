# camera-flasher

Port of sim_telarray `LightEmission/ff-1m.cc`: a flasher on the axis of a single-reflector
telescope, in front of the camera. Command: `dlamb-camera-flasher [options] > out.gst`.

## Options

| Option | Default | Original | Meaning |
| --- | --- | --- | --- |
| `--events N[,...]` | `1` | `ff-1m.cc:196-200` | events per photon level; a short list repeats its last value |
| `--photons N[,...]` | `2.5e6` | `ff-1m.cc:166-167` | photons emitted into 4 pi per flasher |
| `--bunchsize F` (`-b`) | `1` | `ff-1m.cc:168-169` | photons per bunch (at least 1); weights are F |
| `--xy X,Y` | `0,0` | `ff-1m.cc:170-179` | flasher offset from the axis [cm] |
| `--distance D` (`-z`) | `1600` | `ff-1m.cc:180-181` | flasher distance above the telescope centre [cm] |
| `--camera-radius R` (`-r`) | `120` | `ff-1m.cc:182-183` | fiducial camera radius [cm] |
| `--spectrum W` | `400` | `ff-1m.cc:184-185` | wavelength [nm], or `Cherenkov:A:B`; no table files |
| `--lightpulse S` | `Simple:0` | `ff-1m.cc:186-187` | `Simple:W[,D]`, `Gauss:S`, `Gauss-rms:S`, `Gauss-fwhm:F`, `Exp:T`, `Complex:A,B,C[,D]` |
| `--angular-distribution S` | `isotropic` | `ff-1m.cc:188-189` | `isotropic`, `Parallel`, `Cone:DEG`, `FilledCone:DEG`, `Gauss:DEG`, `Gauss-fwhm:DEG`, `Rayleigh` |
| `--run N` | `1` | `ff-1m.cc:190-191` | run number |
| `--altitude M` (`-a`) | `2150` | `ff-1m.cc:164-165` | observation level [m] (the original `--help` says 2100) |
| `--atmosphere N` (`-A`) | `1` | `ff-1m.cc:192-195` | atmospheric profile number (file names not supported) |
| `--seed N` | `1` | new | Philox seed |
| `-o FILE` | stdout | `ff-1m.cc:162-163` | output file |

## Mapping to dlamb

- Frame: `sky`. The emitter sits at `(x, y, distance)` above the telescope centre, which is
  the origin at observation level.
- Source: `Emitter` at that position, pointing `(0, 0, -1)`, with the angular, spectral and
  pulse profiles of the options.
- Photon count: `photons / bunchsize`. For isotropic emission beyond
  `1.02 R` the source uses a `SphereTarget` of radius R (cone with `sin(theta) = R / d`).
  Fixed-angle profiles keep all photons. Other profiles use a `DiscTarget` on the sphere's
  plane (weights; an approximation of the original cut).
- Bypass: meant for sim_telarray `Bypass_Optics=2` (or 1 for a true one-mirror telescope).
  The original propagates the photons to the observation plane, which sim_telarray then
  uses as the camera plane. In the chain: atmosphere stage (G, refraction), then an optics
  stage with no mirror.

## Differences to the original

- Times are emission times. The original writes arrival times at the observation plane,
  i.e. emission time plus `d / c` (53.4 ns for the default distance). The atmosphere stage
  must add the path delay.
- No refraction (applied by the atmosphere stage).
- Positions are emission points, not observation-plane impacts.
- The cut uses `sin(theta) = R / d`: the original's emitted cone matches this (count and
  cz checked, `scripts/compare_camera-flasher.sh`), although its code writes `1.02 R`.
- Output is a gstrickd stream with `zem` (`observation level + z`); no EventIO.

## Checks

- `ctest`: `test_camera-flasher` (count, cone and direction moments, zem, non-default spectrum, pulse,
  bunch size, event repeats).
- Cross-check (default options, 20 events; `scripts/compare_camera-flasher.sh`): photons per event
  3519.7 (original) vs 3520.5 (dlamb); cz mean -0.998593 +- 3e-6 vs -0.998583 +- 3e-6
  (2.2 sigma); wavelength 400 nm both; original time 53.45 ns is the flight time.
