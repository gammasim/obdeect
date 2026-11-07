# rotating-led

Port of sim_telarray `LightEmission/fpls.cc`: a focal-plane LED at 0.7 times the camera
radius above the telescope centre, tilted 15 deg from the vertical and rotated once per
rotation period. Command: `dlamb-rotating-led [options] > out.gst`.

The original has no options. Its constants are options here, with the original defaults.

## Options

| Option | Default | Original | Meaning |
| --- | --- | --- | --- |
| `--altitude M` | `2100` | `fpls.cc:67` | observation level [m] |
| `--camera-radius R` | `1.4` | `fpls.cc:68` | fiducial radius [m]; LED height 0.7 R |
| `--photons N` | `5000` | `fpls.cc:70` | photons into 4 pi per time step |
| `--bunchsize F` | `1` | `fpls.cc:71` | photons per bunch (at least 1) |
| `--rotation-period T` | `60` | `fpls.cc:72` | rotation period [ns] |
| `--time-step DT` | `0.25` | `fpls.cc:87` | time step [ns] |
| `--spectrum W` | `400` | `fpls.cc:82` (`ls_spec.lis`) | wavelength [nm]; the table file is not available |
| `--lightpulse S` | `Simple:0` | `fpls.cc:82` | as in camera-flasher |
| `--angular-distribution S` | `isotropic` | `fpls.cc:82` (`ls_beamed.lis`) | beam profile; the table file is not available |
| `--run N`, `--atmosphere N` | `1`, `1` | `fpls.cc:73, 76` | run keys |
| `--seed N` | `1` | new | Philox seed |
| `-o FILE` | stdout | | output file |

## Mapping to dlamb

- Frame: `sky`. LED at `(0, 0, 0.7 R)`.
- One emission per time step `t = 0, dt, ..., <= period + 0.001` (the loop of
  `fpls.cc:87`), each with pointing `DirectionVect(-75 deg, 0)` rotated by `2 pi t / period`
  about z (`fpls.cc:83-90`), time offset `t`, `photons` photons. 241 emissions by default.
- Target: none (the LED lies inside the fiducial sphere, as in the original).
- Bypass: meant for sim_telarray `Bypass_Optics=2`.

## Differences to the original

- Spectrum and beam come from options (`ls_spec.lis`, `ls_beamed.lis` are not in this
  checkout). The default output is therefore not the original's.
- Times are emission times (no flight time, no refraction). See camera-flasher.
- Cross-check against the original is not feasible: it needs the two table files.

## Checks

- `ctest`: `test_rotating-led` (time steps, rotation of the parallel beam, counts with the default
  beam).
- Original run in the container: not feasible (missing tables and atmosphere path).
