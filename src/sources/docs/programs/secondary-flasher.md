# secondary-flasher

Port of sim_telarray `LightEmission/ff-gct.cc`: four flashers at the corners of the camera
of a dual-mirror telescope (GCT/CHEC), firing down onto the secondary mirror. Command:
`dlamb-secondary-flasher [options] > out.gst`.

## Options

| Option | Default | Original | Meaning |
| --- | --- | --- | --- |
| `--events N` | `1` | `ff-gct.cc:181-182` | events |
| `--photons N[,...]` | `2.5e6` | `ff-gct.cc:159-160` | photons into 4 pi per flasher |
| `--bunchsize F` | `1` | `ff-gct.cc:161-162` | photons per bunch (at least 1) |
| `--flasher-xy X` (`--xy`) | `14.62` | `ff-gct.cc:163-164` | flasher offset on the diagonal [cm] |
| `--flasher-depth D` (`-z`) | `2.46` | `ff-gct.cc:165-166` | flasher depth behind the camera front [cm] |
| `--flasher-inclination A` | `0` | `ff-gct.cc:167-168` | outward inclination [deg] |
| `--mirror-camera-distance D` | `51.06949` | `ff-gct.cc:169-170` | focal front to secondary [cm] |
| `--spectrum W`, `--lightpulse S`, `--angular-distribution S` | `400`, `Simple:0`, `isotropic` | `ff-gct.cc:171-176` | as in camera-flasher |
| `--delays T1,..,T4[,T5]` | `0,0,0,0,0` | `ff-gct.cc:183-185` | flasher delays and the common-flash delay [ns] |
| `--fire all\|alternate\|demo\|1..4` | `all` | `ff-gct.cc:186-208` | which flashers fire |
| `--run N`, `--altitude M`, `--atmosphere N` | `1`, `2100`, `1` | `ff-gct.cc:177-180` | run keys |
| `--seed N` | `1` | new | Philox seed |
| `-o FILE` | stdout | `ff-gct.cc:155-156` | output file |

## Mapping to dlamb

- Frame: `telescope` (T). Origin at the secondary centre (the telescope sphere centre of
  the original), z along the optical axis towards the sky.
- Flashers: positions `(+-xy, +-xy, mirror distance + depth)` (ls1..ls4 of `ff-gct.cc:218-229`),
  pointing `DirectionVect(-90 deg + inclination, atan2(y, x))` (`ff-gct.cc:220`).
- Count: `photons / bunchsize` per flasher. The flashers lie inside the fiducial sphere
  (120 cm), so no target: 4 pi emission, as the original (`frac = 1`).
- Fire modes: `all` (four emissions), `alternate` (flasher `event mod 4`), `demo` (four
  individual flashes, then all with `delays[4]`), `n` (one flasher). `demo` with all-zero
  delays uses the demo values (`ff-gct.cc:197-206`). The original applies them while
  parsing; here they are applied in the program.
- Bypass: meant for sim_telarray `Bypass_Optics=1` (primary bypassed, secondary kept).
  The photons must still hit the secondary, so the optics stage applies its reflection.

## Differences to the original

- Times are emission times plus the delay; the original writes arrival times at the plane
  z = 0 (flight time from each flasher to that plane, not added here).
- Positions are emission points (the original writes impacts on the secondary plane and
  sim_telarray moves them back to the emission point).
- `--delays` with fewer than 4 values is an error; the original keeps the old values.
- No refraction, no EventIO.

## Checks

- `ctest`: `test_secondary-flasher` (flasher positions, direction moments, counts, delays, fire
  patterns, alternate events).
- Cross-check with `--photons 1e5 --events 1` (`scripts/compare_secondary-flasher.sh`): 399939
  bunches (original) vs 400000 (dlamb). Original cz is always negative in eventio (sign
  lost), so the cz comparison is not valid for this isotropic emitter; cx, cy agree within
  1 sigma.
