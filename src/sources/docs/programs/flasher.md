# flasher

One tool for two sim_telarray flasher programs, chosen with `--type`. Command:
`dlamb-flasher [--type dish|secondary] [options] > out.gst`.

- `--type dish` (default): port of `LightEmission/ff-1m.cc`. A flasher on the axis of a
  single-reflector telescope, in front of the camera. Frame `sky`.
- `--type secondary`: port of `LightEmission/ff-gct.cc`. Four flashers at the corners of the
  camera of a dual-mirror telescope (GCT/CHEC), firing down onto the secondary mirror.
  Frame `telescope` (T).

Options that belong to one type are rejected with the other type, e.g.
`--delays applies only to --type secondary`.

## Options

| Option | Type | Default | Original | Meaning |
| --- | --- | --- | --- | --- |
| `--type dish\|secondary` | both | `dish` | | flasher type |
| `--photons N[,...]` | both | `2.5e6` | `ff-1m.cc:166-167`, `ff-gct.cc:159-160` | photons emitted into 4 pi per flasher |
| `--bunchsize F` (`-b`) | both | `1` | `ff-1m.cc:168-169`, `ff-gct.cc:161-162` | photons per bunch (at least 1); weights are F |
| `--spectrum W` | both | `400` | `ff-1m.cc:184-185`, `ff-gct.cc:171-172` | wavelength [nm], or `Cherenkov:A:B`; no table files |
| `--lightpulse S` | both | `Simple:0` | `ff-1m.cc:186-187`, `ff-gct.cc:173-174` | `Simple:W[,D]`, `Gauss:S`, `Gauss-rms:S`, `Gauss-fwhm:F`, `Exp:T`, `Complex:A,B,C[,D]` |
| `--angular-distribution S` | both | `isotropic` | `ff-1m.cc:188-189`, `ff-gct.cc:175-176` | `isotropic`, `Parallel`, `Cone:DEG`, `FilledCone:DEG`, `Gauss:DEG`, `Gauss-fwhm:DEG`, `Rayleigh` |
| `--run N` | both | `1` | `ff-1m.cc:190-191`, `ff-gct.cc:177-178` | run number |
| `--altitude M` (`-a`) | both | `2150` dish, `2100` secondary | `ff-1m.cc:164-165`, `ff-gct.cc:157-158` | observation level [m] |
| `--atmosphere N` (`-A`) | both | `1` | `ff-1m.cc:192-195`, `ff-gct.cc:179-180` | atmospheric profile number (file names not supported) |
| `--seed N` | both | `1` | new | Philox seed |
| `-o FILE` | both | stdout | `ff-1m.cc:162-163`, `ff-gct.cc:155-156` | output file |
| `--events N[,...]` | dish | `1` | `ff-1m.cc:196-200` | events per photon level; a short list repeats its last value |
| `--frame sky\|camera` | dish | `sky` | new | `camera`: photons on the camera plane z = 0, see below |
| `--xy X,Y` (`--flasher-xy X,Y`) | dish | `0,0` | `ff-1m.cc:170-179` | flasher offset from the axis [cm] |
| `--distance D` (`-z D`) | dish | `1600` | `ff-1m.cc:180-181` | flasher distance above the telescope centre [cm] |
| `--camera-radius R` (`-r R`) | dish | `120` | `ff-1m.cc:182-183` | fiducial camera radius [cm] |
| `--events N` | secondary | `1` | `ff-gct.cc:181-182` | events (one value) |
| `--flasher-xy X` (`--xy X`) | secondary | `14.62` | `ff-gct.cc:163-164` | flasher offset on the diagonal [cm] |
| `--flasher-depth D` (`-z D`, `--z D`) | secondary | `2.46` | `ff-gct.cc:165-166` | flasher depth behind the camera front [cm] |
| `--flasher-inclination A` (`--inclination A`) | secondary | `0` | `ff-gct.cc:167-168` | outward inclination [deg] |
| `--mirror-camera-distance D` (`--mirror-distance D`) | secondary | `51.06949` | `ff-gct.cc:169-170` | focal front to secondary [cm] |
| `--delays T1,..,T4[,T5]` | secondary | `0,0,0,0,0` | `ff-gct.cc:183-185` | flasher delays and the common-flash delay [ns] |
| `--fire all\|alternate\|demo\|1..4` | secondary | `all` | `ff-gct.cc:186-208` | which flashers fire |

## Mapping to dlamb

- dish: the emitter sits at `(x, y, distance)` above the telescope centre, which is the
  origin at observation level. Emitter pointing `(0, 0, -1)`. Photon count
  `photons / bunchsize`. For isotropic emission beyond `1.02 R` the source uses a
  `SphereTarget` of radius R (cone with `sin(theta) = R / d`). Fixed-angle profiles keep all
  photons. Other profiles use a `DiscTarget` on the sphere's plane (weights; an approximation
  of the original cut).
- secondary: origin at the secondary centre, z along the optical axis. Flashers at
  `(+-xy, +-xy, mirror distance + depth)` (ls1..ls4 of `ff-gct.cc:218-229`), pointing
  `DirectionVect(-90 deg + inclination, atan2(y, x))` (`ff-gct.cc:220`). Count
  `photons / bunchsize` per flasher. The flashers lie inside the fiducial sphere (120 cm),
  so no target: 4 pi emission, as the original.
- secondary fire modes: `all` (four emissions), `alternate` (flasher `event mod 4`), `demo`
  (four individual flashes, then all with `delays[4]`), `n` (one flasher). `demo` with
  all-zero delays uses the demo values (`ff-gct.cc:197-206`).
- Bypass: dish is meant for sim_telarray `Bypass_Optics=2` (or 1 for a true one-mirror
  telescope); the atmosphere stage (G, refraction) then an optics stage with no mirror.
  secondary is meant for `Bypass_Optics=1` (primary bypassed, secondary kept); the photons
  must still hit the secondary, so the optics stage applies its reflection.

`--frame camera` (dish): frame C. The photons run in a straight line from the flasher to the
plane z = 0 of the camera, without atmosphere. Time is emission time plus path / c. Count and
directions equal the sky-frame output for the same seed; positions lie in a disc of radius
`R / sqrt(1 - (R / d)^2)`. A photon that does not travel towards the plane is an error.
Pipeline: see [](../sim_telarray.md).

## Differences to the originals

- Times are emission times (plus the delay for secondary). The originals write arrival times:
  dish at the observation plane (emission time plus `d / c`, 53.4 ns for the default
  distance), secondary at the plane z = 0 (flight time not added here). The atmosphere
  stage must add path delays.
- No refraction (applied by the atmosphere stage). No EventIO.
- Positions are emission points, not impacts on a plane.
- dish: the cut uses `sin(theta) = R / d`. The original's emitted cone matches this (count and
  cz checked, `scripts/compare_flasher-dish.sh`), although its code writes `1.02 R`.
- secondary: `--delays` with fewer than 4 values is an error; the original keeps the old values.
- The original `--xy` and `--flasher-xy` names are kept for both types with the meaning of
  the type (a pair for dish, one value for secondary).

## Checks

- `ctest`: `test_flasher` (dish: count, cone and direction moments, zem, non-default spectrum,
  pulse, bunch size, event repeats; secondary: flasher positions, direction moments, counts,
  delays, fire patterns, alternate events; type rules).
- dish cross-check (default options, 20 events; `scripts/compare_flasher-dish.sh`): photons per
  event 3519.7 (original) vs 3520.5 (dlamb); cz mean -0.998593 +- 3e-6 vs -0.998583 +- 3e-6
  (2.2 sigma); wavelength 400 nm both; original time 53.45 ns is the flight time.
- secondary cross-check (`--photons 1e5 --events 1`; `scripts/compare_flasher-secondary.sh`):
  399939 bunches (original) vs 400000 (dlamb). Original cz is always negative in eventio (sign
  lost), so the cz comparison is not valid for this isotropic emitter; cx, cy agree within 1 sigma.
