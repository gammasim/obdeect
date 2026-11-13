# flasher

One tool for two flasher types, chosen with `--type`. Replaces the earlier ff-1m and ff-gct
programs; options keep their names. Command: `obdeect sources flasher [--type dish|secondary]
[options] > out.gst`.

- `--type dish` (default): a flasher on the axis of a single-reflector telescope, in front of
  the camera. Frame `sky`.
- `--type secondary`: four flashers at the corners of the camera of a dual-mirror telescope,
  firing down onto the secondary mirror. Frame `telescope` (T).

Options that belong to one type are rejected with the other type, e.g.
`--delays applies only to --type secondary`.

## Options

| Option | Type | Default | Meaning |
| --- | --- | --- | --- |
| `--type dish\|secondary` | both | `dish` | flasher type |
| `--photons N[,...]` | both | `2.5e6` | photons emitted into 4 pi per flasher |
| `--bunchsize F` (`-b`) | both | `1` | photons per bunch (at least 1); weights are F |
| `--spectrum W` | both | `400` | wavelength [nm], or `Cherenkov:A:B`; no table files |
| `--lightpulse S` | both | `Simple:0` | `Simple:W[,D]`, `Gauss:S`, `Gauss-rms:S`, `Gauss-fwhm:F`, `Exp:T`, `Complex:A,B,C[,D]` |
| `--angular-distribution S` | both | `isotropic` | `isotropic`, `Parallel`, `Cone:DEG`, `FilledCone:DEG`, `Gauss:DEG`, `Gauss-fwhm:DEG`, `Rayleigh` |
| `--run N` | both | `1` | run number |
| `--altitude M` (`-a`) | both | `2150` dish, `2100` secondary | observation level [m] |
| `--atmosphere N` (`-A`) | both | `1` | atmospheric profile number (file names not supported) |
| `--seed N` | both | `1` | Philox seed |
| `-o FILE` | both | stdout | output file |
| `--events N[,...]` | dish | `1` | events per photon level; a short list repeats its last value |
| `--frame sky\|camera` | dish | `sky` | `camera`: photons on the camera plane z = 0, see below |
| `--xy X,Y` (`--flasher-xy X,Y`) | dish | `0,0` | flasher offset from the axis [cm] |
| `--distance D` (`-z D`) | dish | `1600` | flasher distance above the telescope centre [cm] |
| `--camera-radius R` (`-r R`) | dish | `120` | fiducial camera radius [cm] |
| `--events N` | secondary | `1` | events (one value) |
| `--flasher-xy X` (`--xy X`) | secondary | `14.62` | flasher offset on the diagonal [cm] |
| `--flasher-depth D` (`-z D`, `--z D`) | secondary | `2.46` | flasher depth behind the camera front [cm] |
| `--flasher-inclination A` (`--inclination A`) | secondary | `0` | outward inclination [deg] |
| `--mirror-camera-distance D` (`--mirror-distance D`) | secondary | `51.06949` | focal front to secondary [cm] |
| `--delays T1,..,T4[,T5]` | secondary | `0,0,0,0,0` | flasher delays and the common-flash delay [ns] |
| `--fire all\|alternate\|demo\|1..4` | secondary | `all` | which flashers fire |

## Geometry

- dish: the emitter sits at `(x, y, distance)` above the telescope centre, which is the
  origin at observation level. Emitter pointing `(0, 0, -1)`. Photon count
  `photons / bunchsize`. For isotropic emission beyond `1.02 R` the source uses a
  `SphereTarget` of radius R (cone with `sin(theta) = R / d`). Fixed-angle profiles keep all
  photons. Other profiles use a `DiscTarget` on the sphere's plane (weights).
- secondary: origin at the secondary centre, z along the optical axis. Flashers at
  `(+-xy, +-xy, mirror distance + depth)` (flashers 1 to 4), pointing
  `DirectionVect(-90 deg + inclination, atan2(y, x))`. Count `photons / bunchsize` per
  flasher. The flashers lie inside the fiducial sphere (120 cm), so no target: 4 pi emission.
- secondary fire modes: `all` (four emissions), `alternate` (flasher `event mod 4`), `demo`
  (four individual flashes, then all with `delays[4]`), `n` (one flasher). `demo` with
  all-zero delays uses the demo values `0,10,20,30,50` ns.

`--frame camera` (dish): frame C. The photons run in a straight line from the flasher to the
plane z = 0 of the camera, without atmosphere. Time is emission time plus path / c. Count and
directions equal the sky-frame output for the same seed; positions lie in a disc of radius
`R / sqrt(1 - (R / d)^2)`. A photon that does not travel towards the plane is an error.
Pipeline: see [Using obdeect sources with sim_telarray](../sim_telarray.md).

## Times and frames

- Times are emission times (plus the delay for secondary). The atmosphere stage adds path
  delays.
- Positions are emission points, not impacts on a plane. No refraction (applied by the
  atmosphere stage). Output is a photon stream.
- `--delays` with fewer than 4 values is an error.
- `--xy` and `--flasher-xy` are both accepted for both types, with the meaning of the type
  (a pair for dish, one value for secondary).

## Checks

- `ctest`: `test_flasher` (dish: count, cone and direction moments, zem, non-default spectrum,
  pulse, bunch size, event repeats; secondary: flasher positions, direction moments, counts,
  delays, fire patterns, alternate events; type rules).
- dish reference run (default options, 20 events; `validation/sources/compare_flasher-dish.sh`): photons
  per event 3519.7 (reference) vs 3520.5 (obdeect sources); cz mean -0.998593 +- 3e-6 vs -0.998583 +- 3e-6
  (2.2 sigma); wavelength 400 nm both; reference time 53.45 ns is the flight time.
- secondary reference run (`--photons 1e5 --events 1`; `validation/sources/compare_flasher-secondary.sh`):
  399939 bunches (reference) vs 400000 (obdeect sources). The reference cz is always negative (sign
  lost), so the cz comparison is not valid for this isotropic emitter; cx, cy agree within
  1 sigma.
