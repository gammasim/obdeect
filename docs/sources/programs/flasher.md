# flasher

Flasher photon stream: a light flasher in front of a telescope.

```console
$ obdeect sources flasher --type dish > flasher.pst
$ obdeect sources flasher --type secondary --fire all --delays 0,0,0,0 > secondary.pst
```

- `--type dish` (default): flasher on the axis of a single-reflector telescope, in front of
  the camera. Frame `sky`.
- `--type secondary`: four flashers at the camera corners of a dual-mirror telescope, firing
  onto the secondary mirror. Frame `telescope` (T).

An option of the other type is rejected, e.g. `--delays applies only to --type secondary`.

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
| `--frame sky\|camera` | dish | `sky` | `camera`: photons on the camera plane z = 0 (see below) |
| `--xy X,Y` (`--flasher-xy X,Y`) | dish | `0,0` | flasher offset from the axis [cm] |
| `--distance D` (`-z D`) | dish | `1600` | flasher distance above the telescope centre [cm] |
| `--camera-radius R` (`-r R`) | dish | `120` | fiducial camera radius [cm] |
| `--events N` | secondary | `1` | events (one value) |
| `--flasher-xy X` (`--xy X`) | secondary | `14.62` | flasher offset on the diagonal [cm] |
| `--flasher-depth D` (`-z D`, `--z D`) | secondary | `2.46` | flasher depth behind the camera front [cm] |
| `--flasher-inclination A` (`--inclination A`) | secondary | `0` | outward inclination [deg] |
| `--mirror-camera-distance D` (`--mirror-distance D`) | secondary | `51.06949` | focal front to secondary [cm] |
| `--delays T1,..,T4[,T5]` | secondary | `0,0,0,0,0` | flasher delays and common-flash delay [ns] |
| `--fire all\|alternate\|demo\|1..4` | secondary | `all` | which flashers fire |

## Dish: notes

- The emitter sits at `(x, y, distance)` above the telescope centre. Pointing `(0, 0, -1)`.
- Photon count `photons / bunchsize`.
- Isotropic emission beyond `1.02 R` uses a sphere target of radius R (cone with
  `sin(theta) = R / d`). Fixed-angle profiles keep all photons. Other profiles use a disc
  target on the sphere's plane (weights).
- `--frame camera`: the photons run in a straight line from the flasher to the plane z = 0
  of the camera, without atmosphere. Time is emission time plus path / c. Count and directions
  equal the sky-frame output for the same seed; positions lie in a disc of radius
  `R / sqrt(1 - (R / d)^2)`. A photon that does not travel towards the plane is an error.
  Pipeline: see [sim_telarray](../../user-guide/pipelines.md).

## Secondary: notes

- Flashers at `(+-xy, +-xy, mirror distance + depth)` (flashers 1 to 4), pointing towards the
  secondary. The flashers lie inside the fiducial sphere (120 cm), so no target: 4 pi emission.
- `--fire`: `all` fires four emissions, `alternate` fires flasher `event mod 4`, `1`..`4` fires
  one flasher, `demo` fires the four flashers singly and then all with `delays[4]`.
  `demo` with all-zero delays uses `0,10,20,30,50` ns.
- Times are emission times plus the delay.

## Times and frames

- Times are emission times (plus the delay for secondary). The atmosphere stage adds path
  delays.
- Positions are emission points, not impacts on a plane. No refraction (applied by the
  atmosphere stage). Output is a photon stream.
