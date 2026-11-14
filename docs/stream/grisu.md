# GrOptics and GrISU

GrOptics traces photons through the telescope optics of
the GrISU simulation chain, and GrISUDet simulates the camera. `obdeect stream to-grisu` writes a
photon stream as the photon text that GrOptics reads, so the obdeect atmosphere and light
sources can feed that chain, in place of sim_telarray.

`obdeect stream to-grisu` writes ground-frame photons as GrISU photon text. The
atmosphere must be applied upstream and the wavelength assigned, for example by
[`obdeect propagation run`](../propagation/propagation.md):

```console
$ zstd -dc run.corsika.zst | obdeect stream from-corsika \
    | obdeect propagation run --transmission atm.dat --site-altitude 2156 --seed 1 \
        --envelope qe=qe_envelope.dat -i - -o - \
    | obdeect stream to-grisu --telescope-config TLLOC.dat > run.grisu
```

Output (one header block, then `S` and `P` lines):

```text
* HEADF  <-- Start of header flag
* Produced by obdeect stream to-grisu 0.1.0 from a photon stream written by obdeect stream from-corsika 0.1.0
* CORSIKA run number: 5
* Energy range [TeV]: 0.1 10
* Spectral slope: -2
* Primary id: 14
* Primary zenith [deg]: 20
* Primary azimuth [deg, CORSIKA]: 0
* Observation height [m]: 2200
* DATAF  <-- end of header flag
R 1
H 2200
S 1 -50 -100 0 -0.34202014 25000 -1 -1 -1
P -20.5 -11.5 -0.5 +0 +20000 +12.5 400 3 7
```

- `R` is the global QE factor (`--qe`, default 1). It is written, not applied.
- `S E_TeV x y dcos dsin h_first -1 -1 -1`: one line per shower position (per array reuse), with
  core position [m] and primary direction cosines in the GrISU frame. `h_first` is the event key
  `first_interaction_m` [m a.s.l.]. `--event-details` adds `C h_first X_first shower_id`;
  `X_first` is -1 (not carried by the stream) and `shower_id` is the event number.
- `P x y dcx dcy zem t lambda 3 tel`: one line per photon. `x`, `y` [m] are the telescope position
  plus the photon position; `zem` the emission height [m]; `t` the arrival time [ns]; `lambda` the
  wavelength [nm, integer]; `3` the emitter type (unknown); `tel` the GrISU telescope id (1-based).
- A bunch with weight above 1 is an error unless `--expand-bunches` writes `round(w)` copies. A
  weight below 1 is an error: thin the photons upstream.
- Numbers use 9 significant digits. Values below 1e-8 in magnitude are written as 0.

Frames: the CORSIKA ground frame (x north, y west) maps to GrISU by `x_g = -y_c`, `y_g = -x_c`.
Azimuths map by `az_g = (3 pi / 2 - az_c) mod 2 pi`, with `az_c = atan2(cy, cx)` for photons; the
written direction cosines are `sin(ze) cos(az_g)` and `sin(ze) sin(az_g)`, with
`sin(ze) = sqrt(cx^2 + cy^2)`. The primary uses the event zenith and azimuth the same way.

`--telescope-config FILE` reads `* TLLOC id x y z offset px py` lines (x east, y north [m]; other
lines are comments). A stream telescope gets the `id` of the closest line within 0.5 m of its
position (transformed as above). Other telescopes are skipped and reported on stderr. Without the
file, `tel` is the stream telescope id.

In Python: `obdeect.stream.stream_to_grisu(run.pst, run.grisu, telescope_config="TLLOC.dat")`.
