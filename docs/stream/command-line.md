# Command line

All commands are subcommands of `obdeect stream`. They read photon streams from files or
standard input (`-`) and write to standard output, so they combine with pipes.

| Command | Purpose |
|---|---|
| `cat` | concatenate streams into one, or print a summary |
| `tee` | copy a stream to standard output and to files |
| `from-corsika` | convert CORSIKA IACT EventIO output to a stream |
| `to-simtel` | convert a stream to sim_telarray input |
| `to-grisu` | convert a stream to GrOptics photon text ([](grisu.md)) |

```console
$ obdeect stream cat run1.pst run2.pst > all.pst
$ obdeect stream cat --summary all.pst
$ obdeect stream cat run.pst | obdeect stream tee copy.pst > out.pst
$ obdeect stream from-corsika run.iact > run.pst
$ obdeect stream to-simtel --camera-map cam.txt run.pst > sim.simtel
$ obdeect stream from-corsika < run.iact | obdeect stream to-simtel > sim.simtel
```

- `cat` writes one stream with a single `stream_end` record. `--summary` prints the records per
  type and the photons and weight per telescope instead.
- `tee` writes the same records to standard output and to each file.
- Errors of the output (a full disk, a closed pipe) end the command with exit status 1, also
  when they come at the final flush or the close of a file.
- `obdeect --version` (the top-level command) prints the installed component versions and the
  photon-stream format version.

Every command accepts `--help`. For complete pipelines see [](../user-guide/pipelines.md).

## Camera map and telescope frames

The camera map lists the camera-frame telescopes, one per line:

```text
# telescope_id mirror_class(1 single | 2 dual) camera_rotation_deg
2 1 0.0
```

Telescopes not listed are ground-frame (sim_telarray `Bypass_Optics=0`). Listed telescopes with
camera-frame photons are converted to the bypass convention and need `Bypass_Optics=1` (class 1)
or `2` (class 2). `obdeect stream to-simtel` prints the value for each telescope. Camera-frame
photons without a `zem` column (e.g. from `obdeect sources run led`) get `zem` = observation level
plus `z`.

TELPOS z of the listed telescopes is set to 0, because their photons are already at the camera.
For array runs, sim_telarray needs the real heights for the plane-wavefront delay between
telescopes; `--keep-telpos-height` writes the TELPOS block unchanged. Default: off.

```console
$ obdeect stream to-simtel --camera-map cam.txt --keep-telpos-height run.pst > sim.simtel
```

## ATMPROF and synthesised headers

`--atmprof FILE` reads a CORSIKA/sim_telarray atmospheric table (altitude [km], density, thickness,
n-1 per line, `#` comments). It is written as ATMPROF (1216) before the first event, with the
profile number that EVTH names (99 for synthesised headers). A stream with its own ATMPROF block
rejects the option.

Photons for sim_telarray need the wavelength column: sim_telarray samples a Cherenkov spectrum for
wavelength 0, so a stream without that column is rejected (a column of zeros is written as it is).
A CORSIKA run must end with its RUNE block before the next RUNH; a file that ends inside a run, or
is empty, is rejected.

A stream without CORSIKA blocks gets synthesised headers. The run end (RUNE) carries the number of
events written, so sim_telarray reports the run as complete. The input block (INPUTCFG, one line
`NSHOW 0`) is written at the start of the run, when the event count is not yet known:

```text
Run(s) completed as expected after 3 events (3 showers).
```

## Telescope-frame photons (dual-mirror telescopes)

A telescope listed with `mirror_class` 2 may carry telescope-frame (T) photons instead of
camera-frame photons, for example four flashers at the camera firing onto the secondary mirror.
The T frame has its origin at the secondary centre and z along the optical axis towards the sky.
Each photon is an emission point (`x`, `y`, `z`) with direction (`cx`, `cy`, `cz`) and emission
time `t`. `obdeect stream to-simtel` writes one bunch per photon at the plane z = 0 of that frame,
and sim_telarray needs `Bypass_Optics=1`:

```text
photon:  x, y, z = 0.1, 0.2, 0.5 m   c = (0.6, 0, -0.8)   t = 3 ns   observation level 2200 m
bunch:   x, y = 47.5, 20.0 cm   cx, cy = 0.6, 0   time = 5.085 ns   zem = 220050 cm
```

The bunch position is where the photon's line crosses z = 0 (path 0.625 m here). The time is `t`
plus that path over the speed of light; `zem` is the observation level plus `z`. Without a `cz`
column the photon moves downwards. The camera rotation is not used, and TELPOS z is set to 0 as for
camera-frame telescopes (unless `--keep-telpos-height` is given).

Ground- and telescope-frame photons must move downwards (`cz` < 0, or no `cz` column):
sim_telarray does not store the sign of `cz`, and the bunch format has no sign either. An upward
or horizontal photon in these frames is an error, never written as a downward one. Camera-frame
photons are already at the camera entrance, where sim_telarray uses only `cx` and `cy`; the sign
of `cz` depends on the side the camera faces and is not checked.

```console
$ cat cam.txt
1 2 0.0
$ obdeect sources run flasher --type secondary --flasher-offset 0.5 --flasher-depth 0.3 \
    --flasher-inclination 2 --mirror-camera-distance 16 --fiducial-radius 12.5 \
    --observation-level 2150 --spectrum fixed:wavelength_nm=400 --photons 1e6 \
    | obdeect stream to-simtel --camera-map cam.txt | sim_telarray -c <config> -C Bypass_Optics=1 -
obdeect stream to-simtel: telescope 1: telescope frame, mirror class 2, Bypass_Optics=1
```

Errors: T photons for a telescope without map entry or with `mirror_class` 1; a photon with
`cz` >= 0 (upward or horizontal) or a non-finite `cz` in the ground or telescope frame; camera-frame
and T photons for the same telescope in one run. The efficiency status of such a telescope is
checked as for camera frame.

## CORSIKA 3D photon bunches (1217)

`obdeect stream from-corsika` converts EventIO 3D photon bunches (block 1217, written by the CORSIKA
IACT library, `write_tel_photons3d` in `io_simtel.c`) to a sphere run. A bunch gives the arrival
point A = (x, y, z) relative to the telescope sphere centre (cm), the direction d = (cx, cy, cz)
and the distance `dist` of the emission point from A. The emission point is E = A - dist d, and the
bunch becomes its G sphere row ([format](format.md), "Photon chunk"), with R the telescope radius from
TELPOS:

- E on or inside the sphere (|E| <= R): the row is E, path 0.
- E outside: the row is the entry point E + s d, where s >= 0 is the smallest root of
  |E + s d| = R, and path = s. The line may enter the sphere after the arrival point A (s > dist).
- No root with s >= 0 (the line misses the sphere, or points away from it): the bunch is counted in
  `misses_sphere` in the run end and not written.

Time at the row: `ctime + (s - dist) / c` (vacuum c; `ctime` is the arrival time at A). Emission
altitude a.s.l.: `zem` = observation level + TELPOS z + E_z (flat geometry). The columns are
`x y z cx cy cz t zem w path wavelength`; `wavelength` is the bunch's wavelength (0 = not set).

Conventions and errors:

- The run key `ground_representation` is `sphere`. A run may not mix 1205 and 1217 blocks; such a
  run is an error. A run with 1205 blocks is a plane run and gets no `ground_representation` key (absent means
  `plane`), so its output is unchanged.
- 1217 needs RUNH (observation level) and a TELPOS entry for each telescope it names.
- `cz` must be consistent with `cx` and `cy` (a unit vector within the float32 tolerance).
  Negative bunch sizes, negative distances and non-finite values are errors, not dropped.
- Heights are not converted for CURVED runs, as for 1205 blocks; the run key `curvature_applied`
  tells whether CORSIKA applied the curvature.
- A top-level 1217 block with indices (999, 999) is kept as a blob, as for 1205.

The bunches of the test fixture `cpp/tests/data/photons3d.iact` (provenance in
`photons3d_provenance.txt` next to it) give these conversions by hand in `test_corsika.cpp`.

## Sphere runs (ground frame)

A run with `ground_representation` `sphere` (format 1.2, [](format.md)) has ground-frame rows on
or inside the telescope sphere instead of at the plane z = 0. `obdeect stream to-simtel` converts
them to sim_telarray's plane rows, which need a downward `cz`:

- A downward photon that starts on the sphere (emitted outside it) is moved along its line to the
  plane z = 0: `x`, `y` at the crossing, `t` plus the path over the speed of light (the path can be
  negative when the entry point is below the plane), `zem` and `cz` unchanged. The `z` and `path`
  columns are dropped.
- A photon that starts inside the sphere (emitted inside it) cannot be represented and is rejected:
  `ground-frame photon emitted inside the telescope sphere: not representable for sim_telarray`.
- Upward and horizontal photons are rejected as in plane runs.

The sim_telarray output is that of a plane run with the converted rows; the run key is not passed
on. A sphere run with rows of the plane form (no `z` column) is rejected. `telescope_r_m` must list
the telescope's sphere radius.

`obdeect stream to-grisu` does not convert sphere runs and rejects them.

## Efficiency checks

`obdeect stream to-simtel` checks each telescope's [efficiency status](efficiency-status.md)
against the frame of its photons (the rules are there). A violation names the factor and exits
with status 1:

```text
obdeect stream to-simtel: efficiency status of telescope 2 (camera frame): factor
'telescope_transmission_constant' is pending, but sim_telarray does not apply it in this frame;
apply it upstream (state applied); pass --allow-pending telescope_transmission_constant to accept it as it is
```

Camera- and telescope-frame photons need the atmosphere applied upstream, for example by
`obdeect propagation run`. A stream whose status lacks `atmosphere` is refused:

```text
obdeect stream to-simtel: efficiency status of telescope 1 (camera frame): factor 'atmosphere' is
missing, but sim_telarray does not apply it in this frame; apply it upstream (state applied), or
pass --allow-pending-atmosphere if sim_telarray is meant to apply it
```

Use `--allow-pending-atmosphere` only when sim_telarray is meant to apply the atmosphere in that
frame; it accepts an atmosphere that is missing or pending. Use `--allow-pending FACTOR` (repeatable)
only when the factor is applied elsewhere in the chain,
for example by the camera module. If the QE envelope was applied upstream (state `envelope`) and
the sim_telarray QE table is divided by that envelope, declare it with `--scaled-factor`:

```console
$ obdeect stream to-simtel --camera-map cam.txt --allow-pending telescope_transmission_constant run.pst > sim.simtel
$ obdeect stream to-simtel --camera-map cam.txt --allow-pending-atmosphere run.pst > sim.simtel
$ obdeect stream to-simtel --camera-map cam.txt --scaled-factor qe run.pst > sim.simtel
```

In Python: `allow_pending=[...]`, `allow_pending_atmosphere=True` and `scaled_factors=[...]` of
`obdeect.stream.stream_to_simtel`.
The factors that sim_telarray applies in each frame are listed in [](efficiency-status.md).
