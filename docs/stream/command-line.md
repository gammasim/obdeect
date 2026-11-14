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
- `obdeect stream --version` prints the package version and the photon-stream format version.

Every command accepts `--help`. For complete pipelines see [](../user-guide/pipelines.md).

## Camera map and telescope frames

The camera map lists the camera-frame telescopes, one per line:

```text
# telescope_id mirror_class(1 single | 2 dual) camera_rotation_deg
2 1 0.0
```

Telescopes not listed are ground-frame (sim_telarray `Bypass_Optics=0`). Listed telescopes with
camera-frame photons are converted to the bypass convention and need `Bypass_Optics=1` (class 1)
or `2` (class 2). `obdeect stream to-simtel` prints the value for each telescope.

## ATMPROF and synthesised headers

`--atmprof FILE` reads a CORSIKA/sim_telarray atmospheric table (altitude [km], density, thickness,
n-1 per line, `#` comments). It is written as ATMPROF (1216) before the first event, with the
profile number that EVTH names (99 for synthesised headers). A stream with its own ATMPROF block
rejects the option.

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

The bunch position is where the photon's line crosses z = 0 (path 0.625 m here, also for a photon
moving away from the plane). The time is `t` plus that path over the speed of light; `zem` is the
observation level plus `z`. Without a `cz` column the photon moves towards the plane. The camera
rotation is not used, and TELPOS z is set to 0 as for camera-frame telescopes.

```console
$ cat cam.txt
1 2 0.0
$ obdeect sources flasher --type secondary --photons 1e6 | obdeect stream to-simtel --camera-map cam.txt | sim_telarray -c <config> -C Bypass_Optics=1 -
obdeect stream to-simtel: telescope 1: telescope frame, mirror class 2, Bypass_Optics=1
```

Errors: T photons for a telescope without map entry or with `mirror_class` 1; a photon with
cz = 0; camera-frame and T photons for the same telescope in one run. The efficiency status of
such a telescope is checked as for camera frame.

## Efficiency checks

`obdeect stream to-simtel` refuses a stream whose [efficiency status](efficiency-status.md) says
that a factor is still to be applied when sim_telarray does not apply it in that frame. The
error names the factor:

```text
obdeect stream to-simtel: efficiency status of telescope 2 (camera frame): factor
'telescope_transmission_constant' is pending, but sim_telarray does not apply it in this frame;
apply it upstream (state applied); pass --allow-pending telescope_transmission_constant to accept it as it is
```

Use `--allow-pending FACTOR` (repeatable) only when the factor is applied elsewhere in the chain,
for example by the camera module. If the QE envelope was applied upstream (state `envelope`) and
the sim_telarray QE table is divided by that envelope, declare it with `--scaled-factor`:

```console
$ obdeect stream to-simtel --camera-map cam.txt --allow-pending telescope_transmission_constant run.pst > sim.simtel
$ obdeect stream to-simtel --camera-map cam.txt --scaled-factor qe run.pst > sim.simtel
```

In Python: `allow_pending=[...]` and `scaled_factors=[...]` of `obdeect.stream.stream_to_simtel`.
The factors that sim_telarray applies in each frame are listed in [](efficiency-status.md).
