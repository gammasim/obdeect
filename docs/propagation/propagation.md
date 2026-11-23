# Atmosphere stage

`obdeect propagation run` (command line) and `obdeect.propagation.propagate` (Python) apply the
atmosphere to a photon stream. They move photons from their emission point to the telescope,
then apply the atmospheric effects of each chunk in a fixed order. Each effect is recorded in the
output [efficiency status](../stream/efficiency-status.md).

Use them for photons in the sky frame: air showers read with `obdeect stream from-corsika` and
the sky-frame [light sources](../sources/index.md). The output is in the ground frame (G).

## Stages

Per telescope and event, without scattering:

| # | Step | Frame | Reference |
|---|---|---|---|
| 1 | Straight line from the emission point to the horizontal plane through the telescope sphere centre; bending (refraction) unless the run says it is applied; air delay | sky → G | [](refraction.md) |
| 2 | Wavelength from 1/lambda² for rows without one (drawn before step 1 with `--dispersion-reference-nm`, which needs the wavelength of each bent photon) | G | [](extinction.md) |
| 3 | Extinction p_atm (plane-parallel airmass by default, zem = emission altitude; curved slant depth with `--curvature`) | G | [](extinction.md), [](curvature.md) |
| 4 | Efficiency bound E = product of the `--envelope` tables; thinning or weighting | G | [](extinction.md) |

With `--scattering`, one more step runs per event, before the steps above:

| # | Step | Frame | Reference |
|---|---|---|---|
| 0 | Single scattering of the primary photons (downwards and upwards, not horizontal; all blocks, or `--scattering-primaries`): wavelengths first, then scattered photons towards every telescope sphere, in any direction (upwards from points below the telescope), in their telescope's block (kept telescopes only with `--telescopes`) | sky | [](scattering.md) |
| 1–4 | Primaries: as above (p = p_atm E) | sky → G | |
| 1s | Scattered photons: sky → G, same run geometry (refraction, air delay) | sky → G | |
| 4s | Scattered photons: p = E exp(-od_focus), thinning or weighting (own random stream, see [](../stream/random-numbers.md)) | G | |

Scattered photons skip steps 2 and 3. Their transmission from the scattering point to the
telescope sphere is already in w (`T(P → C)` of [](scattering.md)), so the atmosphere is not
applied twice. Step 4s adds only the in-telescope part exp(-od_focus(lambda)) of their telescope
(the same `od_focus` as step 4, so the path length applies to them too; with a path of 0 m the
factor is 1). The primaries are scattered before their own extinction and are otherwise unchanged.

An event is held until its `event_end` record, because the scattered photons of a telescope come
from the primaries of all telescopes of the event. Memory: one event of the input.

To add an effect, add a step in `Propagator::photons` (`cpp/include/obdeect/propagation/propagator.hpp`)
and an entry in the list of efficiency factors built there.

## Frames and positions

- Sky: x (north), y (west), z above the observation level; direction cosines (cx, cy, cz). Without
  a cz column a photon moves downwards (cz = -sqrt(1 - cx^2 - cy^2)). With a cz column, cz < 0 is
  downwards, cz > 0 upwards and cz = 0 horizontal; the column must match cx, cy (either sign).
- Telescope i (1-based `telescope_id`) is at the i-th entry of `telescope_x_m`, `telescope_y_m`,
  `telescope_z_m` (z above the observation level).
- G: x, y relative to the telescope centre; z is removed. A cz column of the input is kept (its sign
  follows the direction of the output photon); a stream without cz has none in G.
- Emission altitude a.s.l. = `observation_level_m` + z. Telescope altitude a.s.l. = the same plus
  the telescope z.

A downward photon reaches the plane if it starts at or above it. An upward photon reaches it if it
starts between the observation level and the plane (emitted below the plane, it travels up to it;
see [](refraction.md#upward-photons)). Horizontal photons (cz = 0) never reach the plane. The other
photons are dropped and counted by reason (see [Dropped photons](#dropped-photons)). Their rows stay
in the chunk with w = 0, so the random draws of the other photons do not change.

## Dropped photons

The stage counts the dropped photons by reason. Each dropped photon has exactly one reason:

| Reason | Photons |
| --- | --- |
| `moving_away` | downward photons starting below the plane; upward photons starting above it |
| `horizontal` | cz = 0 (cx^2 + cy^2 = 1) |
| `below_ground` | starting below the observation level (z < 0), either direction |
| `turns_over` | upward photons turned back by refraction (sin theta reaches 1) before the plane |
| `not_bendable` | bending by quadrature failed otherwise (not expected for valid input) |
| `misses_level` | `spherical` mode: the straight line does not reach the telescope level |
| `misses_sphere` | sphere representation: the path does not enter the telescope sphere (also: the integrator leaves the atmosphere or the step limit is reached) |

The summary of a run has `dropped` (the total) and `dropped_by_reason`, a dict with these seven keys.
The sum of the values equals `dropped`. The efficiency status does not record the dropped count; it
is in the summary only (and in the output of the executable as its `dropped` total).

Time: without a profile, t + s/c. With a profile and no refraction in the stream, the air delay
`straight_path_air_delay_ns` is added. With refraction (`refraction_applied` = 0 and a profile),
time and direction come from `bend` at the telescope plane ([](refraction.md)). In `spherical` mode
the bending and air delay are the spherical integrals of [](curvature.md) instead.

## Efficiency status

One `efficiency_status` record follows each `run` record (before the first record that is not an
efficiency-status record). Example:

```text
factors                 atmosphere,wavelength,scattering,qe
atmosphere.state        applied
atmosphere.source       atm_trans_2000_1_10_0_0_2000.dat
atmosphere.pathlen_m    [28, 0]
atmosphere.wavelength_interpolation linear
atmosphere.dispersion_reference_nm 500
wavelength.state        applied
wavelength.range_nm     [250, 700]
scattering.state        applied
scattering.source       scattering_exponential.dat
scattering.points_per_photon 8
qe.state                envelope
qe.source               qe.dat
qe.wavelength_nm        [200, 800]
qe.value                [0.9, 0.5]
```

- `atmosphere`: applied; source = transmission file; `pathlen_m` = path inside each telescope [m],
  in telescope_id order (entry 0 is telescope 1).
- `atmosphere.wavelength_interpolation`: `linear` (default) or `floor`, the mode of the table
  ([](extinction.md)). Not written for geometry-only runs.
- `atmosphere.dispersion_reference_nm`: only with `--dispersion-reference-nm`; the wavelength [nm]
  of the profile's n-1 ([](refraction.md#dispersion)).
- `atmosphere.geometry`, `atmosphere.earth_radius_m`: only with `--curvature corsika` or
  `spherical` (`corsika-curved` or `spherical`, and the radius [m]); absent in mode `none`
  ([](curvature.md)).
- `wavelength`: applied; rows without a wavelength get one from the range.
- `scattering`: only with scattering; applied; source = scattering table; points per photon.
- `scattering.primaries`: only with `--scattering-primaries`; the telescope ids (1-based) whose
  photons are scattered.
- Each efficiency bound: state `envelope` (its value is part of the thinning, not applied in full),
  with its table.

The stage refuses an input stream whose efficiency status shows `atmosphere`, `scattering` or an
efficiency bound as applied (`pending` or absent is fine). Input records are passed through before
the stage's own record. Readers merge the efficiency-status records of a stream in order. The
record format is described in [](../stream/efficiency-status.md).

## Output

- Thinned (default): one row per surviving photon, w = 1 ([](extinction.md)). `first_index` counts
  output rows per telescope and event. Input weights must be at most 2048 photons per row.
- Weighted: all rows with w multiplied by p; rows with w = 0 removed. Any finite weight >= 0 is
  accepted, e.g. a laser pulse of 1e17 photons in 1e4 rows.
- Column p is written only with `--store-p` (`store_p`).
- `run`: `refraction_applied` = 1 when a profile is given.
- `telescope`: frame set to G.
- Scattered photons of telescope t are in its block, after the primaries. A telescope with
  scattered photons but no input block gets its own `telescope` and `telescope_end` records
  before `event_end`.
- `telescope_end`: `photons` is the sum of w of the output.
- Other records (event, blob, end records) are passed through. The output stream record comes
  from the writer; the input one is not copied.
- With a telescope selection (`--telescopes`), the blocks of other telescopes are removed: their
  `telescope` and `telescope_end` records, photons and telescope status. The format allows any
  subset of telescope blocks in an event (`telescope photons* telescope_end`, repeated). Their
  photons are not decoded. Run and event records are unchanged.
- With scattering, the blocks of other telescopes are still read for primaries (when they are
  primaries) but not written, and they are not scattering targets. The summary counts only the
  written photons.

## Command line

```console
$ obdeect propagation run --transmission atm_trans_2000_1_10_0_0_2000.dat --site-altitude 2000 \
    --profile atmprof.dat --envelope qe=qe.dat --envelope mirror=mirror.dat \
    --wavelength-range 250 700 --seed 1 -i in.pst -o out.pst
obdeect propagation run: photons in <n>, out <m> (weight <w>), dropped <d>
```

| Option | Meaning |
|---|---|
| `--transmission FILE` | transmission table, extinction-table format (required, except with `--no-extinction`) |
| `--no-extinction` | geometry only, see below |
| `--site-altitude M` | observation level [m a.s.l.] (required) |
| `--profile FILE` | atmospheric profile (atmprof format): air delay and refraction |
| `--envelope NAME=FILE` | efficiency bound of a later factor, two-column text (nm, value); repeatable |
| `--wavelength-range LO HI` | wavelengths from 1/lambda² [nm] for rows without one; LO = HI is a fixed wavelength; default: run keys |
| `--wavelength-interpolation linear\|floor` | how the transmission table is read at a wavelength (default linear, continuous in lambda); `floor` reproduces sim_telarray; needs `--transmission` ([](extinction.md)) |
| `--dispersion-reference-nm NM` | wavelength [nm], 200 to 2000, of the profile's n-1: bent photons are scaled to their wavelength; needs `--profile`; rejected for a run with `refraction_applied = 1` ([](refraction.md#dispersion)) |
| `--seed N` | Philox seed (required) |
| `--mode thinned\|weighted` | output mode (default thinned) |
| `--store-p` | keep column p |
| `--scattering FILE` | single scattering with this scattering table (off by default) |
| `--scattering-points N` | scattering points per photon, 1 to 1024 (default 8; needs `--scattering`) |
| `--curvature none\|corsika\|spherical` | Earth curvature of the atmosphere (default none), see [](curvature.md) |
| `--earth-radius-m M` | Earth radius for the curved modes [m] (default 6371315, the CORSIKA CURVED radius) |
| `--ground-representation plane\|sphere` | output representation of frame G (default plane); `sphere` rows are points on the telescope sphere ([](#sphere-representation)); needs `telescope_r_m` and `--no-extinction` for now |
| `--pathlen M` | path inside each telescope [m], finite and ≥ 0 (default 0) |
| `--telescope-pathlen ID=M` | path inside telescope ID (1-based) [m]; overrides `--pathlen`; repeatable |
| `--threads N` | threads for the rows of each chunk, 1 to 256 (default 1); the output does not depend on N |
| `--scattering-primaries ID[,ID...]` | only the photons of these telescopes (1-based ids) are scattered; default all blocks; needs `--scattering`; ids must exist in the run |
| `--telescopes ID[,ID...]` | keep only these telescopes (1-based ids); other telescope blocks are dropped and not decoded (with `--scattering`: the kept telescopes are the scattering targets) (default all) |
| `-i FILE\|-`, `-o FILE\|-` | input and output stream (default stdin, stdout) |

Path 28 m for all telescopes, 0 m for telescope 2:

```console
$ obdeect propagation run ... --pathlen 28 --telescope-pathlen 2=0
```

Keep telescopes 1 and 3 only. The summary counts their photons:

```console
$ obdeect propagation run ... --telescopes 1,3 -i in.pst -o out.pst
```

The output for the kept telescopes is the same as the full run with the other blocks removed.
The summary goes to stderr.

Scatter only the photons of telescope 1 (a collector sphere, for example); telescopes 2 and 3 are
targets:

```console
$ obdeect propagation run ... --scattering scattering_exponential.dat --scattering-primaries 1 \
    -i in.pst -o out.pst
```

The scattered photons of telescopes 2 and 3 are the same as in a run whose input holds only the
block of telescope 1.

## Python

```python
import obdeect.propagation

summary = obdeect.propagation.propagate(
    "in.pst",
    "out.pst",
    transmission="atm_trans_2000_1_10_0_0_2000.dat",
    site_altitude_m=2000.0,
    seed=1,
    profile="atmprof.dat",
    envelope={"qe": "qe.dat", "mirror": "mirror.dat"},
    wavelength_range=(250.0, 700.0),
    mode="thinned",
    scattering="scattering_exponential.dat",
    scattering_points=8,  # optional
    scattering_primaries=[1],  # optional, needs scattering
    wavelength_interpolation="linear",  # default; "floor" for sim_telarray
    dispersion_reference_nm=None,  # optional, needs profile; e.g. 500.0
    pathlen_m=28.0,
    telescope_pathlen_m={2: 0.0},  # optional
    telescopes=[1, 3],  # optional
    curvature="spherical",
    earth_radius_m=6371315.0,  # optional, see curvature.md
)
```

`summary` is a dict with `photons_in`, `photons_out`, `dropped`, `dropped_by_reason` (a dict of the
reasons of [Dropped photons](#dropped-photons)) and `weight_out`. Files are paths only (no `-`).

## Geometry only

`--no-extinction` (Python: `extinction=False`) moves photons from the sky frame to G and assigns wavelengths, and
nothing else: no extinction, envelope or thinning. Use it when a later program applies the
transmission itself, for example a light source feeding the camera simulation directly.

```console
$ obdeect propagation run --no-extinction --site-altitude 2000 --profile atmprof.dat --seed 1 \
    -i light.pst -o ground.pst
```

- Positions, directions and times are those of the normal run (same refraction and air delay).
- Photons that reach the telescope plane are kept with their weight w; photons that do not (w = 0)
  are removed, as in weighted mode. There is no column p.
- The efficiency status records `atmosphere.state` as `pending` (no source, no `pathlen_m`), so a later stage
  applies it. `wavelength` is `applied`.
- No transmission table and no `zem` are needed.
- Rejected: `--transmission`, `--envelope`, `--mode thinned`, `--scattering`, `--scattering-points`,
  `--store-p`, `--pathlen`, `--telescope-pathlen`. `--mode weighted` is accepted.

## In-telescope path

The extinction inside the telescope covers the path from the telescope sphere plane to the focus.
With path length l (`pathlen_m`), the optical depth is tau(h_T + l) − tau(h_T), with h_T the
telescope altitude. `--pathlen` sets l for all telescopes; `--telescope-pathlen`
(`telescope_pathlen_m`, 1-based ids) overrides it for single telescopes. Default 0, so existing runs
keep their result. An id that the run does not have is an error.

The path follows from the optical model (`obdeect.propagation.in_telescope_path_m`, from
`obdeect.ray_tracing`): the focal length for a single mirror, `|p0-s0|+|s0-f0|` for a dual mirror
(primary vertex p0, secondary vertex s0, focus f0). The primary vertex is the dish vertex at
`-mirror_offset`: an asphere stores it, and a segmented primary stores it as `vertex_z_m`. Models
compiled before that key was written are rejected; compile them again.

```console
$ obdeect propagation pathlen telescope_model.json
8.000000
```

The value is in metres. Use it as `--pathlen` (or `pathlen_m`) for the telescope:

```console
$ obdeect propagation run ... --pathlen "$(obdeect propagation pathlen telescope_model.json)"
```

## Sphere representation

The default frame G representation is the plane: a row is the photon's position on the telescope
plane (z = 0 relative to the telescope), and the telescope sphere is not used. The representation
`sphere` (`--ground-representation sphere`, Python `ground_representation="sphere"`, run key
`ground_representation = "sphere"`, format 1.2) keeps the three-dimensional position where the
photon's path inside the telescope sphere begins. The sphere has the centre C = (telescope x, y, z)
of the sky frame and the radius R = `telescope_r_m` of the run:

- emitted outside the sphere: the entry point, |r| = R (r is relative to C);
- emitted inside the sphere: the emission point, |r| < R, with path 0.

Columns: `x y z cx cy cz t zem w path`, plus `wavelength` and `p` as before. `path` is the distance
from the emission to the row point [m] (written for every row in sphere mode, 0 for dropped rows);
`cz` is always written. A muon track crossing the sphere gives rows of both kinds in one event.

Sky-frame photons are converted with these rules, in this order (`sky_to_sphere` in `obdeect/propagation/sphere.hpp`):

1. Emitted below the observation level (z < 0): `below_ground`. The ground is opaque, also for
   an emission inside the sphere.
2. Emitted inside the sphere: the emission point, path 0, time and direction unchanged.
3. No refraction (no refraction tables and no spherical bending): the straight line, the smallest
   s > 0 with |P + s u - C| = R. Lines that do not meet the sphere (or meet it only behind the
   emission) are `misses_sphere`. A line that reaches z < 0 before the sphere is `below_ground`.
   The time is t + s / c plus the straight air delay of the profile, if any.
4. Refraction: the hybrid rule ([](refraction.md#hybrid-rule)) when the plane-mode computation
   reaches the telescope plane with |cz| >= 0.1 there and the emission is farther than R from the
   plane; otherwise the 3D ray integrator ([](refraction.md#integrator)).

Rows that do not reach the sphere are dropped and counted in `dropped_by_reason` as above, with
w = 0 and zeros; the row index is kept, as in plane mode.

Limits of this stage (phase 5b-1): extinction and thinning of sphere rows are not implemented yet
(phase 5b-2, docs/propagation/extinction.md). A run in sphere mode with extinction (no
`--no-extinction`) is rejected with an explicit error. Scattering in sphere mode goes through the
same conversion; its photons reach the same rejection. Input runs in the sphere representation
(from sources or CORSIKA) are kept as sphere rows; sphere rows need the `z` and `cz` columns and a
run key `telescope_r_m`. Frame G photons of a plane run cannot be converted to sphere rows (sky-frame
input is needed).

## Limits

- Scattering needs one `array_index` per event (CORSIKA without array reuse, `CSCAT 1`); an event
  with several is rejected.
- Upward photons are bent by quadrature over their path below the plane (refraction, [](refraction.md#upward-photons)) and their extinction is the same integral as for downward photons, over the heights between the emission and the telescope. In `spherical` mode the path below the plane is treated as straight (the G plane is Cartesian); the bending uses the profile and the flat layers. Upward primaries are scattered up to the top of the scattering table, and scattered photons go up to a telescope centre above them ([](scattering.md)).
- Plane-parallel airmass, straight-line geometry and refraction by default (mode `none`). With
  `--curvature spherical` the airmass, bending and air delay are spherical ([](curvature.md)); with
  `corsika` CORSIKA has applied them.
- Scattering is plane-parallel and needs sky-frame input and `telescope_r_m` in the run keys; not
  available for photons already in frame G, and rejected with a curved mode.
- The in-telescope path defaults to 0 m; the optical model gives it ([](propagation.md#in-telescope-path)).
- Transmission tables hold one row per integer nm. The default linear mode interpolates between
  the rows; `floor` uses the row at the integer part of the wavelength.
- Efficiency-bound tables must cover the wavelength range (rounded outwards).
- In mode `none`, photons near the horizon (zenith beyond about 87 degrees, where the series of
  [](refraction.md) is not valid) are bent by quadrature and not dropped ([](refraction.md#near-the-horizon)).
- Dispersion needs photons that are not refracted yet and a wavelength range within 200 to 2000 nm.
  The travel time is scaled by the group factor, so the time at the reference wavelength is not
  unchanged by the option (the group index of the reference is taken from the Ciddor formula).
