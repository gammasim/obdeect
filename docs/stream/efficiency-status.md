# Efficiency status

Many effects reduce the number of photons that reach the camera: the atmosphere, mirror reflectivity,
quantum efficiency and others. Each module records which of these factors it has already applied,
so that no factor is applied twice or forgotten. The record is `efficiency_status`; its metadata are:

- `factors` (string): comma-separated factor names, e.g.
  `atmosphere,qe,mirror_reflectivity,telescope_transmission`.
- per factor `<name>.state` (string): `applied`, `envelope` (only the envelope f^env was
  applied), or `pending`.
- per factor `<name>.source` (string): provenance (file name, model version).
- per factor with a λ table: `<name>.wavelength_nm` and `<name>.value` (float arrays) of the
  applied (envelope) values.

An *efficiency bound* (named `envelope` in options and states) is an upper bound f_k^env of a factor that a module has used early, to thin photons
before the exact factor is known. A module applying factor k with an envelope in the efficiency
status applies f_k / f_k^env and sets `<name>.state` to `applied`. A module never applies a factor whose state is `applied`.

Run status records (before the first event) apply to all telescopes. A telescope status (after
its `telescope` record, before its photons) applies to that telescope and overrides the run
status per factor. Within either scope, a factor may only advance `pending` → `envelope` →
`applied`; a repeated or backward state is invalid. Telescope updates are checked separately from
run updates before the telescope status overrides the run status.

`obdeect stream to-simtel` checks the efficiency status of each telescope against the frame of
its photons. It accepts a factor that is `applied` upstream unless sim_telarray applies it again
in that frame. It accepts `pending` only for factors that sim_telarray applies in that frame.
`envelope` is always an error for those factors: sim_telarray applies the full factor, not f/f_env.

| Factor | Camera frame (C) | Ground frame (G) |
|---|---|---|
| `qe`, `camera_transmission`, `camera_filter` | sim_telarray applies: pending ok | sim_telarray applies: pending ok |
| `mirror_reflectivity`, `telescope_transmission` | not applied by sim_telarray: error | sim_telarray applies: pending ok |
| `atmosphere` | must be applied upstream | any state; with `applied` sim_telarray needs an all-zero optical-depth table |
| `wavelength`, `scattering` | must be applied upstream | must be applied upstream |
| `telescope_transmission_constant` | not applied by sim_telarray: error | not applied by sim_telarray: error |
| any other name | error | error |

Any other factor (or unknown name) raises an error that names it, unless it is given to
`obdeect stream to-simtel --allow-pending FACTOR` or `allow_pending=` in Python.

`envelope` for a factor that sim_telarray applies is valid only with
`obdeect stream to-simtel --scaled-factor FACTOR` (`scaled_factors=` in Python): the sim_telarray
table of the factor then holds f/f_env, e.g. the QE table divided by the QE envelope applied
upstream. Every `--scaled-factor` needs a matching efficiency status in state `envelope`; a
missing status or any other state is an error.

The commands that read and write this record: [](command-line.md#efficiency-checks), [](../propagation/propagation.md#efficiency-status), [](../ray-tracing/photon-stream.md#efficiency-status).
