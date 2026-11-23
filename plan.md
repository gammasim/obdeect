# Plan: photons in all directions

Goal: every stage (stream, sources, propagation, ray tracing) accepts photons
moving in any direction: downward, upward and horizontal.

## Decisions

- D1: horizontal photons are supported. G gets a second representation,
  `sphere` (below); `plane` stays the default. The run header records which
  one a run uses, so columns never depend on photon content.
- D2: sources with `to_plane_z0 = true`: photons that never reach z = 0 are
  dropped and counted (`not_on_plane`, run summary; no error).
- D3: sim_telarray cannot take upward or horizontal photons. `to-simtel`
  rejects them with an explicit error; it never writes them as downward.
- D4: ray tracing phase 1, then phase 2 (non-sequential); both done.
- Status: phases 1-4 in PR 80, phases 5-7 in PR 81 (validation PRs 17, 18).

## Conventions

- Direction `(cx, cy[, cz])`; without `cz` a photon moves downward. Upward or
  horizontal photons need the `cz` column. Stages keep `cz` if the input had it.
- Reachability is geometric (does the bent path reach the telescope sphere),
  not the sign of `cz`. With refraction a near-horizontal ray can turn over.
- The ground is opaque: emission or a path below the observation level is lost.
- Loss counts in the run summaries (not the efficiency status, which holds
  factors only): drop reasons, `not_on_plane`, ray-tracing statuses.

## G frame, sphere representation

- Row = the point where the photon starts its path inside the telescope
  sphere, relative to the sphere centre, with its direction and the time there.
- Emitted outside the sphere: the entry point, |r| = R.
- Emitted inside the sphere: the emission point itself, |r| < R. A muon track
  crossing the sphere gives both kinds of rows in one event.
- Run key `ground_representation` = `plane` (default) or `sphere`; format 1.2
  adds column `path` (distance from emission to the row point, m).
- Bending: if the ray crosses the plane with |cz| >= 0.1, plane-mode bending
  and the integrator from the crossing to the sphere; else a 3D ray integrator
  (RK4 on d(n d)/ds = grad n, flat and spherical; turning points included).
- Extinction: plane formulas for |u| >= 0.1, else the chord integral of k.
- Plane mode unchanged (CORSIKA, sim_telarray); `to-simtel` converts sphere
  to plane for downward photons emitted outside the sphere only.

## Phases

1. Stream: helpers, format 1.1, explicit rejections in `to-simtel` (also fixes
   the lost `cz` sign for T-frame export).
2. Propagation, plane mode, upward photons:
   - `sky_to_ground`, `ground_to_sky`, `validate_ground_chunk`.
   - Extinction over [min, max] of emission and telescope altitude.
   - Refraction: upward rows (short paths below the plane) by Snell quadrature.
   - Scattering: upward primaries; scattering points below the telescope
     sending light upward; keep `cz`.
3. Sources: `to_plane_z0` (D2), `plane_wave` any non-horizontal direction,
   presets (laser, lidar, beam scan) any beam direction.
4. Ray tracing phase 1 (see below).
5. All directions, sphere mode:
   - 5a stream (format 1.2, run key, sphere rows, `to-simtel` conversion);
     5b propagation (sphere rows, integrator, extinction, scattering);
     5c sources (G sphere rows, horizontal plane waves and beams);
     5d ray tracer starts at the sphere row; the sphere must enclose all
     optics and structure (checked when the optical model is compiled).
6. CORSIKA 3D bunches (1217) imported as sphere runs; fixture written with
   the bernlohr library (provenance file).
7. Ray tracing phase 2 (required).

## Ray tracing phase 1

Downward photons keep the current code path. Photons from below start at an
input plane below the optics; first hit among all surfaces and obscurers:

- reflecting face of the primary: continue the sequential path;
- mirror back, obscurer, camera housing: `blocked`, surface recorded;
- nothing: `missed`;
- reflecting face out of sequence (secondary first, focal surface directly):
  counted loss with its own status (`non_sequential`, `direct_focal`).

Phase 1 counts give the light that phase 2 recovers.

## Ray tracing phase 2

- After every interaction search all surfaces again until detected, absorbed
  or escaped; fixed bounce limit, counted when reached.
- Random numbers keyed additionally on the bounce number.
- Acceleration: per-surface bounds and a facet grid in the compiled optical
  model; no allocation or virtual dispatch per photon.
- The sequential path is a performance hint only: test the expected next
  surface first, then the bounds of the others; results must be identical.

## Validation

In obdeect: unit tests per change (analytic first), full checks, FMA and
sanitizer builds, Towncrier fragments. In `../obdeect-validation` (cases,
adapters, suites, changelog), with each phase, against independent references:

- Regression: all existing suites unchanged (downward photons bit-identical
  where the code path is unchanged); performance suite shows no slowdown.
- `stream.direction-roundtrip`: upward and horizontal `cz` survive `.pst`
  write and read exactly; `to-simtel` fails on them with the expected message.
- `propagation.upward-geometry`: straight-line position and time from emission
  below the plane, exact to float32 rounding.
- `propagation.upward-refraction`: bending and delay against a scipy Snell
  integral (and reversed downward ray); 1e-6 m, 1e-6 ns.
- `propagation.upward-extinction`: plane and curved transmission against an
  independent table interpolation and slant integral (reciprocity).
- `propagation.direction-losses`: mixed population; exact counts per loss
  category; photons in = photons out + losses.
- `propagation.upward-lidar`: upward laser, single scattering into an
  off-axis telescope against the numerical lidar equation (total and per
  time bin, within Monte Carlo errors).
- `determinism.mixed-directions`: same output for chunk size, threads, order.
- `sources.isotropic-to-plane`: kept fraction 1/2 (binomial), counted rest,
  radial density H / (r^2 + H^2)^(3/2).
- `optics.upward-absorption`: uniform light from below on LST, MST, SST:
  `blocked` fraction against the projected mirror area (independent geometry
  from the optical model), status counts sum to the input, none detected.
- Phase 5: horizontal light (sphere entry points against analytic
  intersections; extinction L * k(h); turning points against scipy), muon
  crossing the sphere (inside and outside rows, total yield vs analytic).
- Phase 7: phase 1 and 2 agree for sequential photons; stray light analytic.
