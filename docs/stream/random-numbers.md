# Random numbers

Every obdeect module draws random numbers from the same generator, with the same rules. This
page is the one place where they are described.

## What you can rely on

- **Same seed, same result.** A run with the same input, options and `--seed` gives the same
  output, bit for bit on one platform.
- **No dependence on how the work is split.** The result does not depend on chunk size, photon
  order or the number of threads.
- **Independent modules.** Each module uses its own stage id, so the atmosphere and the optics
  never draw the same numbers, even with the same seed.
- **Independent emissions.** Each emission of a run has its own seed, `source_seed(run seed, k)`
  for emission `k` in event order, so emissions never share draws. Presets draw their random
  choices from `stage_programs` (see below).

## Generator

Philox-4x32-10 (`obdeect/stream/random.hpp`) is a counter-based generator: every draw is a pure
function of a key, with no state to share between photons. The key combines:

| Field | Bits | Meaning |
|---|---|---|
| seed | 64 | the `--seed` of the command; for sources, the emission seed |
| stage | 32 | the module (table below) |
| photon index | 40 | index of the photon within its telescope and event |
| block | 24 | block of draws for one photon |
| telescope | 16 | telescope id |
| stream | 16 | purpose of the draw within the module |
| event | 32 | event number |

Out-of-range fields are rejected. Other generators agree in the probabilities, not in the
individual draws.

## Stage ids

Stage ids are defined once in `obdeect/stream/stages.hpp`; a `static_assert` checks that they are
distinct. The values are fixed: changing one changes the random numbers.

| Constant | Value | Module |
|---|---|---|
| `stage_sources` | `0x444C0001` | obdeect sources: photon and count draws, keyed on the emission seed |
| `stage_programs` | `0x444C0002` | obdeect sources presets: random choices (pixels, core positions) |
| `stage_ray_tracing` | `0x444C0003` | obdeect ray tracing |
| `stage_atmosphere` | `0x444C5546` | obdeect propagation: atmosphere |

Presets draw their random choices with numpy's Philox generator, keyed on (run seed,
`0x444C0002`), not with the generator above. Photon draws never use this stage.

A new module takes the next free value and adds it to `stage_ids`.

## Streams of each module

| Module | Streams | Notes |
|---|---|---|
| propagation | 1 wavelength, 2 thinning, 3 scattering | Thinning draws one Bernoulli trial per photon, so bunches are limited to 2048 photons; split larger bunches upstream. |
| sources | 1 direction, 2 target, 3 wavelength, 4 pulse, 5 count | The assignment is part of the [specification](../sources/spec/conventions.md). |
| ray tracing | one seed per (event, telescope) block, derived from `--seed` | The photon id is `(event_number << 32) \| photon_index`. |
