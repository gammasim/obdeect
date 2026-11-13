# obdeect project analysis

Analysis of obdeect v0.7.2 (ray tracing only), kept for reference; names are those of v0.7.2.

Analysis of obdeect v0.7.2 (`main` at `e4ade6a`, 2026-11-06), seen as a
gamma-ray astronomer, an optics expert, a high-performance-computing (HPC)
expert and a release manager. Goal: the most precise and fastest code that is
easy to maintain, with smooth CI and good documentation.

This is a review, not a change: no tracked file was modified (`README.md`
included). Statements are marked **[verified]** when I ran or read the code to
confirm them, and **[to check]** when they follow from reading and need a
measurement or a second opinion before anyone acts on them.

## 1. Summary

obdeect is in very good shape for a project that is seven weeks old (251
commits since 2026-09-18). The architecture follows its own rules: a
standard-library-only header core, immutable compiled optical models,
counter-based random numbers (Philox) that make results independent of thread
count, block size and photon order, and a thin nanobind boundary. The
validation against sim_telarray is thorough (D80 within 0.1-1.6 %, effective
area within 0.3 %, air-shower photoelectrons within 0.5 %), and the
documentation is unusually honest about limitations.

The weak points are not in the core numerics. They are:

1. **Validation is single-reference and not reproducible from the repository.**
   Everything is compared with sim_telarray; scripts, containers and result
   files live in an external workspace.
2. **One known, uncorrected geometric inconsistency in the mirror-scatter
   code** (non-orthogonal rotation axes on curved facets).
3. **CI covers the happy path only:** no ThreadSanitizer for the hand-written
   thread pool, no clang-tidy despite the guide, no Windows or Intel-Mac
   builds, no performance guard, no coverage, no type checking.
4. **Release hygiene:** four different version numbers in the tree, three
   releases in three days, broken README links, no dependency-update
   automation.
5. **Performance headroom of 10-30 %** is available cheaply (serial
   compaction, compiler flags, bunch-path duplicate work) before any SIMD work.

### Scorecard

| Viewpoint | Grade | One-line reason |
| --- | --- | --- |
| Gamma-ray astronomer | A- | Production-grade parity with sim_telarray; accuracy limited by shared assumptions and missing as-built effects |
| Optics expert | B+ | Sound geometry; scatter-axis flaw, root-finding corner cases, no refraction or polarisation |
| HPC expert | A- | Deterministic, allocation-free hot path, near-linear scaling; compaction and flags leave easy wins |
| Release manager | B | Good tooling (towncrier, trusted publishing, versioned docs); drift and gaps in the matrix |
| Maintainability / docs | B+ | Clear code and docs; a few oversized headers, 12-argument kernel call, dead links |

### Top 12 actions, in order

| # | Action | Viewpoint | Effort | Gain |
| --- | --- | --- | --- | --- |
| 1 | Add a ThreadSanitizer CI job (pool, `trace_bunches`) | HPC / release | S | Catches races in the only hand-written concurrency code |
| 2 | Make validation reproducible: move scripts, model fixtures and reference numbers into the repo or a DOI-versioned data release | Astronomer | M | Turns claims into regression tests |
| 3 | Remove mirror-realisation noise from the sim_telarray comparison by feeding identical per-mirror realisations to both codes | Astronomer / optics | M | D80 agreement from about 1 % to about 1e-4 |
| 4 | Orthogonalise the scatter axes against the hit normal (`trace.hpp:134`) | Optics | S | Removes a known 1.9e-2 rad geometric inconsistency |
| 5 | Parallelise or remove the serial `detected_only` compaction | HPC | S | Measured 30 % slower than writing all rows |
| 6 | Try `-fno-math-errno` (and `-fno-trapping-math`) on GCC/Linux wheels | HPC | S | Likely 5-15 %, results unchanged [to check] |
| 7 | One version source of truth; fix CMake 0.2.0, CITATION 0.1.0, fallbacks | Release | S | Correct citations and metadata |
| 8 | Fix README links to `USER_GUIDE.md` / `PHYSICS.md` (suggested, not applied) | Docs | S | Dead links on the front page |
| 9 | Replace `trace_photon`'s 12 positional parameters with an options struct | Maintainability | S | Removes a bool-ordering bug class |
| 10 | Add Dependabot or Renovate for Actions and pre-commit; pin Actions consistently | Release | S | Supply-chain hygiene |
| 11 | Add Windows and Intel-macOS to the test or wheel matrix, or state them unsupported | Release | M | Honest portability claim |
| 12 | Add a deterministic work-counter performance test (Newton steps, node visits) per real model | HPC | M | Catches algorithmic regressions on any runner |

Effort: S under a day, M a few days.

## 2. What I ran

| Check | Result |
| --- | --- |
| `cmake --preset release` + build (12 targets) | OK, 2.7 s wall on 12 cores |
| `ctest` (6 native test executables) | 6 of 6 pass, 4.2 s |
| `pytest python/tests` (scratch venv, fresh non-editable install) | 99 of 99 pass, 8.9 s |
| `ruff format --check python`, `ruff check python` | clean |
| `clang-format --dry-run --Werror` over `cpp/` | clean |
| `clang-tidy` | **not available here**, not run |
| Sphinx docs build, pre-commit, sanitizer and FMA presets | **not run** |
| Throughput probe (section 5.1) | run, numbers below |

**Environment finding [verified]:** the first `pytest` run failed at collection
with `CalledProcessError: cmake --build` and `No module named nanobind`. An
editable install made earlier in another virtual environment (with
`editable.rebuild = true`) tried to rebuild on `import obdeect` and could not
find nanobind, because pip's isolated build environment no longer exists.
`AGENTS.md` and `docs/installation.md` give `pip install -e '.[dev]'`; with
rebuild-on-import that needs `--no-build-isolation` and the build requirements
installed in the environment. Either document that, or make the failure
legible: the `execute_process(... -m nanobind --cmake_dir)` at
`CMakeLists.txt:47` only prints a CMake traceback.

## 3. Gamma-ray astronomer

### 3.1 What is strong

- Quantified agreement with sim_telarray over twelve telescope models, 0-4 deg
  field, plus measured-PSF checks (LSTN-01 D80 3.050 vs 3.052 cm).
- End-to-end air-shower check (1 TeV gammas, zenith 20/40/60 deg, five
  telescopes): all nine points within 0.5 % and 1.5 sigma, with timing (KS test)
  and image-shape agreement.
- `docs/limitations.md` names what is not modelled and sizes the effect (for
  example the camera-window focus shift, 1-3 mm for 3-10 mm PMMA).
- Every result carries the optical model's SHA-256; the model is rejected if
  its content does not match its hash (`python/obdeect/model.py:34`).

### 3.2 Findings

**A1. Validation is equivalence, not correctness.** All twelve comparisons use
sim_telarray as the reference. A convention shared by both (for example the
scatter definition, section 4) cannot show up. The comparison already found one
such convention (point-source distance measured from the dish vertex), which
shows the method works, but only for differences. Add reference cases whose
answer is known independently: a perfect parabola on axis (zero spot), the
Davies-Cotton coma and astigmatism formulas at small field angles, a
Schwarzschild-Couder design with a published spot diagram. Some of this exists
as "photons in analytic optics" in `cpp/tests/test_trace.cpp`; extend it to
quantitative aberration coefficients.

**A2. The quoted D80 spread is realisation noise, not code difference.** The
1.6 % D80 deviation for MST and 0.6 % for LST comes from each backend drawing
its own random mirror errors (`python/obdeect/compiler.py:364`, Python's
`random.Random` seeded per facet). The documentation says so, but it also means
the comparison cannot detect a 0.5 % error. If sim_telarray's realised mirror
tilts and focal lengths are exported and imported as explicit per-mirror data,
both codes trace the same telescope. Only the scatter draw and photon sampling
then differ, and with 10^7 photons the agreement should be limited by about
1e-4. This is the single most valuable precision upgrade of the validation.

**A3. Not reproducible from the repository.** `docs/VALIDATION.md` cites commit
hashes, a container (`obdeect-val`), scripts under `val/` and a "validation
workspace (`phase3`)", none of which are in the tree. `AGENTS.md` asks for
provenance on fixtures and cross-tool comparisons. Put the scripts and the
result tables in a versioned `validation/` directory, or publish them as a
DOI-referenced data release (the Zenodo record already exists for the code),
and make the tables in `VALIDATION.md` generated from it.

**A4. As-built effects missing.** Already listed in `limitations.md`:
misalignment of dual-mirror segments and the secondary (SST, SCT), camera-window
refraction (focus shift, section 4), and structure shadowing for LST, SST and
SCT (only the MST has a `telescope_structure` table; others use the empirical
`telescope_transmission`). For absolute effective-area work on the LST,
shadowing by the real masts and camera frame is the largest of these. This is
where obdeect can become more accurate than its reference, as the MST structure
already does.

**A5. Atmospheric-transmission workaround in the shower route.** With bypassed
optics, sim_telarray derives the airmass from the camera-frame direction, so
the route rescales every photon weight by the ratio of two transmissions
(`VALIDATION.md`, "Air-shower photons"). Without it, 7 % (40 deg) and 33-39 %
(60 deg) too many photoelectrons. This is correct but fragile: a change in
sim_telarray's formula or a flag change silently double-counts. Add a regression
test that fails when the product of the two transmissions deviates from the
native chain by more than 0.1 %, and push the sim_telarray option request.

**A6. Wavelength quantisation.** `lambda = floor(exact)` copies sim_telarray
(`bunches.hpp`, `wavelength_law`). It shifts the mean wavelength by about
0.5 nm. Keep it, but make it an explicit option (`wavelength_quantisation="sim_telarray" | "none"`) so the physics-exact setting is available and the parity
setting is a named choice.

**A7. Timing.** `propagation_group_index` is a single scalar
(`model.hpp:75`); there is no wavelength dependence and no transit time in
glass. The documented effect is 5-17 ps for 3-10 mm windows, so this is fine
for current needs. It matters only if sub-100 ps pulse shapes become a
requirement.

## 4. Optics expert

### 4.1 What is strong

- Double precision throughout; sphere intersection uses the numerically
  safe form and rejects back faces explicitly (`surfaces.hpp:79-131`).
- Aspheres evaluated in `r^2` with Horner, a rigorous slope bound
  (`slope_bound`, `surfaces.hpp:285`) used to prove a single root, and work done
  near the vertex plane to avoid precision loss for sources 10^7 m away
  (`surfaces.hpp:347`).
- Grid acceleration structures are cross-checked against brute force in the
  tests, and ties are broken by lowest id so that acceleration never changes the
  answer.
- Rounding-level changes are measured against a replay with explicit flip
  and difference limits (`VALIDATION.md`, "Rounding-level changes"); the limits
  are not loosened to pass.

### 4.2 Findings

**O1. Scatter rotation axes are not orthonormal on curved facets [verified in
code and in the project's own note].** `reflect_at` (`trace.hpp:133-138`) uses
the facet's `u` as first rotation axis and `e2 = cross(n, u)` as second. On a
curved facet, `u` is the tangent at the vertex and is not perpendicular to the
hit normal `n`; `VALIDATION.md` records angles up to 1.9e-2 rad. `rotate()` is
Rodrigues' formula, which requires unit axes, but `|cross(n, u)| = cos(angle)`
is less than 1 (relative error about 2e-4 at the worst facets), and the two
axes are not perpendicular. The effect on the scatter width is about 1e-4
relative to sigma, below the validation noise, which is why it was left alone
("This property is not changed by the round"). Fix: Gram-Schmidt `u` against
`n` before use, which costs one dot product. It changes results at that
level, so run it through the existing replay gate and note it in the changelog.

**O2. `surface_slopes` scatter is parity, not physics [to check].** For
aspheric surfaces the code perturbs the *graph slopes* by half the angle
(`trace.hpp:127-132`). Slope space is not angle space: a slope change `ds`
changes the normal angle by `ds * cos^2(theta)` radially and `ds * cos(theta)`
tangentially, so on a surface tilted by 40 deg (SCT secondary) the effective
reflected-ray width is about 0.6 radially and 0.77 tangentially of the nominal
value. If the SCT/SST scatter parameters were fitted with this convention it
may be self-consistent, but then the parameter's meaning is "sim_telarray
class-2 convention", not "RMS reflected-ray angle", and the documentation
should say so. If it was only matched to sim_telarray, measure the PSF
anisotropy before relying on it for anything but parity.

**O3. Root selection when the single-root proof does not hold
(`surfaces.hpp:456-509`) [to check].** In the monotone case (`single_root`)
the first accepted Newton root is the only root. Otherwise the code returns the
first *accepted* root from the first guess that converges, which need not be the
nearest, and the 64-sample bracketing fallback can miss a double root or a
grazing pair. A non-converging Newton iteration also falls through silently
to the fallback and, if that finds nothing, to a "miss"; `numerical_failure` is
reserved for non-finite residuals. Impact is limited to steep rays on strongly
curved secondaries, but a miss is silent. Add a property test that compares
against a dense bracketing search for random rays on the real SST and SCT
primaries and secondaries, including rays at the largest field angle, and
count how often the non-monotone branch is taken.

**O4. Absolute tolerances.** `kMinDistanceM = 1e-9 m` rejects self-hits;
`kTolerance = 1e-10 m` is the Newton residual; verdict margins are `1e-9 m`.
All are far above double rounding for telescope-sized coordinates, but they are
absolute. A model in other units or a very large dish would need them scaled.
Collect them in one documented header with the reasoning (section 7, M5).

**O5. Not traced, by design.** Refraction (camera windows, lenses, Winston
cones) and polarisation (Fresnel s/p split at the 16-50 deg incidence angles
of the focal plane; the reflectivity table depends on angle but is
unpolarised). Reported `focal_incidence_deg` lets the downstream camera
simulation apply cone acceptance. This is a sound boundary; the risk is the
documented one: `focus_offset` may already contain the window shift, so adding
the window to the trace could count it twice. The audit of each telescope's
`focus_offset` origin that `limitations.md` demands should be tracked as an
issue.

**O6. Scatter statistics.** Two independent Gaussians per hit, uncorrelated
between photons, the same everywhere on the mirror. Real mirrors have spatially
correlated slope errors and a wide scattering halo. For D80 and image
parameters this is fine; for the far PSF tails (D95 and beyond, NSB-like
stray light) it is a model limit worth one sentence in `limitations.md`.

**O7. Facet geometry.** `Facet::intersect` rejects early with a plane test
whose bound is `reach + sag*sqrt(1-dn^2)/|dn|`; correct and tested. The grid
visits all facets for `|dz| < 0.05` and uses the bounding box of the segment
from the lowest to the highest facet height otherwise
(`surfaces.hpp:199-218`). For the telescope fields of view this is efficient;
for a very wide-field design the box grows and cost approaches brute force.
Not a defect, but state the assumption.

## 5. HPC expert

### 5.1 Measured [verified]

Apple M-series, 12 logical cores, Release build, 4 million star photons, best of
three, `outputs=["status"]`, small synthetic telescopes from
`python/tests/telescopes.py` (not LST-sized):

| Model | 1 thread | 2 | 4 | 8 |
| --- | --- | --- | --- | --- |
| single-mirror (ns/photon) | 104 | 53 | 27 | 17 |
| dual-mirror (ns/photon) | 260 | 131 | 67 | 42 |

Scaling is linear up to 4 threads and 6.2x at 8 (the host mixes performance and
efficiency cores). Pool and chunking work as intended.

All 14 columns versus `detected_only=True`, 8 threads, same photons:

| Model | all rows | `detected_only` |
| --- | --- | --- |
| single-mirror | 26 ns/photon | 34 ns/photon |
| dual-mirror | 47 ns/photon | 52 ns/photon |

`detected_only` costs 10-30 % **more** than writing everything (see H1). The
documented per-core numbers (LSTN-01 0.26 us/photon, SSTS-17 0.50 us/photon)
are on larger real models and I did not repeat them.

### 5.2 What is strong

- No allocation, I/O or virtual dispatch in the per-photon path; outputs are
  allocated once, first touched by the worker threads.
- Counter-based RNG: determinism is by construction, not by scheduling.
- Dynamic chunk scheduling (4096 photons) copes with the very different cost of
  photons; `ThreadPool` is owned by the caller, no global state, exceptions
  propagate (`trace.hpp:382-401`).
- GIL released around every long call in the bindings.

### 5.3 Findings

**H1. Serial compaction in `detected_only` [verified by measurement,
`trace.hpp:478-554`].** After the parallel trace, `compact_chunks` runs on the
calling thread, calling `memmove` for every chunk and every one of up to 15
columns. It is bandwidth-bound but single-threaded. Options: compact the
columns in parallel (columns are independent); or have each chunk write its
rows into a per-chunk scratch and copy once to the final offset after a prefix
sum (the same pattern `trace_bunches` already uses); or skip compaction when
the detected fraction exceeds a threshold and return a mask. Expected gain: the
whole 8-34 ns penalty.

**H2. Compiler flags on Linux wheels [to check].** `-ffp-contract=off` is set
for determinism (good), but no `-fno-math-errno`. GCC on Linux sets `errno` for
`sqrt`, `log`, `acos`; the compiler then keeps a branch and a libm call after
each `sqrtsd`. The tracer calls `sqrt` in almost every step. `-fno-math-errno`
and `-fno-trapping-math` do not change results (they remove errno and
exception side effects, not rounding), so determinism is preserved. macOS
libm does not set errno, which is why a local macOS benchmark would not show
it. Measure with `bench_trace` on the manylinux build before and after. Wheels
come from GCC, so this affects what most users run.

**H3. Bunch path repeats work (`bunches.hpp`) [verified by reading].** For a
bunch of wavelength 0, Philox is evaluated for the wavelength in the thinning
pass and again in the tracing pass; the `rows_max` pass walks the input a third
time only to validate and count. `keep_flags` is one byte per candidate photon
and `roulette_draws` is a `std::vector<double>` per chunk grown by
`push_back` in the worker. Because the generator is counter-based, the
roulette draw and the wavelength can be recomputed (one Philox, tens of
ns) instead of stored, removing 8 B per survivor and the per-chunk heap
growth. Measure first; a Philox call is a few percent of a trace.

**H4. Thread count default [to check].** `thread_count(0)` is
`std::thread::hardware_concurrency()` (`core.cpp:47`). In containers, cgroup-
limited CI and batch jobs with CPU affinity this can exceed the allotted cores
and oversubscribe. In the Python layer use `os.process_cpu_count()` (Python
>=3.13, respects affinity) as default, and honour an `OBDEECT_NUM_THREADS`
environment variable for cluster users.

**H5. Two `parallel_for` implementations.** One spawns `std::thread`s on every
call (`trace.hpp:324`, used by `photon_arrays` in the bindings) and one uses
the pool. Both are tested, but thread creation per call is a latent cost for
small blocks. Prefer the pool everywhere and delete the spawning version, or
document why both exist.

**H6. Scalar per-photon pipeline.** Each photon runs the full chain in a
branchy, scalar function that returns a ~200-byte `TraceResult`. SIMD across
photons would help the facet-sphere test and the focal-plane step, but the
asphere Newton iteration and the obscurer traversal diverge, and the
standard-library-only rule excludes `std::experimental::simd` and intrinsics.
Treat this as an experiment after H1-H3: a stage-wise "wavefront" over a block
of 8-16 photons with `#pragma omp simd` on the facet loop, kept bit-identical by
`-ffp-contract=off`. I would not start here; the cheap items above are worth
more.

**H7. No performance guard in CI.** The benchmark tools are built only
manually (`OBDEECT_BUILD_BENCH=OFF`), and wall-clock benchmarks on shared
runners are noise. The project already has the right idea in
`test_asphere_steps` (counts of Newton steps with `OBDEECT_BENCH_COUNTERS`).
Extend to deterministic work counters per real model: Newton steps per photon,
obscurer node visits and exact tests, segment `contains` calls. Fail the test
if a counter rises by more than 2 %. These are exact on any machine.

**H8. Memory scale.** The default output has 14 columns, about 130 bytes per
photon; 10^9 photons would need 130 GB. `trace_many`, `trace_emitter_blocks`
and `outputs=` exist; make streaming blocks the documented default for large
runs and add a guard that warns above a configurable output size.

**H9. Release flags.** CMake `Release` gives `-O3 -DNDEBUG`. LTO is
irrelevant for a header-only core in one translation unit. PGO could give 5-10 %
on the wheel but complicates the build; skip unless the cheap items fall short.
`-march=native` is correctly excluded from wheels.

## 6. Release manager

### 6.1 What is strong

- towncrier fragments, enforced by a workflow and a pre-commit check; changelog
  in Keep-a-Changelog form.
- PyPI trusted publishing (OIDC), TestPyPI by manual dispatch, wheel smoke test
  with the full pytest suite inside cibuildwheel.
- Versioned documentation site that keeps `dev` and release pages together; PR
  builds only build `dev`.
- CI concurrency groups cancel superseded runs; sanitizer and FMA jobs exist.
- Version derived from git tags (setuptools-scm).

### 6.2 Findings

**R1. Four different version numbers [verified].**

| Place | Value |
| --- | --- |
| git tag / `CHANGELOG.md` | v0.7.2 |
| `CMakeLists.txt:2` | `VERSION 0.2.0` |
| `CITATION.cff:4` | `version: 0.1.0`, `date-released: 2026-09-18` |
| `pyproject.toml` `fallback_version`, `__init__.__version__` fallback | `0.1.0`, `0.0.0` |

A user citing from `CITATION.cff` cites a version that does not exist. Use
`SKBUILD_PROJECT_VERSION` for the CMake project version (or drop it), update
`CITATION.cff` as a step of the release (a small script run in the release
branch, or let the Zenodo integration own the metadata), and make the two
Python fallbacks one constant that fails loudly in CI.

**R2. Release process is manual and visibly rushed.** v0.7.0 and v0.7.1 were
released on 2026-11-04 and v0.7.2 on 2026-11-06, with "-rc" branches. A
`release.yml` would remove the guesswork: on a `v*` tag check that the
changelog has an entry for the tag, build and test wheels, publish to TestPyPI
automatically for `-rc` tags, publish to PyPI for final tags, and create the
GitHub Release from the changelog section. Add a documented release checklist
to `docs/` or `AGENTS.md`.

**R3. Platform matrix does not match the claims.** `AGENTS.md` says
"portable". CMake has MSVC branches (`/W4 /WX`), but nothing builds on
Windows. `MACOSX_DEPLOYMENT_TARGET = 13.3` is set, but the only macOS runner
(`macos-15`) is arm64, so no Intel wheel exists. No free-threaded
(`cp314t`) wheel, although the GIL is already released around all long calls.
The sanitizer job and the FMA job run on one OS only. Either add Windows
(MSVC), Linux/Clang and a minimum-compiler job (the oldest GCC and Clang that
the C++20 code needs), or write the supported platforms into the docs and
classifiers. Also consider nanobind's stable-ABI mode: with
`requires-python >= 3.14` one `abi3` wheel per platform would replace the
cp314 and cp315 pair and cover future versions, at a small speed cost to
measure.

**R4. Test CI does the build twice and has no cache.** `tests.yml` runs
`pip install ".[plot]"` (a full wheel build), then configures and builds
the C++ tree again with CMake, for six matrix cells. Add `ccache` (or
`actions/cache` on the build directory), or build once and reuse. Python 3.15
is final as of this writing, so `allow-prereleases: true` can go.

**R5. Missing CI gates.**

| Gap | Why it matters |
| --- | --- |
| ThreadSanitizer job | `ThreadPool` uses a hand-written mutex/condition-variable protocol with a stack-allocated `Job` reachable from workers (`trace.hpp:382-445`); TSan is the right tool |
| Sanitizer run of the Python tests | Bindings own raw buffers via capsules (`core.cpp:37`) |
| `clang-tidy` | `AGENTS.md` asks for it; `.clang-tidy` enables 9 checks; CI never runs it |
| Coverage report | Cannot see untested branches (`classify_polygon`, `tangent_frame` edge cases) |
| Type checking and `py.typed` | The Python layer is annotated but never checked |
| Docs link check | The README dead links below would be caught |
| Fuzzing of `json.hpp` / `model_json.hpp` | Parser of user-supplied files |
| Linux/Clang, MSVC, minimum compiler | Only GCC (Linux) and AppleClang (macOS) are exercised |

**R6. Supply chain.** `changelog.yml` pins `actions/checkout` by commit SHA;
every other workflow uses floating major tags (`@v7`, `@v8`) and
`pypa/gh-action-pypi-publish@release/v1`. There is no Dependabot or Renovate
configuration, and `.pre-commit-config.yaml` pins `ruff-pre-commit` at v0.12.7
while `pyproject.toml` enables `preview = true` and targets py314, so hook and
local ruff versions can drift. Pin all Actions by SHA (keep the tag in a
comment, as `changelog.yml` does), add `.github/dependabot.yml` for
`github-ecosystem` and `pip`, and run `pre-commit autoupdate` through
pre-commit.ci or a scheduled job. The `publish` job is the one place where a
compromised action matters most.

**R7. Packaging metadata.**
- `Development Status :: 2 - Pre-Alpha` while a production-grade validation
  exists; `3 - Alpha` or `4 - Beta` is the honest classifier.
- `license = { file = "LICENSE" }` is the legacy table form; PEP 639 prefers
  `license = "BSD-3-Clause"` with `license-files`.
- No `[project.urls]` (Homepage, Documentation, Issues, Changelog): the PyPI
  page has no links.
- The sdist includes `cpp/**` with tests, benchmarks and fixtures; fine for
  reproducibility, but `docs/`, `CITATION.cff` and `CHANGELOG.md` are absent.
- No `py.typed`, no `Typing :: Typed`.

**R8. CMake as a library.** No `install()` or exported target, no
`obdeect::core` alias, `OBDEECT_BUILD_TESTING` defaults to ON even when
obdeect is a subproject (use `PROJECT_IS_TOP_LEVEL`), and
`find_package(Python 3.14 ...)` hard-codes the version. Standard-library-only
C++ is a selling point for embedding (for example in a CORSIKA or
sim_telarray reader); give it a three-line `find_package` story.

**R9. README links are dead [verified].** `README.md:14` links to
`docs/USER_GUIDE.md` and `docs/PHYSICS.md`; the files are `docs/user-guide.md`
and `docs/optics.md`. As instructed, I did not touch `README.md`. Suggested
change:

```text
See the [user guide](docs/user-guide.md), the [physics](docs/optics.md) of the
photon path, and the [validation](docs/VALIDATION.md) against sim_telarray.
```

Also unify the file naming in `docs/` (`VALIDATION.md` upper-case, all others
lower-case) so that links survive case-sensitive file systems.

## 7. Maintainability and documentation

### 7.1 What is strong

- Comments explain reasoning and bounds (for example the plane pre-test and the
  segment verdict proof) rather than restating code.
- Tests follow the stated order: analytic kernels (162 checks in
  `test_kernels.cpp`), then traced analytic optics, then models, then
  cross-tool; sanitizer-clean in CI.
- The project documents rejected experiments (Estrin evaluation, transmission
  table, closed-form scatter) with numbers. This saves future maintainers days.
- `AGENTS.md` is specific and enforceable.

### 7.2 Findings

**M1. `trace_photon` has 12 parameters, four of them adjacent bools
(`trace.hpp:175-179`).** `apply_responses, direct, incidence_angles,
classify_losses` can be swapped without a compiler diagnostic, and
`BlockTask` and `BlockRequest` repeat the same flags. A
`struct TraceOptions { bool apply_responses = true; ... }` passed by value
costs nothing in the hot path (it inlines) and gives names at call sites. Do it
before adding the next flag.

**M2. Large mixed-purpose headers.** `surfaces.hpp` (931 lines) holds apertures,
facets, the facet grid, aspheres, segments and the segment grid;
`trace.hpp` (666) holds the physics, the thread pool and block I/O;
`obscurers.hpp` (763) holds primitives and the bounding hierarchy. Split into
`aperture.hpp`, `facets.hpp`, `asphere.hpp`, `segments.hpp`, `parallel.hpp`,
keeping `trace.hpp` as the physics. This is mechanical, helps review diffs, and
the build is fast enough (2.7 s) that compile time is not a reason to wait.

**M3. Repeated literals.** `6.283185307179586476925` appears in at least five
places, `0.86602540378443864676` in three, and `kDegree` is local to
`trace.hpp`. Put `kTwoPi`, `kPi`, `kSqrt3Over2`, `kDegree` in `vec.hpp` or a
`constants.hpp`.

**M4. Numeric thresholds are scattered.** `1e-6` (unit direction), `1e-9`,
`1e-10`, `1e-12`, `0.05` (grid dz), `0.5` (`kGridMinDz`), `4096` (chunk), `64`
(samples, depth). Gather the ones that carry physical or proof meaning in one
header with a sentence each, and add a test that perturbs each by a factor of
10 and verifies the results are unchanged (or that the change is flagged).

**M5. JSON number parsing is locale-dependent [verified in code,
`json.hpp:144`].** `std::strtod` honours `LC_NUMERIC`. The grammar check before
it accepts only `.` as separator, but under a process locale such as `de_DE`
`strtod("1.5")` parses `1` and stops. Python does not set `LC_NUMERIC`, but
GUI toolkits, some C extensions and embedding hosts do, and the user base is
European. Use `std::from_chars` (GCC 11+, MSVC 2019+; AppleClang from Xcode 16
supports `double`; check the minimum), or parse integer and fractional parts
with exact scaling. A two-line test with `setlocale(LC_NUMERIC, "de_DE.UTF-8")`
guards it.

**M6. Fixtures are synthetic only.** The Python tests use a small single-mirror
and an SST-like dual-mirror telescope; the native tests use analytic cases.
Real models (LST, MST with 203 structure elements, SST, SCT) exist only in
external simulation-models checkouts, so `bench_trace` and `replay` cannot be
run from a clean clone. Add compiled optical models (the format is stable
`obdeect.optical-model.v2`, hash-addressed) as release assets or a small data
package with provenance, then run replay against stored statistics in CI.

**M7. Documentation structure.** The docs set is complete (installation, user
and Python guides, optics, models, shadowing, light sources, batch tracing,
air-shower photons, limitations, validation). Improvements: a single page
"supported parameters of the simulation models and what obdeect does with each"
(accepted, ignored, rejected), because the design rule says never silently
discard fields and users will ask; a generated C++ API reference (Doxygen or
a hand-kept header map); and `VALIDATION.md` split into a stable method page and
a dated results page so that old commit hashes do not read as current.

**M8. Small items.**
- `ctest` registers one test per executable, so a failure in `test_kernels` is
  one red line for 162 checks; registering sections with `add_test` filters (or
  using doctest/Catch2 only under `tests/`) would make failures local. Optional;
  the custom `check.hpp` is a fair choice for a no-dependency core.
- The `OBDEECT_BENCH_COUNTERS` macro and static counter singletons are fine for
  single-thread tests, but make the dependence explicit in the type
  (`CountedAsphere`) if they spread.
- `numerical_failure` is documented as "did not converge" but, as noted in
  O3, only fires on non-finite values. Align name and behaviour.

## 8. Suggested order of work

**Week 1 (hygiene, low risk, no change in results).**
R1 (versions), R9 (README links, applied by the owner), R6 (Dependabot, SHA
pins), R4 (cache, drop prereleases), R5 first row (TSan job), M5 (locale test
and `from_chars`), the environment note in section 2.

**Week 2-3 (performance without changing results).**
H1 (parallel compaction), H2 (flags, measured on the manylinux build), H3
(bunch pass cleanup), H4 (thread default), H7 (deterministic work counters),
M1 (options struct).

**Month 1-2 (accuracy and trust).**
A2 with A3 (shared mirror realisations, reproducible validation), O1 (scatter
axes, with the replay gate), O3 (root-selection property test), A1 (analytic
aberration references), O2 (decide and document the slope convention).

**Later, only if needed.**
Windows and Intel-macOS wheels, stable-ABI wheels, free-threaded wheels, LST
structure shadowing, SIMD wavefront experiment (H6), PGO.

## 9. Open questions for the maintainers

1. Is `surface_slopes` meant to be the physical RMS reflected-ray angle, or the
   sim_telarray class-2 convention? (O2)
2. Can sim_telarray print or accept per-mirror tilts and focal lengths, so both
   codes trace one realisation? (A2)
3. Are Windows and Intel macOS supported targets, or is "portable" meant at
   source level only? (R3)
4. Where should the validation scripts and results live long term: this
   repository, a sibling repository, or a Zenodo data record? (A3)
5. What was the origin of each telescope's `focus_offset`, window included or
   not? (O5)
6. Is the nearest-root guarantee required for the non-monotone aspheric case, or
   is "an accepted root" enough? (O3)
