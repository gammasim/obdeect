# Python API

Package `obdeect.sources`: `__version__` (distribution version or `"0.0.0"` if not installed);
extension `obdeect.sources._core` (nanobind) with `core_version()` (three-tuple), `stream_contract()`
(dict: `format_major`, `format_minor`, `max_chunk_rows`, `sources_format_major`, `frames`
{sky 0, telescope 2, camera 3}, `emitted_columns` list in obdeect stream order), `make_spectrum`,
`make_pulse`, `make_angular` (kind, params, x, y) and classes `Emitter(position, pointing,
angular)`, `Star(towards_source, distance_m, disk_radius_m, disk_z_m, start_z_m)`, `Source(frame
int, geometry, spectrum, pulse, photons)` with `set_disc_target`, `set_sphere_target`,
`validate`, `expected_count`, `photon_count(seed, event, telescope)`, `emit(seed, event,
telescope, first_index, rows)` returning a dict of float32 numpy arrays (columns in
`source_columns` order). Core `invalid_argument` becomes `ValueError`; vectors must have
three finite components; frame int must be 0, 2 or 3.

## `obdeect.sources.source`

Constants: `CHUNK_ROWS = 65536`, `FRAMES = {"sky":0,"telescope":2,"camera":3}`, `COLUMNS =
("x","y","z","cx","cy","cz","t","w","wavelength")` (as in obdeect stream).

Source dictionary (unknown keys at any level are `ValueError("unknown <label> key(s): ...")`;
`kind` allowed only in profile and target sections):

- top level: `frame` (default `sky`), `photons` (required), exactly one of `emitter` / `star`,
  `spectrum` (default fixed 400), `pulse` (default fixed 0), `target` (emitters only).
- `emitter`: `position_m`, `pointing`, `angular`, all required.
- `star`: `towards_source` (3 numbers, not all zero; normalised by the API), `distance_m`
  (default inf), `disk_radius_m`, `disk_z_m`, `start_z_m` required.
- spectrum kinds: `fixed{wavelength_nm}`, `uniform{low_nm,high_nm}`,
  `inverse_square{low_nm,high_nm}`, `table{wavelength_nm[],density[]}`.
- pulse kinds: `fixed{time_ns}`, `tophat{width_ns, delay_ns=0}` = UniformRange(delay -
  width/2, delay + width/2), `gauss{sigma_ns, mean_ns=0}`, `exponential{decay_ns, mean_ns=0}`,
  `complex{flat_ns,gauss_rms_ns,decay_ns,delay_ns=0}`, `table{time_ns[],rate[]}`.
- angular kinds: `isotropic`, `cone{angle_deg}`, `filled_cone{angle_deg}`, `lambertian`,
  `gauss{sigma_deg}`, `rayleigh`, `table{angle_deg[],value[]}` (angles finite in `[0,180]`,
  converted to `u = 1 - cos`). Cone angles must be finite in `[0,180]` (message contains
  `finite and in`).
- target: `disc{centre_m, normal, radius_m}`, `sphere{centre_m, radius_m}`.
- A missing parameter is `ValueError("missing '<name>'")`; missing `kind`: `the <label> needs a
  'kind'`; unknown kind: `unknown <label> kind '<k>'`; both/neither emitter-star: message
  contains `exactly one`; target on a star: `only emitters take a target`; fixed-angle with disc
  target message contains `fixed-angle`; zero `towards_source`: message contains `zero`; frame
  `ground`: message contains `frame`.

Functions: `build_source(spec)` (validates), `expected_photons(spec)`, `photon_count(spec, *,
seed, event=0, telescope=0)`, `emit_chunk(spec, *, seed, event=0, telescope=0, first_index=0,
rows=65536)` (rows in `1..65536`, message contains `rows`; `first_index >= 0`, message
contains `first_index`), `iter_chunks(spec, *, seed, event, telescope)` (chunks of at most
65536, `ceil(count/65536)` chunks), `emit(...)` (concatenation; empty float32 arrays when the
count is 0). Semantics identical to [source-core](source-core.md) (the dictionary maps one to
one; `photons` is the unnormalised source mean).

## `obdeect.sources.plot` (optional matplotlib)

`profile_figure(spec, *, seed=1, photons=200000)` -> figure with three axes whose x labels
start with `angle`, `wavelength`, `emission time`; histograms (60 bins) of the first
`min(photons, count, 65536)` photons weighted by `w`; angle between the direction and the
emitter pointing (or `-towards_source` of a star). Zero photons: `ValueError("the source emits
no photons")`.
