# muon-track: specification

Options and examples: [muon-track](../programs/muon-track.md). Source: [track](track.md).
Writer C, producer `obdeect sources muon-track`. Namespace `obdeect::sources::muon_track`.

## Option parsing

- Aliases and values are in the [program options](../programs/muon-track.md#options). An option
  without a value: `option <arg> needs a value`. Unknown option: `unknown option <arg>`.
- `parse_command_line(args)` returns `{options, output, help}` and throws `invalid_argument`.
- Errors on stderr as `obdeect sources muon-track: <message>` plus usage, exit 1.

## Parameter rules

- Theta: `0 <= theta < 90` deg. Phi: finite. Impact `x, y`: finite. Start altitude `> 0`.
- Refraction `v`: if `1 < v < 1.0003`, `n-1 = v-1`; if `0 < v < 0.0003`, `n-1 = v`; else
  `invalid refraction index: <v>`. Final `n-1` must lie in `(0, 0.0003)`.
- Atmosphere: only `const` or `-1`; atmosphere tables are not supported.
- Wavelengths: `200 <= wlmin <= wlmax <= 2000` nm. Charge `!= 0`. Mass in `[0.0005, 100]` GeV.
- Momentum `-p` `> 0`; `beta = p / sqrt(p^2 + m^2)` with the final mass, applied after all
  options. `--beta` in `(0, 1]`; `--beta` with `-p`: `--beta and --momentum are exclusive`.
- `--primary ID|NAME`: a text starting with a digit is a CORSIKA id, else a name. Sets mass and
  charge at that point; later `-m` or `-q` override. Names: `gamma` 1, `positron` 2,
  `electron` 3, `muon`/`mu-` 6, `mu+` 5, `proton` 14; other names: `unknown primary name: <name>`.
- `mass_charge(id)` (mass [GeV], charge): 3: -1, 0.000511; 2: +1, 0.000511; 6: -1, 0.105658;
  5: +1, 0.105658; 8: +1, 0.13957; 11: +1, 0.493677; 14: +1, 0.938272; 13 and 25: 0, 0.939565;
  9: -1, 0.13957; 12: -1, 0.493677; 10 and 16: 0, 0.497648; 15: -1, 0.938272; 66-69, 133, 134:
  0, 1e-11; 130: +1, 1.777; 131: -1, 1.777; 7: 0, 0.134; 17: 0, 0.54751; 50: 0, 0.78265;
  51: 0, 0.7755; 52: +1, 0.7755; 53: -1, 0.7755; other ids in `(1,100)` or `(110,199)`: 0, 1.2;
  ids `100..110` and `>= 199` (nuclei): charge `id mod 100`, mass
  `(0.931494 - 0.000511) floor(id/100)` (id 402: `0.931845 x 4`); ids `<= 1`: 0, 0.
- Scatter (`--enable-scatter`) and `-I...`: rejected. `--rscat` must be 0; `--nscat` must be 1.
- `validate(Options)`: finite impact position; ranges above; a non-empty profile needs equal
  sizes, at least two nodes, strictly increasing altitudes and positive `n-1`; bunch `> 0`.
- Index profile file: whitespace-separated numbers. Errors: `index profile <name>: expected
  'altitude n-1' lines`, `has no data`, `a value is missing its pair`.

## Geometry

- `direction = (-sin(theta) cos(phi), -sin(theta) sin(phi), -cos(theta))`, `dist = h_start / cos(theta)`.
- Impact point `I = (x, y, 0)`.
- Constant index: `n = 1 + n1`, `cos_ct = 1/(n beta)`. If `cos_ct > 1`, there are no segments
  (events contain empty telescope blocks). Otherwise one `Track` with
  `start = I - direction dist`, `length = dist`, `time = -dist/(beta c)`, `cos_theta = cos_ct`.
- Profile: `steps = ceil(dist/10)`, `step = dist/steps`. Walking up from `I`: `upper = lower -
  direction step`, `mid = lower - direction step/2`, `travelled += step`, `n = 1 + index_at(mid.z)`,
  `cos_ct = 1/(n beta)`; if `cos_ct <= 1` emit a Track (`start = upper`, `length = step`,
  `time = -travelled/(beta c)`). Segments are in order, lowest first.
- `index_at(h)`: linear interpolation of `n-1`, constant beyond the first and last node.
- Every Track: spectrum Fixed(wlmin) if `wlmin == wlmax`, else InverseSquare(wlmin, wlmax);
  target `SphereTarget((0,0,0), radius)`; mean photons (Frank-Tamm)
  `2 pi alpha q^2 (1 - cos_ct^2) (1/wlmin - 1/wlmax) 1e9 x length`, `alpha = 0.0072973526`,
  wavelengths in nm.

## Stream

- Run: `wavelength_min/max` = wlmin/wlmax, telescope id 1 at (0, 0, 0), radius `r`.
- Events `1..N`. Event record: `primary_id`, `energy_gev` (`gamma x mass`, with
  `gamma = 1/sqrt(1 - beta^2)` if `beta < 1` and `mass > 0`, else `1e-9`), `core_x_m = [x]`,
  `core_y_m = [y]`.
- In each event, segment `k` in order: `emit_track(segment, source_seed(seed, k), event, 1, ...)`.
- Stderr after the stream: `obdeect sources muon-track: expected weight <%.6g>, wrote <rows>
  photons, weight <%.6g>` (expected = sum of `expected_track_weight` over segments times events).
