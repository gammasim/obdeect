# obdeect sources muon-track

Cherenkov light of one muon-like track, filtered to a telescope sphere at the origin, sky
frame. Writer C, producer `obdeect sources muon-track`; source in [track](track.md). Namespace
`obdeect::sources::muon_track`.

## Options (all take a value unless stated)

| Option (aliases) | Default | Rules |
| --- | --- | --- |
| `-r`, `--radius`, `--telescope-radius` | 15 | sphere radius [m], `> 0` finite |
| `-x`, `-y` | 0 | impact position [m], finite |
| `-h`, `--h-start` | 3000 | start altitude above the observation level [m], `> 0` |
| `--theta`, `--zenith-angle` | 0 | `0 <= theta < 90` deg |
| `--phi`, `--orientation`, `--azimuth` | 0 | deg, finite |
| `-n`, `--n1`, `--refraction` | `0.000292` (n-1) | value `v`: if `1 < v < 1.0003` then `n-1 = v-1`; else if `0 < v < 0.0003` then `n-1 = v`; else error `invalid refraction index: <v>`; final `n-1` in `(0, 0.0003)` |
| `-A`, `--atmosphere` | | only `const` or `-1` accepted, else error (atmosphere tables not supported) |
| `--index-profile FILE` | none | altitude [m] and `n-1` pairs |
| `-a`, `--altitude`, `--obslev`, `--observation-level` | 2150 | written to run record |
| `-b`, `--bs`, `--bunch-size` | 5 | `> 0`; no effect on photons |
| `-l`, `--wlmin` / `-L`, `--wlmax` | 240 / 1000 | `200 <= wlmin <= wlmax <= 2000` nm |
| `-q`, `--charge` | -1 | `!= 0` |
| `-m`, `--mass` | 0.105658 | GeV, in `[0.0005, 100]` |
| `-p`, `--momentum` | | GeV, `> 0`; sets `beta = p / sqrt(p^2 + m^2)` with the final mass, after all options (order independent) |
| `--beta` | 1 | `0 < beta <= 1`; with `-p`: error `--beta and --momentum are exclusive` (either order) |
| `--primary ID\|NAME` | 6 | a text starting with a digit is the CORSIKA id, else a name; also sets `mass` and `charge` via `mass_charge` at that point (later `-m`/`-q` override) |
| `-N`, `--events` | 1 | integer `1..4294967295` (`events must be between 1 and 4294967295`) |
| `--run` | 1 | run number |
| `--seed` | 1 | |
| `-o`, `--output` | stdout | |
| `--enable-scatter` (`--scatter`), `--disable-scatter` (`--no-scatter`) | off | flags without value; enabled -> error `multiple scattering is not supported ...` at validation |
| `--rscat` | 0 | must be 0, else `--rscat is not supported` |
| `--nscat` | 1 | must be 1, else `--nscat must be 1 (each event is used once)` |
| `-I...` | | any option starting `-I`: `-I (atmosphere table search path) is not supported` |
| `--help` | | usage to stdout, exit 0 |

An option other than the flags or `--help` without a following value: `option <arg> needs a
value`. Unknown option: `unknown option <arg>`. `parse_command_line(args)` (vector of strings
without argv[0]) returns `{options, output, help}` and throws `invalid_argument`.

Names: `gamma` 1, `positron` 2, `electron` 3, `muon`/`mu-` 6, `mu+` 5, `proton` 14; other name:
`unknown primary name: <name>`. `mass_charge(id)` returns mass [GeV] and charge: id 3: -1,
0.000511; 2: +1, 0.000511; 6: -1, 0.105658; 5: +1, 0.105658; 8: +1, 0.13957; 11: +1, 0.493677;
14: +1, 0.938272; 13 and 25: 0, 0.939565; 9: -1, 0.13957; 12: -1, 0.493677; 10 and 16: 0,
0.497648; 15: -1, 0.938272; 66-69, 133, 134: 0, 1e-11; 130: +1, 1.777; 131: -1, 1.777; 7: 0,
0.134; 17: 0, 0.54751; 50: 0, 0.78265; 51: 0, 0.7755; 52: +1, 0.7755; 53: -1, 0.7755; other ids
in `(1,100)` or `(110,199)`: 0, 1.2; ids `100..110` and `>= 199` (nuclei): charge `id mod 100`, mass
`(0.931494 - 0.000511) * floor(id/100)` (id 402: `0.931845 * 4`); ids `<= 1`: charge 0, mass 0.

Further validation (`validate(Options)`): finite impact position; the `beta` and `mass`
ranges above; non-empty profile needs equal sizes, at least two nodes, strictly increasing
altitudes and positive `n-1` values; `bunch > 0`. Messages are free text; the type is
`invalid_argument`.

Index profile file: whitespace separated numbers; non-numeric content: `index profile <name>:
expected 'altitude n-1' lines`; no data: `has no data`; odd count: `a value is missing its pair`.

## Geometry

`theta, phi` in rad: `direction = (-sin(theta) cos(phi), -sin(theta) sin(phi), -cos(theta))`.
`dist = h_start / cos(theta)`. Impact point `I = (x, y, 0)`.

Constant index (no profile): `n = 1 + n1`, `cos_ct = 1/(n beta)`. If `cos_ct > 1`: no track
segments (events contain empty telescope blocks). Otherwise a single `Track` with
`start = I - direction dist`, `length = dist`, `time = -dist/(beta c)`, `cos_theta = cos_ct`.

Profile: `steps = ceil(dist / 10)`, `step = dist/steps`. Walking from the impact point up
(`lower = I`, `travelled = 0`): `upper = lower - direction step`, `mid = lower - direction
step/2`, `travelled += step`, `n = 1 + index_at(mid.z)`, `cos_ct = 1/(n beta)`; if
`cos_ct <= 1` emit a Track (`start = upper`, `length = step`, `time = -travelled/(beta c)`,
`cos_theta = cos_ct`); `lower = upper`. Segments are in this order (lowest first).
`index_at(h)`: linear interpolation of `n-1` between nodes, constant below the first and above
the last node.

Every Track: `beta`, spectrum = Fixed(wlmin) if `wlmin == wlmax` else InverseSquare(wlmin,
wlmax), target `SphereTarget((0,0,0), radius)`, and mean photons (Frank-Tamm)
`photons = 2 pi alpha q^2 (1 - cos_ct^2) (1/wlmin - 1/wlmax) 1e9 * length`, `alpha =
0.0072973526`, `q` the charge, wavelengths in nm. Default: one segment, `photons = per_m *
3000`.

## Stream

Writer C: run `wavelength_min/max` = wlmin/wlmax, telescope id 1 at (0,0,0), radius `r`. Events
`1..N` (64-bit loop counter). Event record: `primary_id`, `energy_gev` (= `gamma * mass` with
`gamma = 1/sqrt(1 - beta^2)` if `beta < 1` and `mass > 0`, else `1e-9`), `core_x_m = [x]`,
`core_y_m = [y]`. In each event, for segment `k` in order: `emit_track(segment,
source_seed(seed, k), event, 1, ...)`; rows are appended in order (frame sky). Total weight
sum ~ `expected_track_weight` summed over segments times `N`.

After the stream, stderr gets one line: `obdeect sources muon-track: expected weight <%.6g>, wrote
<rows> photons, weight <%.6g>` (expected = sum of `expected_track_weight` over segments times
events). On error: `obdeect sources muon-track: <message>` plus the usage text on stderr, exit 1.

## Acceptance (`test_muon-track`)

Arc kernel versus brute force; on-axis arcs; default: one segment with the exact Frank-Tamm
photons and `time = -3000/c`, all weights in `[0,1]`, all `t <= 0`, weight sum within 4 sigma
of the expectation; non-default (theta 20, phi 30, impact (40,-30), h 2000, band 300-600, n-1
0.00028): length `2000/cos 20deg`, track ends at the impact point, emission angle exact,
1/lambda^2 moments, times between `-dist/c` and 0; profile: 300 segments of 10 m, first uses
`n-1` at 5 m altitude, `time = -10/c`; scatter and `--rscat` rejected; stream round trip (2
events, `refraction_applied = 0`, `r_m = 15`, sky frame chunks, `telescope_end` sums); `zem -
z = obslev`; magnitude limits above 1e30 and counts above 2^53 rejected; `event_count`,
primary names (`mu+` charge +1, `muon` charge -1, `tau` rejected); index profile reader errors;
momentum option order; `-p 0` rejected.
