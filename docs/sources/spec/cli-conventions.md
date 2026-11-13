# CLI conventions

## Process behaviour

- Output: binary stream to stdout, or to `-o FILE` (`--output` where listed). `-` or an empty
  value means stdout in flasher and led. The file is opened (truncated) before the stream is
  planned in the sky tools; on error it may be left empty.
- `--help` (`-h` in flasher and led): usage text to stdout, exit 0, nothing else.
- Error: one line to stderr, exit 1, no usage dump (muon-track also prints its usage to stderr).
  Prefix: flasher and led `<argv0>: <message>`; sky tools `obdeect sources <tool>: <message>`.
- Unknown option: flasher/led `unknown option <arg>`; sky tools `unknown option <arg> (see
  --help)`; muon-track `unknown option <arg>`. Missing value: flasher/led `option <arg> needs a
  value`; sky tools `<arg>: missing value`; muon-track `option <arg> needs a value`.
- A negative `--seed` converts to `uint64` modulo 2^64.
- Options are processed left to right; later repeats overwrite. Value-less flags take no value.

## Numbers

- Real: the whole text must be a decimal/exponent number (C `strtod` grammar); NaN, inf and
  out-of-range values are errors. Integer: optional sign and decimal digits only; trailing text
  and fractions (`1.5`, `1e3`) are errors. Messages: flasher/led `<opt>: not a number:
  '<text>'` / `not an integer`; sky tools `<opt>: not a number: <text>` / `not an integer`;
  muon-track `invalid number for <opt>: '<text>'` / `invalid integer for <opt>`.
- Event counts of sky tools (`to_count`): integer `0..4294967295`; message `<opt>: must be an
  integer in 0..4294967295: <text>`.
- Lists: comma separated, each item a real; an empty item is an error. Vector `x,y,z` needs
  exactly three items (`<opt>: needs three comma-separated values`).
- Table files: two columns, `#` starts a comment, lines that do not parse as two numbers are
  skipped; no rows is an error (`<opt>: no two-column rows in <path>`); unreadable file:
  `<opt>: cannot open table <path>`. Consecutive rows with an equal first column merge into one
  node with the mean second column (so abscissae are strictly increasing).
- `is_number(text)`: `strtod` consumes the whole text (so `4e2`, `-5`, `1.5E-3` are numbers,
  `4e2x`, `1e`, `abc`, empty are not). A text that is a number is never treated as a file name.

## Spectrum syntax

Sky tools (`parse_spectrum`, case-insensitive prefix): a number = fixed wavelength, must be in
`200..2000` nm (message contains `200..2000`); `Cherenkov:a:b` with exactly two fields, both in
`200..2000`, `a <= b` (`a == b` gives Fixed a, else InverseSquare); anything else = table file
(`Tabulated`). Flasher/led (`programs::parse_spectrum`): same numbers and `Cherenkov:` rules
(errors: `spectrum '<text>': ...`), but no table files (other text is an error naming the
supported forms).

## Pulse syntax

Names are case-insensitive prefixes followed by `:` and comma-separated reals; `d` = delay
(default 0). Result per form:

| Text | Pulse |
| --- | --- |
| `Simple:w[,d]` or bare `w` (flasher/led, laser-beam) | `w > 0`: UniformRange(d - w/2, d + w/2); else Fixed(d) |
| `Gauss:s[,d]`, `Gauss-rms:s[,d]` | GaussianShape(mean d, sigma s); flasher/led: `s <= 0` gives Fixed(d) |
| `Gauss-fwhm:f[,d]` | GaussianShape(mean d, sigma f/2.35482004503095) |
| `Exp:tau[,d]` | Exponential(mean d, tau); flasher/led: `tau <= 0` gives Fixed(d) |
| `Complex:a,b,c[,d]` | ComplexPulse(flat a, rms b, decay c, delay d) |
| file name | Tabulated(time ns, value) (sky tools only) |

Sky tools: wrong parameter count is an error (`Simple`, `Gauss*`, `Exp`: 1..2; `Complex`: 3..4);
unknown name is `<opt>: unknown pulse <text>`; a bare number in point-source (`bare_is_gauss`)
is Gauss sigma with mean 0, elsewhere a Simple width (`w <= 0` Fixed 0). Flasher/led: bare text
must start with a digit; `Complex:` takes 1..4 values (missing ones 0) with the first
non-empty; widths `< 0` in `Complex:` are errors (`widths must be >= 0`); other unsupported
text is an error naming the supported forms; no table files.

## Angular syntax

`isotropic`, `Parallel` (= FixedAngle(0)), `Cone:a` (FixedAngle, `u0 = 1 - cos(a deg)`),
`FilledCone:a`, `Gauss:s` / `Gauss-rms:s` (2-D r.m.s. in deg, `gaussian_angular_table`),
`Gauss-fwhm:f` (sigma = f/2.35482004503095), `Rayleigh`, case-insensitive; trailing garbage in
the number (`Cone:30abc`) is an error. Sky tools additionally accept a table file (angle deg,
value; `u = 1 - cos(angle)`); flasher/led do not.

## Fiducial-sphere rule (flasher, led)

`fiducial_target(position, R, angular)` with the sphere centred on the origin: if `|position| <=
1.02 R`: no target (full 4 pi emission). Else Isotropic -> `SphereTarget(origin, R)` (note:
radius `R`, cone from `sin theta = R/d`); FixedAngle -> no target (all photons kept); other
profiles -> `DiscTarget(origin, normal (0,0,1), R)`, and if `position.z == 0` an error
(`a non-isotropic source in the camera plane outside the fiducial sphere is not supported`).
Point-source uses the 1.02 rule with `rmax = fiducial_super_radius` and sphere radius `1.02
rmax` (see its page).

## Telescope list files (point-source)

`--telpos-file`: lines `x y z r` [cm]; `--corsika-telpos-file`: only lines whose first word is
`TELESCOPE`, followed by `x y z r`; other lines skipped. `#` starts a comment, blank lines
skipped. Error (message prefixed `<opt> line N: `): a telescope line without exactly four
words (`needs four values x y z r, got K`), a non-number (also NaN/inf), `r <= 0`; no
telescope records at all: `<opt>: no telescope records [(TELESCOPE lines)] in the file`;
unopenable: `<opt>: cannot open telescope file <path>`. There is never a fallback telescope.
