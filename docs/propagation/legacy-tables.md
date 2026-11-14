# Legacy extinction tables

`obdeect propagation convert-table` (Python: `obdeect.propagation.convert_table`) converts older
extinction tables to the format above:

```console
$ obdeect propagation convert-table --format kascade us76.ext -o atm_trans_1800.dat --observation-level 1800
$ obdeect propagation convert-table --format modtran5 us76.ext -o out.dat --observation-level 1800 --min-wavelength 210
```

All legacy layouts hold the cumulative vertical optical depth $c(\lambda, h)$, which rises with
altitude. The optical depth of a vertical path between the observation level and an emission
height is $c(\lambda, h_e) - c(\lambda, h_{\text{obs}})$, with $c$ linear in height and in
wavelength. The converter writes $\tau(\lambda, h) = c(\lambda, h) - c(\lambda, h_{\text{obs}})$ for
every integer wavelength from the first to the last of the source (up to 999 nm), and for heights
above `--observation-level` (= H2): every altitude of the source, eighth steps in the lowest 6 km
above the observation level, and points 10 m to 500 m above it. The absolute offset of $c$ does not
matter. At the grid heights $\tau$ reproduces the source exactly; between them obdeect propagation
interpolates in $\log_{10} h$ instead of $h$, which differs by less than 1 % of $\tau$. If the grid
would exceed the 100 heights (H2 included), coarser steps are used, or the table is rejected. The
observation level must lie inside the altitude range of the source and above 0 m.

| `--format` | layout | wavelengths |
|---|---|---|
| `corsika-atmabs` | blocks: a line with one wavelength [nm], then 51 values $c$ for 0, 1, ..., 50 km (any number per line); lines containing `EXTINCTION` are headers | in the file, 5 nm steps (from 180 nm) |
| `kascade` | one row per altitude, 0 km and 1 km steps; one column per wavelength | from 180 nm, 5 nm steps |
| `modtran4` | as `kascade` (`us76.*.ext`) | from 200 nm |
| `modtran5` | as `kascade` | from 205 nm |

`--min-wavelength NM` overrides the first wavelength of the matrix formats. Comments start with `#`.
Files that are not rectangular, contain non-finite or non-numeric values, have a wavelength step other
than 5 nm, or whose optical depth decreases with altitude above the observation level are rejected
with `ValueError` (command line: message and exit status 1).

The converted table is the input of [](extinction.md).
