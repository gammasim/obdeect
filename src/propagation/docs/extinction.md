# Atmospheric extinction

Direct atmospheric transmission of Cherenkov and artificial-source photons, from a table of
optical depth against height and wavelength. Also: wavelength assignment and early thinning
with the efficiency envelope.

## Inputs

| Quantity | Unit | Source |
|---|---|---|
| Transmission table | nm, km a.s.l., optical depth | extinction-table format (`# H2=`, `H1=` header) |
| Site altitude | m a.s.l. | must match H2 within 25 m |
| Telescope altitude | m a.s.l. | telescope position |
| Path inside telescope | m | telescope optics |
| Envelope E(lambda) | nm, [0, 1] | product of later efficiencies, upper bound |

## Formulas

Airmass of a downward photon with direction cosines (cx, cy), plane-parallel:

$$
X = \frac{1}{|c_z|}, \qquad c_z = -\sqrt{1 - c_x^2 - c_y^2}, \qquad X = 10^{10} \text{ if } c_z = 0.
$$

Optical depth at height h (metres a.s.l.) and integer wavelength lambda (nm), linear in
log10 of the height in cm and clamped to the table edges:

$$
\tau(\lambda, h) = \operatorname{rpol}\left(\log_{10} h_{\text{cm}},\; \tau_i(\lambda)\right).
$$

Telescope optical depths, with h_T the telescope altitude and l the path inside the
telescope:

$$
\tau_{\text{tel}} = -\ln T(\lambda, h_T), \qquad
\tau_{\text{focus}} = -\ln\frac{T(\lambda, h_T + l)}{T(\lambda, h_T)},
$$

where T is the vertical transmission $T = e^{-\tau}$. If T(400 nm, h_T) gives
$\tau_{\text{tel}} > 1000$ (the table does not fit the site), both are set to 0.

Acceptance of a photon at emission altitude z with wavelength lambda:

$$
T_{\text{atm}} = \exp\!\left[-\left(\left|\tau(\lambda, z) - \tau_{\text{tel}}\right| X + \tau_{\text{focus}}\right)\right],
$$

with $T_{\text{atm}} = 0$ below the observation level, for $\tau_{\text{tel}} > 1000$, or for an
optical depth above 100. Photons with lambda >= 1000 nm get 0.

Wavelength: for photons without a wavelength, density proportional to $1/\lambda^2$ on
[lo, hi], by inverse CDF, $\lambda = 1 / (1/\lambda_{lo} - u\,(1/\lambda_{lo} - 1/\lambda_{hi}))$
with u uniform in (0, 1). The value is truncated to an integer nm.
The range needs finite $1 \le \lambda_{lo} < \lambda_{hi} \le 2^{24}$ nm. Input wavelengths
must be finite, non-negative integers (0 = not assigned); other values are rejected.

Thinning: the acceptance is $p = T_{\text{atm}}(\lambda)\,E(\lambda)$, where E is linear
between the envelope nodes. Each bunch of weight w gives

$$
n = \mathrm{Binomial}(\lfloor w \rfloor, p) + \mathrm{Bernoulli}\big((w - \lfloor w \rfloor)\,p\big),
$$

independent Bernoulli trials per photon of the bunch. Thinned output has one row per surviving
photon with w = 1. Weighted output keeps every row with $w \leftarrow w\,p$. Both write the
column p.

## Legacy tables

`dlufd-convert-table` (Python: `dlufd.convert_table`) converts older extinction tables to the
format above, so that existing setups can be used:

```console
$ dlufd-convert-table --format kascade us76.ext -o atm_trans_1800.dat --observation-level 1800
$ dlufd-convert-table --format modtran5 us76.ext -o out.dat --observation-level 1800 --min-wavelength 210
```

All legacy layouts hold the cumulative vertical optical depth $c(\lambda, h)$, which rises with
altitude. The optical depth of a vertical path between the observation level and an emission
height is $c(\lambda, h_e) - c(\lambda, h_{\text{obs}})$, with $c$ linear in height and in
wavelength. The converter writes $\tau(\lambda, h) = c(\lambda, h) - c(\lambda, h_{\text{obs}})$ for
every integer wavelength from the first to the last wavelength of the source (up to 999 nm) and
for heights above `--observation-level` (= H2): every altitude of the source, eighth steps in
the lowest 6 km above the observation level and points 10 m to 500 m above it. The absolute
offset of $c$ does not matter. At the grid heights $\tau$ reproduces the source exactly; between
them dlufd interpolates in $\log_{10} h$ instead of $h$, which differs by less than 1 % of
$\tau$. If the grid would exceed the 100 heights a table can hold, coarser steps are used, or the table is rejected. The observation level must lie inside the altitude range of the source and above 0 m.

| `--format` | layout | wavelengths |
|---|---|---|
| `corsika-atmabs` | blocks: a line with one wavelength [nm], then 51 values $c$ for 0, 1, ..., 50 km (any number per line); lines containing `EXTINCTION` are headers | in the file, 5 nm steps (from 180 nm) |
| `kascade` | one row per altitude, 0 km and 1 km steps; one column per wavelength | from 180 nm, 5 nm steps |
| `modtran4` | as `kascade` (`us76.*.ext`) | from 200 nm |
| `modtran5` | as `kascade` | from 205 nm |

`--min-wavelength NM` overrides the first wavelength of the matrix formats. Comments start
with `#`. Files that are not rectangular, contain non-finite or non-numeric values, have
a wavelength step other than 5 nm, or whose optical depth decreases anywhere with altitude above
the observation level are rejected with `ValueError` (command line: message and exit status 1).

## Random numbers

Philox (gstrickd `uniform4`), key (seed, stage `DLUF`, event, telescope, photon index,
stream). Wavelength draws use stream 1, thinning stream 2. The binomial is a sum of
Bernoulli trials, so bunches are limited to 2048 photons (split larger bunches upstream).
Results are reproducible for a given seed; other generators agree in the probabilities, not
in the individual draws.

## Example

```python
import dlufd

table = dlufd.TransmissionTable("atm_trans_2000_1_10_0_0_2000.dat", site_altitude_m=2000.0)
out = dlufd.apply_extinction(
    columns, table, telescope_altitude_m=2030.0,
    envelope_wavelength=[200.0, 800.0], envelope_value=[0.9, 0.3],
    seed=1, event=1, telescope=1, first_index=0,
    weighted=False, wavelength_range=(250.0, 700.0), pathlen_m=28.0,
)
```

## References

- Test reference values: `python/tests/test_extinction.py`, `cpp/tests/test_transmission.cpp`.
- Photon stream and random numbers: gstrickd, `docs/format.md`.
