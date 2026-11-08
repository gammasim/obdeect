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

Optical depth at height h (metres a.s.l.) and integer wavelength lambda (nm): the table
gives $\tau_i(\lambda)$ at heights $h_i$ above the observation level H2, where $\tau = 0$. Between
the points $\tau$ is linear in $\log_{10} h$; it is 0 at and below H2 and constant above the
last height:

$$
\tau(\lambda, h) = \tau_i + \frac{\log_{10} h - \log_{10} h_i}{\log_{10} h_{i+1} - \log_{10} h_i}\,(\tau_{i+1} - \tau_i).
$$

Wavelengths are truncated to integer nm and clamped to the first and last table row. The
table rows are consecutive integer wavelengths below 1000 nm; the heights must rise above H2.

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
