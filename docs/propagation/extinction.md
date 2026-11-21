# Atmospheric extinction

How much light the atmosphere absorbs and scatters out of the beam on its way to the telescope,
from a table of optical depth against height and wavelength. The page also covers the wavelength
assignment and the thinning of photons with an efficiency bound. The order of all steps is in
[](propagation.md); random numbers are in [](../stream/random-numbers.md).

## Inputs

| Quantity | Unit | Source |
|---|---|---|
| Transmission table | nm, km a.s.l., optical depth | extinction-table format (`# H2=`, `H1=` header) |
| Site altitude | m a.s.l. | must match H2 within 25 m |
| Telescope altitude | m a.s.l. | telescope position |
| Path inside telescope | m | telescope optics |
| Efficiency bound E(lambda) | nm, [0, 1] | upper bound of later efficiencies (product of the `--envelope` tables) |

## Formulas

Airmass of a downward photon with direction cosines (cx, cy), plane-parallel (the default):

$$
X = \frac{1}{|c_z|}, \qquad c_z = -\sqrt{1 - c_x^2 - c_y^2}, \qquad X = 10^{10} \text{ if } c_z = 0.
$$

The plane-parallel airmass is an approximation: its relative error is about $(H/R)\tan^2\theta$,
about 2 % at 75 degrees and above 10 % near 85 degrees. With the `curvature` option the airmass
is replaced by the curved slant optical depth in a spherical atmosphere; see [](curvature.md).
The formulas below for the depth at height h stay the same, since the slant depth uses the same
table.

Optical depth at height h (m a.s.l.) and wavelength lambda (nm). The table gives
$\tau_i(\lambda)$ at heights $h_i$ above the observation level H2, where $\tau = 0$. Between points
$\tau$ is linear in $\log_{10} h$; it is 0 at and below H2 and constant above the last height:

$$
\tau(\lambda, h) = \tau_i + \frac{\log_{10} h - \log_{10} h_i}{\log_{10} h_{i+1} - \log_{10} h_i}\,(\tau_{i+1} - \tau_i).
$$

The table has one row per integer wavelength. Table rows are consecutive integer wavelengths below
1000 nm; there are at most 100 heights including H2, and the H1 heights rise above H2. Two
wavelength modes select the rows (`WavelengthInterpolation`):

- `linear` (default; `--wavelength-interpolation linear`, Python `wavelength_interpolation="linear"`):
  $\tau$ is interpolated linearly in wavelength between rows $n$ and $n+1$, with
  $f = \lambda - n$, after the height interpolation of each row. Every lookup is continuous in
  $\lambda$, including the efficiency bound (below).
- `floor`: the row at $n = \lfloor\lambda\rfloor$. This is the sim_telarray behaviour; the
  reference values of the tests are for this mode, which the tests set explicitly.

$$
\tau(\lambda, h) = (1 - f)\,\tau_n(h) + f\,\tau_{n+1}(h).
$$

For example, at 400.25 nm, $\tau = 0.75\,\tau_{400} + 0.25\,\tau_{401}$. In both modes a wavelength
below the first row uses the first row and one at or above the last row uses the last row (no
extrapolation). At integer wavelengths the two modes are equal.

Telescope optical depths, with h_T the telescope altitude and l the path inside the telescope:

$$
\tau_{\text{tel}} = -\ln T(\lambda, h_T), \qquad
\tau_{\text{focus}} = -\ln\frac{T(\lambda, h_T + l)}{T(\lambda, h_T)},
$$

where $T = e^{-\tau}$ is the vertical transmission. If T(400 nm, h_T) gives $\tau_{\text{tel}} > 1000$
(the table does not fit the site), both are set to 0. They are tabulated per integer nm; in linear
mode they are blended in the same way, $\tau_{\text{tel}}(\lambda) = (1-f)\,\tau_{\text{tel},n} +
f\,\tau_{\text{tel},n+1}$, so $\tau_{\text{tel}}(\lambda) = \tau(\lambda, h_T)$ in both modes. The
curved slant depth (see [](curvature.md)) blends the row slopes the same way, which gives the
blend of the floor depths at $n$ and $n+1$.

Acceptance of a photon at emission altitude z with wavelength lambda:

$$
T_{\text{atm}} = \exp\!\left[-\left(\left|\tau(\lambda, z) - \tau_{\text{tel}}\right| X + \tau_{\text{focus}}\right)\right],
$$

with $T_{\text{atm}} = 0$ below the observation level, for $\tau_{\text{tel}} > 1000$, or for an optical
depth above 100. Photons with lambda ≥ 1000 nm get 0.

Wavelength: for photons without one, density ∝ $1/\lambda^2$ on [lo, hi], by inverse CDF,
$\lambda = 1 / (1/\lambda_{lo} - u\,(1/\lambda_{lo} - 1/\lambda_{hi}))$ with u uniform in (0, 1). The
value is continuous, not rounded. The range needs finite $1 \le \lambda_{lo} < \lambda_{hi} \le 2^{24}$ nm.
Input wavelengths must be finite and non-negative (0 = not assigned); other values are rejected.

Thinning: the acceptance is $p = T_{\text{atm}}(\lambda)\,E(\lambda)$, with E linear between the
nodes of the efficiency bound. Each bunch of weight w gives

$$
n = \mathrm{Binomial}(\lfloor w \rfloor, p) + \mathrm{Bernoulli}\big((w - \lfloor w \rfloor)\,p\big),
$$

independent Bernoulli trials per photon of the bunch. Thinned output has one row per surviving
photon with w = 1. Weighted output keeps every row with $w \leftarrow w\,p$. Both write the column p
only with `--store-p`.

## Python example

```python
import obdeect.propagation

table = obdeect.propagation.TransmissionTable(
    "atm_trans_2000_1_10_0_0_2000.dat", site_altitude_m=2000.0
)
out = obdeect.propagation.apply_extinction(
    columns,
    table,
    telescope_altitude_m=2030.0,
    envelope_wavelength=[200.0, 800.0],
    envelope_value=[0.9, 0.3],
    seed=1,
    event=1,
    telescope=1,
    first_index=0,
    weighted=False,
    wavelength_range=(250.0, 700.0),
    pathlen_m=28.0,
)
```

Older extinction tables are converted with [](legacy-tables.md).
