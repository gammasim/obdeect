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

### Upward photons

An upward photon (cz > 0) emitted between the observation level and the telescope ($z < h_T$) has the
same formula: the absolute value is the optical depth between emission and telescope, and the
airmass is $1/|c_z|$ with $|c_z| = \sqrt{1 - c_x^2 - c_y^2}$ (the sign of cz does not enter). For
a photon from height $h_e$ to the telescope at $h_t$ in a plane-parallel atmosphere,

$$
T_{\text{atm}} = \exp\!\left[-\,\frac{|\tau(\lambda, h_t) - \tau(\lambda, h_e)|}{|c_z|}\right]
$$

(tested against the table at both heights, `test_upward.cpp`). In the spherical mode the
slant depth is the integral over the heights between $h_e$ and $h_t$ of
$|d\tau/dh|/\cos\theta(h)$, where $\cos\theta_t = u$ is the cosine of the upward direction to the
local vertical at the telescope (the cz column gives it, $u = d \cdot n > 0$). An upward photon moving
down against the local vertical ($u < 0$) gets 0. The spherical slant depth is checked against a
fine integral of the table and, for $R \to \infty$, against the plane-parallel value.

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

## Chord integral

The plane formula $|\tau(h_1) - \tau(h_2)| / |u|$ of a straight path with cosine $u$ to the local
vertical is ill-conditioned for $u \to 0$. For a straight chord the same optical depth is the
integral of the extinction coefficient $k = |d\tau/dh|$ along the path,

$$
\tau_{\text{chord}} = \int_0^L k\big(h(s)\big)\,\mathrm{d}s,
$$

with $L$ the chord length and $h(s)$ its altitude. The table's coefficient is the exact derivative
of its interpolation (`TransmissionTable::extinction_coefficient`): in each interval, with $y =
\log_{10} h$,

$$
k(h) = \frac{|\tau_{i+1} - \tau_i|}{(\log_{10} h_{i+1} - \log_{10} h_i)\, h \ln 10},
$$

zero at and below H2 and above the last height. The linear blend in wavelength blends the
coefficients as it blends the depths. The integral is 8-point Gauss-Legendre on equal pieces, with a
piece count of $\lceil \Delta h / 250\,\mathrm{m} \rceil$, where $\Delta h$ is the altitude change over
the piece: $k$ depends only on the altitude, so a horizontal chord is a single piece. A chord that dips below the observation level (possible only for a refracted row, whose bent path the
straight chord approximates) needs no cut: k = 0 below H2, so its depth is that of the part above
the ground. The chord is split where it crosses a table height (the coefficient has a kink there), and at the lowest point
of a spherical chord (the altitude is monotone between the splits). Any chord length is accepted.

Flat geometry: $h(s)$ is linear, and the integral equals the plane formula exactly
($|\tau(h_1) - \tau(h_2)| / |u|$ with $u = dh/ds$), so the chord is the same function, stable as
$u \to 0$ (for a constant height it is $k(h)\,L$).

Spherical geometry (Earth radius R): for a straight line from the Earth-centred point $v$ with unit
direction $d$, $r(s)^2 = b^2 + (s + c)^2$ with $c = v \cdot d$ and $b = |v \times d|$, and the
altitude $h(s) = r(s) - R$. The line dips to $h = b - R$ at $s = -c$ when this lies inside the
chord. The slant integral of [](curvature.md) is the same integral written in the height.

Near-horizontal use: the chord is used when $|u| < 0.1$ (the constant `near_horizontal_cos`). At
$|u| = 0.1$ the plane (or slant) form and the chord agree to $10^{-9}$ in depth for chords up to
100 km, in both wavelength modes (test: `test_chord`).

## Sphere rows

Rows of the sphere representation (`ground_representation = sphere`, [](propagation.md#sphere-representation))
have the row point $Q = C + r$ on the telescope sphere, the direction $d$ at Q, the emission altitude
$z_{em}$ and the chord length $L$ = `path`. The emission point is $E = Q - L\,d$. For a refracted row
the bent path is approximated by this straight chord. The chord formulas need a unit direction: the
float32 components are normalised in double before E = Q - L d is formed (their norm is 1 only to
about 6e-8, which would move a spherical chord by about 0.3 m in altitude, a depth error of 3e-5
for a horizontal chord at 2 km).

- $h_E = z_{em}$ (the emission altitude as written).
- $h_Q$: the altitude of Q. Flat: observation level + $Q_z$ (sky frame, z above the observation
  level, $Q_z = z_C + r_z$). Spherical: $|Q - \text{Earth centre}| - R$, from the sky position with
  the Earth centre at $(0, 0, -(R + H_2))$.
- $u$: the cosine between d and the local vertical at Q (flat: $c_z$; spherical: $d \cdot n_Q$).

The extinction of a row with $L > 0$ and $z_{em} \ge H_2$ is:

| Case | Flat | Spherical |
|---|---|---|
| $\lvert u \rvert \ge 0.1$ | $\exp\!\big(-\lvert\tau(h_E) - \tau(h_Q)\rvert / \lvert u\rvert - \tau_{\text{focus}}\big)$ | the slant integral between $h_Q$ and $h_E$ with cosine $|u|$ at Q (`transmission_curved`, or `_upward` if $h_E < h_Q$), times $\exp(-\tau_{\text{focus}})$ |
| $\lvert u \rvert < 0.1$ | $\exp\!\big(-\int_E^Q k\,\mathrm{d}s - \tau_{\text{focus}}\big)$ | the same integral along the straight line from E to Q (chord integral) |

$\tau_{\text{focus}}$ is the in-telescope depth of the run (`od_focus`: $\tau(h_T + \ell) - \tau(h_T)$
for the telescope path $\ell$), as for the plane rows. The depth at Q is taken as the table value at
$h_Q$ in place of the telescope altitude; the limits are those of the plane rows: zero below H2, zero
for $\tau(h_Q) > 1000$ and for a total depth above 100. A row with a cosine clearly pointing away from
the ground at Q (spherical) gets $p = 0$, as in the spherical plane mode.

Rows emitted inside the sphere ($L = 0$) have no atmospheric depth: $p = \exp(-\tau_{\text{focus}})$
with the same limits. Rows with a wavelength at or above 1000 nm get $p = 0$.

The Q altitude replaces the telescope altitude in the formulas, so a row point on the sphere
above or below the telescope centre is attenuated over its own path. The sphere rows of a photon
that crosses the telescope centre height (Q at the centre height) agree with the plane rows of the
same photon to $10^{-6}$ relative (test: `test_chord`, flat and spherical).

The sphere rows take the same wavelength assignment, efficiency bounds and thinning as the plane
rows; the efficiency status records the same factors (`atmosphere`, `wavelength`).

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
