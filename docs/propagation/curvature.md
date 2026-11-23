# Earth curvature

The plane-parallel airmass $1/|c_z|$ used by default ([](extinction.md)) assumes a flat atmosphere.
For inclined photons this error is small, but it grows quickly towards the horizon. This page
describes the optional spherical atmosphere of the [atmosphere stage](propagation.md): the
curvature modes, the slant optical depth that replaces the airmass, and how the modes relate to
CORSIKA's CURVED option.

## Why

The plane-parallel airmass overestimates the path through the air by a relative amount of about
$(H/R)\tan^2\theta$, with scale height $H \approx 8$ km, Earth radius $R$ and zenith angle $\theta$:

| zenith angle | relative error of the plane-parallel airmass |
|---|---|
| 45 deg | about 0.1 % |
| 60 deg | about 0.4 % |
| 75 deg | about 2 % |
| 85 deg | more than 10 % |

The plane-parallel airmass diverges at the horizon. The curved slant depth stays finite.

## Modes

The option `curvature` (CLI `--curvature none|corsika|spherical`) selects the mode. The optional
`earth_radius_m` (CLI `--earth-radius-m`) sets the radius, by default 6371315.0 m: the radius of
CORSIKA's CURVED option, constant `obdeect::stream::earth_radius_m`.

| Mode | Input requirement (run `curvature_applied`) | What the stage does | Frames |
|---|---|---|---|
| `none` (default) | 0, or absent (1 is rejected) | Plane-parallel airmass and straight-line geometry | sky → G |
| `corsika` | 1 (required) | CORSIKA has already applied curvature to emission geometry, positions, directions, timing and refraction. The stage only replaces the airmass by the curved slant depth. | G only |
| `spherical` | 0 (required) | The stage applies curvature: emission altitude from the sphere, per-telescope local vertical and telescope altitude, curved extinction, and curved refraction (bending and air delay, see [](#refraction-in-the-spherical-mode)) when an atmospheric profile is given. With `refraction_applied` = 1 in the input only the air delay is applied. The output run record gets `curvature_applied` = 1. Photons already in frame G get the curved extinction only. | sky or G |

The run flag `curvature_applied` is read from the EventIO input: the EventIO importer sets it to 1
when the CURVED flag of the CORSIKA event header (EVTH word 79) is 1 or 2 (CORSIKA documents 0 = standard, 1 = curved with downgoing primaries, and 2 = curved with upgoing primaries also allowed), and to 0 otherwise. The same information is present in bit `0x40` of the
interface option word (EVTH word 77), which sim_telarray reads; the importer uses word 79, and the
two indicators agree in a CORSIKA run. Producers that are not CORSIKA, such as the sources package,
write 0. The stage requires the flag in `corsika` mode, requires 0 in `spherical` mode (the key
must be present in both) and rejects 1 in mode `none`, so curvature is neither applied twice nor
forgotten.

## Slant depth

Photons travel on straight lines. For a telescope at altitude $h_t$ and a photon with zenith
angle $\theta_t$ at the telescope, the invariant of the straight line in a spherical atmosphere is

$$
p = (R + h_t)\sin\theta_t, \qquad
\cos\theta(h) = \sqrt{1 - \frac{p^2}{(R+h)^2}},
$$

and the slant optical depth from the telescope up to the emission altitude $z_{em}$ is

$$
\tau_{\text{slant}} = \int_{h_t}^{z_{em}} \frac{1}{\cos\theta(h)}\,\frac{d\tau}{dh}\,dh \qquad (z_{em} \ge h_t,\ \text{downward}).
$$

For a downward photon below the telescope ($z_{em} < h_t$) the slant depth is zero (the transmission is
$\exp(-\tau_{\text{focus}})$). An upward photon ($z_{em} < h_t$, cz > 0) crosses the slab between its
emission and the telescope: `slant_optical_depth_upward` integrates $|d\tau/dh|/\cos\theta(h)$ over
$[z_{em}, h_t]$, with the same $\cos\theta(h)$ (the straight line with the cosine $\cos\theta_t$ of the
upward direction at the telescope), and `transmission_curved_upward` uses it. The first segment is at
the lower end of that integral, where $\cos\theta$ is smallest.

The table is piecewise: $\tau$ is linear in $\log_{10} h$ between the table heights (see
[](extinction.md)), so on segment $i$

$$
\frac{d\tau}{dh} = \frac{a_i}{h \ln 10}, \qquad
a_i = \frac{\tau_{i+1} - \tau_i}{\log_{10} h_{i+1} - \log_{10} h_i}.
$$

- Each segment, clipped to $[h_t, z_{em}]$, is integrated exactly in the table form above, with
  8-point Gauss-Legendre quadrature per segment.
- For $\cos\theta_t < 0.2$ the integrand behaves like $(h - h_t)^{-1/2}$ at the telescope (downward
  photons; for upward photons at the emission level, where $\cos\theta$ is smallest). The first
  segment at the lower end then uses $h = h_{\text{low}} + u^2$, which removes the singularity.
- Above the last table height $d\tau/dh = 0$; there is no contribution.

The transmission of a photon is

$$
T = \exp\!\left[-\left(\tau_{\text{slant}} + \tau_{\text{focus}}\right)\right],
$$

with $\tau_{\text{focus}}$ the optical depth over the path inside the telescope, as in the
plane-parallel case. It stays plane-parallel: the path is a few metres inside the telescope.

The zenith cosine $\cos\theta_t$ is taken in the telescope frame from the photon direction $\mathbf{d}$
and the telescope's local vertical $\mathbf{n}$:

$$
\cos\theta_t = -\,\mathbf{d}\cdot\mathbf{n} \quad\text{(downward)}, \qquad
\cos\theta_t = +\,\mathbf{d}\cdot\mathbf{n} \quad\text{(upward, } \mathbf{d} = (c_x, c_y, c_z), c_z > 0\text{)}.
$$

Photons moving away from the ground ($\cos\theta_t$ clearly negative) have $T = 0$; an upward
photon with $\mathbf{d}\cdot\mathbf{n} < 0$ is moving down against the local vertical and cannot
arrive from below.

In the plane-parallel limit $R \to \infty$, $\cos\theta(h) \to \cos\theta_t$ and
$\tau_{\text{slant}} \to |\tau(z_{em}) - \tau(h_t)|/\cos\theta_t$, the plane-parallel result.

### Emission altitude in spherical mode

In `spherical` mode the emission altitude is taken from the sphere. With the Earth centre
$C = (0, 0, -(R + h_{obs}))$ in the sky frame (z above the observation level, observation level
$h_{obs}$),

$$
z_{em} = |\mathbf{p} - C| - R,
$$

where $\mathbf{p}$ is the emission point. The local vertical $\mathbf{n}$ of each telescope is the
unit radial vector at the telescope, and $h_t$ its altitude a.s.l.

The mode therefore assumes the geometry of the sky frame: $x = y = 0$ is the array reference, $z$
is the Cartesian height above the observation level (the plane tangent to the sphere at the array
reference), and the Earth centre lies on the vertical through the array reference at
$z = -(R + h_{obs})$. A producer of photons in the sky frame must use this frame.

## Refraction in the spherical mode

In `spherical` mode the ray is bent in spherical shells and the air delay is taken along the bent
ray. Without an atmospheric profile neither is applied. The ray and the straight line are compared
between the same two levels: the telescope level $r_t = R + h_t$ and the emission level
$r_e = R + h_e$, with $r = R + h$ and $n(r) = 1 + \nu(h)$.

A photon is emitted at altitude $h_e$ with a given direction. Its vacuum straight line passes the
telescope level at zenith angle $\theta_s$, with impact parameter $b$ about the Earth centre $C$:

$$
\sin\theta_s = \frac{b}{r_t}, \qquad K = n_e\, b, \qquad n_e = 1 + \nu(h_e).
$$

Along the vacuum line $\sin\theta_{vac} = b/r$. The bent ray follows Bouguer's law
$n(r)\, r \sin\theta_{ref} = K$, with $\theta$ measured from the local radial direction. The
integrals over $r$ from $r_t$ to $r_e$ are

$$
\begin{aligned}
\text{shift} &= r_t \int_{r_t}^{r_e} \frac{\tan\theta_{ref} - \tan\theta_{vac}}{r}\, dr, \\
\Delta t_{\text{bent}} &= \frac{1}{c} \int_{r_t}^{r_e} \left(\frac{n}{\cos\theta_{ref}}
  - \frac{1}{\cos\theta_{vac}}\right) dr, \\
\Delta t_{\text{straight}} &= \frac{1}{c} \int_{r_t}^{r_e} \frac{\nu}{\cos\theta_{vac}}\, dr,
\end{aligned}
$$

with $\sin\theta_{ref} = K/(n r)$ and $\cos\theta_{vac} = \sqrt{(r^2 - r_t^2 + r_t^2\cos^2\theta_s)}/r$.

- `shift` is the displacement of the arrival point along the horizontal direction of motion. It is
  zero or negative: the ray ends nearer the point below the emission.
- $\Delta t_{\text{bent}}$ is the travel time of the bent ray minus that of the vacuum straight line
  between the same levels. It can be negative at large zenith angles, where the bent ray is steeper
  and its shorter path outweighs the slower speed in air. This is the same sign convention as the
  plane-parallel `bend` ([](refraction.md)).
- $\Delta t_{\text{straight}}$ is the air delay along the straight line, used when the input has
  `refraction_applied` = 1.

The arrival direction needs no integral: $\sin\theta_t'/\sin\theta_s = n_e/n_t$ at the telescope.
The horizontal direction cosines $(c_x, c_y)$ are scaled by $f = (1 + \nu_e)/(1 + \nu_t)$, as in the
plane-parallel case.

Because $n$ decreases with height, there is no total reflection for emission from above. For
emission heights above the profile top the integrals take their values at the top, so they are
constant.

### Tables

The integrals are evaluated by quadrature and stored in tables over $(\theta_s, h_e)$, one per level
of a 10 m altitude grid. Each table is built once for its level and looked up per photon with a
$4 \times 4$ Lagrange stencil, with no allocation.

- Quadrature: 8-point Gauss-Legendre on sub-intervals of at most 500 m in $r$. Near the telescope,
  for $\cos\theta_s < 0.2$, the substitution $r = r_t + u^2$ removes the square-root singularity of
  $1/\cos\theta_{vac}$ for grazing rays.
- Zenith grid: steps of 0.25 deg up to 80 deg, then geometric steps in $\cos\theta_s$ down to
  $10^{-4}$, then 90 deg (698 nodes). The interpolation is in $\ln(\cos\theta_s + 10^{-4})$.
- Height grid: geometric steps of 5 % up to 20 km above the telescope, then 1 km steps up to
  $\min(\text{top}, 120\,\text{km})$ (350 nodes).
- Telescope altitude: a telescope at altitude $h_t$ between two grid levels uses the two tables
  that bracket it, and the integrals are interpolated linearly in $h_t$ between them. Spherical
  telescope altitudes differ from the level of the telescope plane by about lateral$^2/(2R)$; the
  second-order error of the linear interpolation over the 10 m spacing is about $3\times10^{-6}$
  relative in the grazing-ray delay and about $10^{-3}$ m in position at 80 deg. At 85 deg with
  emission at 20 km it is 0.02 to 0.03 m and 0.06 to 0.12 ns on a delay of $3.6\times10^4$ ns.
- Build time is about 3 to 5 s per table. Tables are built on first use, when a telescope receives
  photons in the sky frame, and cached by level; a telescope altitude uses two of them.
- CORSIKA's CURVED build also tabulates integrated rays, on a 0.25 deg by 0.25 km grid. The grid
  here is finer than that below 20 km above the telescope and near the horizon, and coarser above
  20 km. The integrals are evaluated directly, not by ray tracing.

### Accuracy and limits

- Quadrature against an independent composite Simpson reference: agreement within the $10^{-6}$
  relative tolerance asserted by the tests.
- Table against direct quadrature: about $10^{-6}$ relative. The tests check 2000 random points over
  all zenith angles and emission heights: `shift` within the larger of 0.01 m and $10^{-3}$
  relative, and each time within 0.01 ns plus $10^{-3}$ relative. Near the horizon, with emission
  close to the telescope, the table is the weakest part of the approximation; the tolerance above
  applies there too.
- Flat limit ($R \to \infty$) against `bend()` and `straight_path_air_delay_ns()`: agreement to
  $10^{-4}$ relative, or 1 mm and $10^{-3}$ ns.
- Physical effect: at zenith 80 deg and emission at 20 km the curved shift is 7 % smaller than the
  plane-parallel one (-234 m against -252 m, ratio 0.93).
- Output precision: the time column of the stream is float32. Its spacing near $3.8\times10^5$ ns is
  0.031 ns (rounding error up to 0.016 ns), so very long grazing-ray delays are limited by the
  resolution of the column.
- A table is referenced to its own level, so it cannot serve a photon on its own when the
  telescope altitude lies between two levels. Where the linear blend is not possible (a photon
  emitted less than one 10 m step above the telescope, a grazing line that does not reach the lower
  level, or a telescope within 10 m of the profile top where no upper level fits inside the
  profile) the integrals are computed directly at the telescope altitude by the reference
  quadrature. This is exact and rare. A photon is dropped (counted in `dropped`) only if the line
  does not reach the telescope level or the emission is below it.
- A table level within a fraction of a metre of the profile top has fewer natural height nodes
  than the interpolation stencil needs; its height grid is then filled up uniformly to eight nodes.
- The shifted arrival point is placed on the telescope-level sphere, with the time referenced to
  that sphere, not on the tangent plane of the frame G. For a shift $s$ the sphere lies
  $s^2/(2R)$ below the plane, so a ray arriving at zenith angle $\theta$ crosses the plane about
  $(s^2/2R)\tan\theta$ away, a time difference of about that distance divided by $c$. The effect
  is not applied; it is below 0.05 m and 0.15 ns at 80 deg and about 2 m and 7 ns at 85 deg with
  emission at 20 km (more at higher emission altitudes), and negligible below 75 deg.
- The zenith cosine of a photon is taken from the stored $c_x$, $c_y$ in double precision (a
  stored $c_z$ is only checked for consistency). With float32 columns the relative error of
  $\cos\theta$ is about $1.4\times10^{-4}$ at 89.5 deg, which is about $7\times10^{-4}$ in the
  transmission there. CORSIKA input has the same float32 limit.
- The telescope plane is the Cartesian plane $z$ = telescope $z$. The local vertical tilt of
  laterally displaced telescopes (at most lateral$/R$) is neglected in the direction of the
  displacement.

## Spherical integrator (sphere representation)

With `--curvature spherical` and the sphere representation ([](propagation.md#sphere-representation))
the ray integrator of [](refraction.md#integrator) uses the spherical geometry: the altitude of a
point is |X - C_E| - R with the Earth centre C_E = (0, 0, -(R + observation level)) and the
gradient of n points along the radial unit vector from C_E. The straight-line air delay of rows
without refraction takes the altitude along the line from |r| (`straight_line_air_delay_sky_ns`).
A horizontal photon in vacuum (no profile) at the origin passes a telescope at (L, 0, 0), whose
spherical altitude is obs + sqrt(L^2 + (R + obs)^2) - R - obs = L^2 / (2 (R + obs)) to second order;
the sphere row of such a photon is the straight entry, with no curvature in the path (the flat
frame is the Earth-centred frame translated, so straight lines stay straight). Only this vacuum
case is tested for the spherical geometry (`test_sphere.cpp`); the spherical refraction of the
integrator is implemented but has no independent reference yet.

## Efficiency status

The efficiency status of the [atmosphere](../stream/efficiency-status.md) records the geometry:

| Key | Value |
|---|---|
| `atmosphere.geometry` | `corsika-curved` (mode `corsika`) or `spherical` (mode `spherical`); absent for `none` |
| `atmosphere.earth_radius_m` | Earth radius used [m]; absent for `none` |

## What CORSIKA does

The provenance of the `corsika` mode is the CORSIKA 7 source (`corsika7/src/corsika.F`, CURVED
option):

- The Earth radius is `DATA C /6371315.D2`: 6371.315 km in CORSIKA's centimetre units.
- The CURVED option uses a "sliding planar atmosphere": the shower is tracked in local flat frames,
  which are rotated by the Earth angle each time the horizontal displacement exceeds a limit of
  6 to 20 km (depending on altitude). Cherenkov output is converted to the detector frame.
- Bending and time of flight come from precomputed two-dimensional tables (cos θ against height,
  steps of 0.25 deg and 0.25 km; routines `DISTIP` and `TOFIP`) of numerically integrated rays.
- The CORSIKA documentation (Bernloehr, `bernlohr/README`) states that with CURVED the emission
  altitude is the altitude above sea level, and that the transmission no longer scales exactly
  with the secant of the zenith angle.

## Relation to sim_telarray

sim_telarray does not model curvature in the atmosphere, also not for CORSIKA CURVED input. In
`sim_telarray/common/sim_telarray.c` the airmass of every photon bunch is $1/|c_z|$ ("assume a
plane atmosphere", for the bunches with stored emission altitude and for the 3D bunches), and the
transmission is $\exp[-(|\tau(z_{em}) - \tau_{tel}|\,\text{airmass} + \tau_{focus})]$
(`sim_absorb.c`, `atmospheric_transmission2`), the plane-parallel formula of mode `none`. A CORSIKA
CURVED run only makes sim_telarray print that CORSIKA was compiled with the CURVED option.

The curved modes of this stage therefore differ from the atmosphere of sim_telarray at large
zenith angles by the curved-minus-plane-parallel correction of the table in [](#why), by design. To
use the curved extinction in a chain, apply it in this stage and give sim_telarray an all-zero
optical-depth table (`atmosphere` is `applied` upstream, see [](../stream/efficiency-status.md)).

## Choice of radius

The local radius of curvature varies with latitude and azimuth by at most 0.5 %. For WGS84 it is
about 6366.6 km at La Palma and about 6364.2 km (Gaussian radius) at Paranal. This changes the
curvature correction, which is at most a few per cent, by less than $10^{-3}$ relative. The
default is therefore one constant, the CORSIKA radius. The user can override it with
`earth_radius_m`.

## Limitations

- Single scattering is rejected in both curved modes; it is plane-parallel ([](scattering.md)).
- The in-telescope path is not curved ([](propagation.md)).
- The Earth is a sphere, not an ellipsoid.
- The curved extinction integrates the slant depth per photon; its cost compared with the
  plane-parallel airmass has not been measured.
- The curved modes are validated against analytic limits and against integrals computed inside the
  tests. A comparison with CORSIKA CURVED output is not part of these tests; the validation chain
  of the project covers it.
- The `corsika` mode does not apply refraction itself; CORSIKA has already applied it
  ([](refraction.md)).

## Validation

The tests are in `src/propagation/cpp/tests/`:

- `test_curvature.cpp` (slant depth and extinction):
  - The limit $R \to \infty$ equals the plane-parallel result (relative difference below
    $10^{-6}$), for the fixture table and for a synthetic table.
  - The slant depth agrees with a brute-force midpoint integration of the same piecewise integrand.
  - At large zenith angles the curved depth is finite and below the plane-parallel depth.
  - Monotonicity, zero below the observation level, and input validation.
  - The per-row acceptance kernel and the extinction step are independent of row order and
    chunking.
- `test_curved_refraction.cpp` (spherical refraction): the flat limit against `bend()`, an
  independent fine Simpson integration at the real radius, the physical sign and size of the
  integrals, the tables against direct quadrature, tables near the profile top, and input
  validation.
- `test_curvature_modes.cpp` (the stage): the input checks of the three modes, the efficiency status
  keys, the CORSIKA mode extinction, the spherical mode against `none` in the limit
  $R \to \infty$, the emission altitude against $d^2/(2R)$, telescope altitudes between table
  levels, dropped photons, and independence of chunk size and photon order.

The validation against independent references and CORSIKA CURVED showers (suite `curvature`) is
described in the `obdeect-validation` repository (`docs/curvature-validation.md`).
