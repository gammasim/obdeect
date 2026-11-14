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
| `spherical` | 0 (required) | The stage applies curvature: emission altitude from the sphere, per-telescope local vertical and telescope altitude, curved extinction. The output run record gets `curvature_applied` = 1. | sky → G |

The run flag `curvature_applied` is read from the EventIO input: the EventIO importer sets it to 1
when EVTH word 79 equals 2 (CORSIKA CURVED). The stage requires the flag in `corsika` mode and
rejects it in the other two modes, so curvature is neither applied twice nor forgotten.

## Slant depth

Photons travel on straight lines. For a telescope at altitude $h_t$ and a photon with zenith
angle $\theta_t$ at the telescope, the invariant of the straight line in a spherical atmosphere is

$$
p = (R + h_t)\sin\theta_t, \qquad
\cos\theta(h) = \sqrt{1 - \frac{p^2}{(R+h)^2}},
$$

and the slant optical depth from the telescope up to the emission altitude $z_{em}$ is

$$
\tau_{\text{slant}} = \int_{h_t}^{z_{em}} \frac{1}{\cos\theta(h)}\,\frac{d\tau}{dh}\,dh.
$$

The table is piecewise: $\tau$ is linear in $\log_{10} h$ between the table heights (see
[](extinction.md)), so on segment $i$

$$
\frac{d\tau}{dh} = \frac{a_i}{h \ln 10}, \qquad
a_i = \frac{\tau_{i+1} - \tau_i}{\log_{10} h_{i+1} - \log_{10} h_i}.
$$

- Each segment, clipped to $[h_t, z_{em}]$, is integrated exactly in the table form above, with
  8-point Gauss-Legendre quadrature per segment.
- For $\cos\theta_t < 0.2$ the integrand behaves like $(h - h_t)^{-1/2}$ at the telescope. The first
  segment then uses $h = h_t + u^2$, which removes the singularity.
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
\cos\theta_t = -\,\mathbf{d}\cdot\mathbf{n}.
$$

Photons moving away from the ground ($\cos\theta_t$ clearly negative) have $T = 0$.

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

## Choice of radius

The local radius of curvature varies with latitude and azimuth by at most 0.5 %. For WGS84 it is
about 6366.6 km at La Palma and about 6364.2 km (Gaussian radius) at Paranal. This changes the
curvature correction, which is at most a few per cent, by less than $10^{-3}$ relative. The
default is therefore one constant, the CORSIKA radius. The user can override it with
`earth_radius_m`.

## Limitations

- Curved refraction (bending) and curved air delay are not implemented. The `spherical` mode
  rejects an atmospheric profile; the `corsika` mode relies on CORSIKA having applied them
  ([](refraction.md)).
- Single scattering is rejected in both curved modes; it is plane-parallel ([](scattering.md)).
- The in-telescope path is not curved ([](propagation.md)).
- The Earth is a sphere, not an ellipsoid.
- Planned: two-dimensional ray tables for curved bending in a follow-up.

## Validation

The tests are in `src/propagation/cpp/tests/test_curvature.cpp`:

- The limit $R \to \infty$ equals the plane-parallel result (relative difference below $10^{-6}$),
  for the fixture table and for a synthetic table.
- The slant depth agrees with a brute-force midpoint integration of the same piecewise integrand.
- At large zenith angles the curved depth is finite and below the plane-parallel depth.
- Monotonicity, zero below the observation level, and input validation.
- The per-row acceptance kernel and the extinction step are independent of row order and chunking.
