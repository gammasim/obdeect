# Single scattering

Light that a source emits can scatter in the air towards the telescope. obdeect computes this
optional single scattering (off by default) for light-source photons. Each primary photon is
scattered at N_s points along its path. Each scattering point gives one photon per telescope
sphere, aimed at the sphere. The primary photon is unchanged.

Primaries may move downwards or upwards. A downward primary is scattered on its path down to the
observation level, an upward primary on its path up to the top of the scattering table (no scattering
above it). A horizontal primary (cz = 0) is not scattered. Scattered photons move in any direction:
a scattering point below the telescope centre sends light upwards to the sphere.

Assumptions: single scattering, no polarisation (the Rayleigh depolarisation enters only the phase
function), straight lines, plane-parallel atmosphere (no Earth curvature), no refraction. Timing is
in vacuum, or with the air delay of the straight path when the stage passes an air profile (see
Output). Scattering is rejected in the curved
modes of the atmosphere stage ([](curvature.md)); a spherical treatment of the scattering geometry
is not implemented.

## Scattering table

Plain text, `#` comments, one row per (wavelength, altitude):

```text
# wavelength_nm altitude_km tau_mol tau_aer ssa g
340 0.1 6.74679963e-01 1.74821941e-01 0.9 0.65
340 0.5 6.41775433e-01 1.33901079e-01 0.9 0.65
```

| column | meaning |
|---|---|
| wavelength_nm | integer nm, consecutive, at least 1 |
| altitude_km | a.s.l., positive, strictly ascending; the same list for every wavelength |
| tau_mol | molecular (Rayleigh) vertical optical depth from the altitude to the top, ≥ 0 |
| tau_aer | aerosol (Mie) vertical optical depth from the altitude to the top, ≥ 0 |
| ssa | aerosol single-scattering albedo, [0, 1] |
| g | aerosol asymmetry parameter, (-1, 1) |

Rows are grouped by wavelength (all altitudes of one wavelength, then the next). Optical depths must
not increase with altitude. Interpolation: each quantity is linear in log10(altitude) and linear in
wavelength between the integer rows (a wavelength between two rows is interpolated, not truncated);
both are clamped to the table edges. The table is site independent. The depth above the top node is taken as the top value,
so the top node must be high enough that tau there is negligible.

### Optional lines

Keyword lines may be added to the table. Existing tables without them stay valid.

```text
angles_deg 0 90 180
depolarization 0.03
phase 340 0.1 0.0795774715 0.0795774715 0.0795774715
phase 340 0.5 0.0795774715 0.0795774715 0.0795774715
...
```

- `angles_deg`: at most once; 2 to 721 values in degrees, strictly ascending, first 0, last 180.
- `phase wavelength_nm altitude_km p_1 ... p_K`: one row per (wavelength, altitude) pair of the
  main rows; every pair is required once any phase row is present; K is the number of angles; each
  value is finite and > 0 in sr^-1, in the order of `angles_deg`. The normalisation
  $2\pi \int p\, d(\cos\theta)$ (trapezoid rule in $\cos\theta$ over the nodes) must be within 5 % of 1;
  the row is then rescaled to exactly 1 in that rule. Phase rows need `angles_deg`, and `angles_deg`
  needs phase rows. The example is the isotropic function $1/(4\pi)$.
- `depolarization rho`: at most once; 0 <= rho < 0.5; without it rho = 0.
- Any other keyword is an error. The main rows keep g (required); with phase rows g is kept for
  reference and is not used for the phase function.

`cpp/tests/data/scattering_tabulated.dat` has the phase rows of a Henyey-Greenstein function with
g = 0.65 (1 degree nodes) at 354 and 355 nm; its header gives the provenance.

`cpp/tests/data/scattering_exponential.dat` is an exponential example table, not a fitted model:
tau_mol(h) = 0.00865 (lambda/um)^-4.05 exp(-h/8 km), tau_aer(h) = 0.1 (lambda/550 nm)^-1.3
exp(-h/1.5 km), ssa = 0.9, g = 0.65, 340 to 400 nm, nodes 0.1 to 60 km.

## Phase functions (normalised over 4π)

$$
p_R(\theta) = \frac{3}{8\pi(2 + \gamma)}\left[(1 + \gamma) + (1 - \gamma)\cos^2\theta\right], \qquad
\gamma = \frac{\rho}{2 - \rho}, \qquad
p_{HG}(\theta) = \frac{1 - g^2}{4\pi (1 + g^2 - 2g\cos\theta)^{3/2}}.
$$

For $\rho = 0$ the Rayleigh function is $\frac{3}{16\pi}(1 + \cos^2\theta)$. The aerosol function $p_{aer}$
is the tabulated function (optional phase rows, interpolated in $\cos\theta$, log altitude and
wavelength) if the table has phase rows, otherwise $p_{HG}$. Both are normalised over $4\pi$.

$\cos\theta$ is between the primary's propagation direction and the direction towards the telescope.
The mixed phase function of a segment with optical depths $\tau'_{mol}$, $\tau'_{aer}$ is

$$
p_{mix} = \frac{\tau'_{mol} p_R + \omega \tau'_{aer} p_{aer}}{\tau'_{mol} + \omega \tau'_{aer}}.
$$

## Next-event estimation

Input: primary photons of frame sky (position, direction, wavelength, weight w at emission). The
direction is the signed `cz` column if present (negative downward, positive upward, zero horizontal
and not scattered); without `cz` a photon moves downward. The path from the emission altitude to the
end of the scattering (the observation level for downward photons, the top of the table for upward
photons) has slant extinction optical depth $D_{tot}$ (tau_mol + tau_aer, divided by $|c_z|$; for an
upward photon $\tau(h_0) - \tau(h_{top})$). Split it into N_s segments of equal depth
$\Delta D = D_{tot}/N_s$, boundaries $D_j = j\Delta D$. In segment $j$ the scattering point has density
$\propto e^{-D}$ (drawn by inverse CDF, along the path of the photon's direction). The probability of scattering in the segment is

$$
P_{scatter} = e^{-D_j}\,(1 - e^{-\Delta D})\,\omega_{mix}, \qquad
\omega_{mix} = \frac{\tau'_{mol} + \omega \tau'_{aer}}{\tau'_{mol} + \tau'_{aer}},
$$

with the segment's vertical depths $\tau'$ (differences of the table at the boundaries) and $\omega$, g
(and the tabulated phase, if any) evaluated at the sampled point. For each point and telescope (centre C, radius R, distance D from the
point P) the scattered weight is

$$
w_s = w\, P_{scatter}\, p_{mix}(\theta)\, \Omega_{tel}\, T(P \to C), \qquad
\Omega_{tel} = 2\pi\left(1 - \sqrt{1 - (R/D)^2}\right),
$$

$$
T = \exp\!\left[-\left|\tau_{tab}(h_P) - \tau_{tab}(h_C)\right| \frac{D}{|h_P - h_C|}\right],
$$

where $\tau_{tab}$ is the optical depth of the transmission table (from the observation level, growing
with altitude). The scattered photon is aimed at a uniform point of the sphere's projected disc (radius
R, perpendicular to the line of sight through C), so it lands on the sphere. The scattered photon
travels down if $h_P > h_C$ and up if $h_P < h_C$; the same solid angle, aim point and phase function
apply in both cases. A point with $h_P = h_C$ would emit horizontal photons, which are not part of the
sphere's phase function: it gives nothing. Points with $D \le R$ give nothing.

The angle $\theta$ is between the primary's direction and the direction to the aim point, so
$\cos\theta = c_x n_x + c_y n_y + c_z n_z$ with the signed $c_z$.

The sum over $j$ of $w_s$ is the expected scattered weight, unbiased for $N_s \to \infty$. Approximations:

- the aim-point disc sampling and the exact sphere solid angle differ at order $(R/D)^2$;
- the segment depths and $\omega$, g are taken from the segment ends and the sampled point;
- the tabulated aerosol phase is linear in $\cos\theta$ between the angle nodes (test table: 4e-4
  relative to Henyey-Greenstein at worst, 1 degree nodes), and linear in log altitude and wavelength
  between the table rows; its normalisation is the trapezoid rule in $\cos\theta$;
- the air delay (with an air profile) uses the straight path at the primary's zenith angle, not a
  refracted path. It is evaluated from $|c_z|$ directly, so it stays finite for any nonzero $|c_z|$;
- the discretisation error falls as $N_s^{-2}$ (test: 3e-4 relative at 32 points, 1.5e-5 at 128).

### Upward primaries

An upward primary (cz > 0) is scattered on its path from the emission altitude up to the highest
altitude of the scattering table. As for downward primaries, a primary emitted at or below the
observation level (z <= 0) is not scattered: a laser on the ground must start above it (for example
at z = 1 m), and the path below that start point is not simulated. The boundaries are at equal slant depth from the emission, with the
depth $\tau(h_0) - \tau(h)$ (the table's optical depth counts from the altitude to the top). The
scattering point $P$ is at $s = (h_P - h_0)/|c_z|$ along the path. The time at $P$ is
$t_{in} + s / c$, plus the air delay of the straight path from $h_0$ to $h_P$ with an air profile.
The scattered light from a point below the telescope centre travels upwards to the sphere; this is
the only way light from a source below the telescope reaches it.

### Scattering below the telescope (changed behaviour)

Before this change, a scattering point at or below the telescope centre gave nothing (the scattered
photon had to travel downwards). Now the scattered photon travels upwards to the sphere when
$h_P < h_C$, with the transmission and airmass of the path between the two heights. Light is added,
nothing is removed: for downward primaries the scattered weight of a run with telescopes above the
observation level grows by the upward part of the path. Telescopes on the ground (telescope_z_m = 0
for the centre at the observation level) are not changed, since no point lies below their centre.

The photons that point to a part of the sphere below the telescope centre travel downwards at a
point below the centre; in the plane representation they never reach the telescope plane, and the
propagation stage counts them as `moving_away`.

## Sphere runs

With `ground_representation = sphere` the scattering uses the same kernel with two additions
(`ScatteringConfig::sphere_geometry`; plane runs are unchanged):

- Scattered photons whose path to the telescope centre is near horizontal,
  $|h_S - h_C| / d < 0.1$, take the transmission of the chord integral
  ([](extinction.md#chord-integral)) from the scattering point to the centre, $T = \exp(-\int k\,\mathrm{d}s)$.
  The plane form $|\tau(h_S) - \tau(h_C)|\, d / |h_S - h_C|$ is singular at $h_S = h_C$; in sphere
  runs a point at the telescope's height is kept (its photons are horizontal), where the plane
  kernel skips it. Other scattered photons are unchanged.
- Primaries with $|c_z| < 0.1$ (including $c_z = 0$, horizontal) are scattered along their straight
  path. The primary is not bent: the scattering point is at $x_0 + c_x s$, $y_0 + c_y s$, altitude
  $h_0 + c_z s$, and the time $t_0 + s/c$ plus the straight air delay with an air profile. The
  cumulative depth from the emission is

$$
D(s) = \frac{|\tau(h_0 + c_z s) - \tau(h_0)|}{|c_z|},
$$

  the exact depth of a straight path with altitude linear in $s$ (with $\tau$ the scattering table's
  optical depth above an altitude). For $|c_z| < 10^{-9}$ the altitude is constant to the path's
  precision and $D(s) = k(h_0)\,s$, with $k$ from a central difference of the table over 1 m; the
  inverse is then exact. Otherwise $D$ is tabulated at 257 path lengths, uniform in $s$ up to the end
  of the path, and inverted by linear interpolation ($D$ is monotone). The path ends where the
  altitude leaves [observation level, top of the table], where $D$ reaches $-\ln 10^{-9} = 20.72$
  (the rest of the path is not sampled: its scattering probability is below $10^{-9}$), or at
  $10^9$ m. Equal-depth segments are as for the other primaries: the boundaries are the inverse of
  $D$ at $jD_{\text{end}}/N_s$, and the point of segment $j$ is drawn by inverting $e^{-D}$ inside it.
  The molecular and aerosol mix of a segment is taken at its middle (central differences of the
  table), since a horizontal path has zero vertical depth per segment.

Steep primaries ($|c_z| \ge 0.1$) use the downward and upward branches unchanged. A plane run keeps
its behaviour: horizontal primaries (cz = 0) are not scattered and $h_S = h_C$ points are skipped.

Test (`test_scattering`): a horizontal beam at the telescope's height, against the integral
$W = \int_0^{s_{\text{end}}} \alpha e^{-\alpha s}\,P(\cos\theta)\,\Omega(d)\,e^{-\alpha_T d}\,\mathrm{d}s$
with P the Rayleigh phase and $\alpha_T d$ the chord to the centre. The per-primary totals of 100 000
primaries agree with W to 1 % (within the Monte Carlo error, about 1 %).

## Output

Scattered photons, one list per telescope, frame sky, `telescope_id` set, columns x, y, z (position =
scattering point, z above the observation level), cx, cy, cz (towards the aim point), t, zem (altitude
of the scattering point, a.s.l.), w, wavelength. Time is the time at the scattering point,
`t_in + (path from emission to the scattering point) / c`, c in vacuum. With an air profile the time
adds the air delay of the straight path from the scattering point to the emission height at the
primary's zenith angle (`straight_path_air_delay_ns`, [](refraction.md)); the test geometry (emission at
50 km, zenith 12.8 deg, points from 2 to 50 km) gives 0 to 6.4 ns. The path from the point to the
telescope plane is not included; the atmosphere stage adds it ([](propagation.md)).

```text
t_in = 0 ns, emission -> point 15000 m  =>  t = 15000 / c = 50034.6 ns
```

- Run this step on primary photons before the extinction step. The factor T is in w, so the scattered
  photons must not go through the atmospheric acceptance again. Their in-telescope part exp(-od_focus)
  is applied instead ([](propagation.md)).
- Scattered photons need wavelengths; assign them before this step.
- The primary photon is not changed; its extinction is applied by the extinction step.

## Python example

```python
import obdeect.propagation

table = obdeect.propagation.ScatteringTable("scattering_exponential.dat")
site = obdeect.propagation.TransmissionTable(
    "atm_trans_2000_1_10_0_0_2000.dat", site_altitude_m=2000.0
)
out = obdeect.propagation.scatter(
    columns,  # x, y, z, cx, cy, t, w, wavelength (float32), z above the observation level
    table,
    site,
    telescopes=[(1, 3000.0, 0.0, 0.0, 12.5)],  # id, x, y, z, radius [m]
    seed=1,
    event=1,
    first_index=0,
    n_points=32,
)
photons = out[1]  # dict of float32 arrays
```

Random numbers: [](../stream/random-numbers.md). The transmission table and the optical depths are described in [](extinction.md).

Scattered photons keep the signed `cz` of their direction (the aim point direction), so near-horizontal
and upward directions are preserved. The propagation stage (sky to G, [](ground-to-sky.md)) takes the
magnitude of `cz` from this column. An upward scattered photon from below the telescope plane reaches
the plane and is written with `cz > 0`; one from above the plane moves away from it and is not
written.
