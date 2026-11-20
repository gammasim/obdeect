# Single scattering

Light that a source emits can scatter in the air towards the telescope. obdeect computes this
optional single scattering (off by default) for light-source photons. Each primary photon is
scattered at N_s points along its path. Each scattering point gives one photon per telescope
sphere, aimed at the sphere. The primary photon is unchanged.

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

Input: primary photons of frame sky (position, direction, wavelength, weight w at emission). The path
from the emission altitude down to the observation level has slant extinction optical depth
$D_{tot}$ (tau_mol + tau_aer, divided by $|c_z|$). Split it into N_s segments of equal depth
$\Delta D = D_{tot}/N_s$, boundaries $D_j = j\Delta D$. In segment $j$ the scattering point has density
$\propto e^{-D}$ (drawn by inverse CDF). The probability of scattering in the segment is

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
T = \exp\!\left[-\left(\tau_{tab}(h_P) - \tau_{tab}(h_C)\right) \frac{D}{h_P - h_C}\right],
$$

where $\tau_{tab}$ is the optical depth of the transmission table (from the observation level). The
scattered photon is aimed at a uniform point of the sphere's projected disc (radius R, perpendicular
to the line of sight through C), so it lands on the sphere. Points with $h_P \le h_C$ or $D \le R$ give
nothing (photons travel downwards only).

The sum over $j$ of $w_s$ is the expected scattered weight, unbiased for $N_s \to \infty$. Approximations:

- the aim-point disc sampling and the exact sphere solid angle differ at order $(R/D)^2$;
- the segment depths and $\omega$, g are taken from the segment ends and the sampled point;
- the tabulated aerosol phase is linear in $\cos\theta$ between the angle nodes (test table: 4e-4
  relative to Henyey-Greenstein at worst, 1 degree nodes), and linear in log altitude and wavelength
  between the table rows; its normalisation is the trapezoid rule in $\cos\theta$;
- the air delay (with an air profile) uses the straight path at the primary's zenith angle, not a
  refracted path. It is evaluated from $|c_z|$ directly, so it stays finite for any nonzero $|c_z|$;
- the discretisation error falls as $N_s^{-2}$ (test: 3e-4 relative at 32 points, 1.5e-5 at 128).

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
