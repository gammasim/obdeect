# Single scattering

Optional single scattering of light-source photons (plan step 6b, off by default). Each
primary photon is scattered at N_s points along its path; each scattering point gives one
photon per telescope sphere, aimed at the sphere. The primary photon itself is unchanged.

Assumptions: single scattering, no polarisation, straight lines, plane-parallel atmosphere
(no Earth curvature), no refraction. Light travels in vacuum for timing.

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
| tau_mol | molecular (Rayleigh) vertical optical depth from the altitude to the top, >= 0 |
| tau_aer | aerosol (Mie) vertical optical depth from the altitude to the top, >= 0 |
| ssa | aerosol single-scattering albedo, [0, 1] |
| g | aerosol asymmetry parameter, (-1, 1) |

Rows are grouped by wavelength (all altitudes of one wavelength, then the next). Optical
depths must not increase with altitude. Interpolation as in the transmission table: the
wavelength is truncated to integer nm (clamped to the table), each quantity is linear in
log10(altitude) and clamped to the edge altitudes. The table is site independent.
The depth above the top node is taken as the top value, so the top node must be high
enough that tau there is negligible.

The example `cpp/tests/data/scattering_exponential.dat` is an exponential atmosphere:
tau_mol(h) = 0.00865 (lambda/um)^-4.05 exp(-h/8 km), tau_aer(h) = 0.1 (lambda/550 nm)^-1.3
exp(-h/1.5 km), ssa = 0.9, g = 0.65, 340 to 400 nm, nodes 0.1 to 60 km. It is an example, not
a fitted model.

## Phase functions (normalised over 4 pi)

$$
p_R(\theta) = \frac{3}{16\pi}(1 + \cos^2\theta), \qquad
p_{HG}(\theta) = \frac{1 - g^2}{4\pi (1 + g^2 - 2g\cos\theta)^{3/2}}.
$$

$\cos\theta$ is between the propagation direction of the primary and the direction towards
the telescope. The mixed phase function of a segment with optical depths $\tau'_{mol}$,
$\tau'_{aer}$ is

$$
p_{mix} = \frac{\tau'_{mol} p_R + \omega \tau'_{aer} p_{HG}}{\tau'_{mol} + \omega \tau'_{aer}}.
$$

## Next-event estimation

Input: primary photons of frame sky (position, direction, wavelength, weight w at emission).
The path from the emission altitude down to the observation level has slant extinction
optical depth $D_{tot}$ (tau_mol + tau_aer, divided by $|c_z|$). Split it into N_s segments of
equal depth $\Delta D = D_{tot}/N_s$, boundaries $D_j = j\Delta D$. In segment $j$ the
scattering point has density $\propto e^{-D}$ (drawn by inverse CDF). The probability of
scattering in the segment is

$$
P_{scatter} = e^{-D_j}\,(1 - e^{-\Delta D})\,\omega_{mix}, \qquad
\omega_{mix} = \frac{\tau'_{mol} + \omega \tau'_{aer}}{\tau'_{mol} + \tau'_{aer}},
$$

with the segment's vertical depths $\tau'$ (differences of the table at the boundaries) and
$\omega$, g evaluated at the sampled point. For each point and telescope (centre C, radius R,
distance D from the point P) the scattered weight is

$$
w_s = w\, P_{scatter}\, p_{mix}(\theta)\, \Omega_{tel}\, T(P \to C), \qquad
\Omega_{tel} = 2\pi\left(1 - \sqrt{1 - (R/D)^2}\right),
$$

$$
T = \exp\!\left[-\left(\tau_{tab}(h_P) - \tau_{tab}(h_C)\right) \frac{D}{h_P - h_C}\right],
$$

where $\tau_{tab}$ is the optical depth of the transmission table (from the observation level).
The scattered photon is aimed at a uniform point of the sphere's projected disc (radius R,
perpendicular to the line of sight through C), so it lands on the sphere. Points with
$h_P \le h_C$ or $D \le R$ give nothing (photons travel downwards only).

Sum over $j$ of $w_s$ gives the expected scattered weight, unbiased for $N_s \to \infty$.
The approximations are:

- the aim-point disc sampling and the exact sphere solid angle differ at order $(R/D)^2$;
- the segment depths and $\omega$, g are taken from the segment ends and the sampled point;
- the discretisation error falls as $N_s^{-2}$ (test: 3e-4 relative at 32 points, 1.5e-5 at 128).

## Output

Scattered photons, one list per telescope, frame sky, `telescope_id` set, columns
x, y, z (position = scattering point, z above the observation level), cx, cy, cz (towards the
aim point), t, zem (altitude of the scattering point, a.s.l.), w, wavelength. Time is the time
at the scattering point, `t_in + (path from emission to the scattering point) / c`, c in vacuum.
The path from the point to the telescope plane is not included; the atmosphere stage adds it
(see [](propagation.md)).

```text
t_in = 0 ns, emission -> point 15000 m  =>  t = 15000 / c = 50034.6 ns
```

Attention:

- Run this step on primary photons before the extinction step. The factor T is in w, so the
  scattered photons must not go through the atmospheric acceptance (extinction) again. The
  in-telescope part exp(-od_focus) is applied to them (see [](propagation.md)).
- Scattered photons need wavelengths; assign them before this step.
- The primary photon is not changed; its extinction is applied by the extinction step.

## Random numbers

Philox (gstrickd `uniform4`), stage `DLUF`, stream `stream_scattering` (3), keyed on
(seed, event, photon index = index of the primary photon in the event). Deviates 0..N_s-1
place the scattering points; deviates N_s + 2 (j n_tel + t) and +1 are the aim-point offsets
of point j and telescope t. The output does not depend on chunk splits when the chunks are
passed in order.

## Example

```python
import dlufd

table = dlufd.ScatteringTable("scattering_exponential.dat")
site = dlufd.TransmissionTable("atm_trans_2000_1_10_0_0_2000.dat", site_altitude_m=2000.0)
out = dlufd.scatter(
    columns,  # x, y, z, cx, cy, t, w, wavelength (float32), z above the observation level
    table, site,
    telescopes=[(1, 3000.0, 0.0, 0.0, 12.5)],  # id, x, y, z, radius [m]
    seed=1, event=1, first_index=0, n_points=32,
)
photons = out[1]  # dict of float32 arrays
```

## Checks

- Phase functions: normalisation and mean cosine (`test_scattering.cpp`).
- Vertical beam, one telescope at 3 km, ground level, 355 nm, molecular only: Monte Carlo
  per-primary mean 1.736507e-8 +- 1.5e-12 against the numerical integral 1.736294e-8
  (1.4 sigma). The numerical integral without beam attenuation is 1.755344e-8; the ls-beam
  formula gives 1.755214e-8 (7e-5 relative), so the beam attenuation (1.1 %) is the only
  difference, as ls-beam neglects it.
- Example table, aerosol and molecular, 355 nm, 128 points: MC 4.041856e-7 +- 4.7e-12
  against the integral 4.041797e-7 (1.3 sigma).
- Chunk splits give bit-identical output; zero optical depth gives no photons.

## ls-beam correspondence

`sim_telarray/LightEmission/ls-beam.cc`: the emitted photon density per step is
`nphot_per_cm * hscale * steplength_cm` with `hscale = exp(-z/7500 m)` (`ls-beam.cc:314`),
emitted with the Rayleigh angular distribution (`IactLightEmission.cc:1818-1836`, the
tabulated form of (3/16 pi)(1 + cos^2) integrated over azimuth), and each photon reaches the
sphere with the sphere solid angle. ls-beam neglects the attenuation of the laser beam before
the scattering point. The test ports this sum (1 m steps) and matches the analytic integral
without attenuation to 1e-4.

## References

- `cpp/include/dlufd/scattering.hpp`, `python/dlufd/scattering.py`.
- sim_telarray `LightEmission/ls-beam.cc`, `LightEmission/IactLightEmission.cc`.
- Transmission table and optical depths: [extinction](extinction.md).
