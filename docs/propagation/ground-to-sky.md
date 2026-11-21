# Ground frame to sky frame

CORSIKA photons from `obdeect stream from-corsika` are in frame G (ground). Single scattering
(`obdeect propagation run --scattering`) needs photons in the sky frame. `obdeect propagation to-sky`
converts G to the sky frame. It is the inverse of the sky-to-G step of the
[atmosphere stage](propagation.md).

## Conversion

Per photon, with telescope position $(X_t, Y_t, Z_t)$ [m], observation level $h_\mathrm{obs}$
[m a.s.l.], $h_t = h_\mathrm{obs} + Z_t$, $c_z = \sqrt{1 - c_x^2 - c_y^2}$ and
$\zeta = \operatorname{atan2}(\sqrt{c_x^2 + c_y^2},\, c_z)$:

$$
s = \frac{z_\mathrm{em} - h_t}{c_z}, \qquad
x_\mathrm{sky} = X_t + x - s\,c_x, \qquad
y_\mathrm{sky} = Y_t + y - s\,c_y, \qquad
z_\mathrm{sky} = z_\mathrm{em} - h_\mathrm{obs}
$$

$$
t_\mathrm{sky} = t - \left(\frac{s}{c} + \Delta t\right), \qquad
\Delta t = \mathrm{straight\_path\_air\_delay\_ns}(\text{profile}, h_t, z_\mathrm{em}, \zeta)
$$

$\Delta t$ is the straight-path air delay of [timing](refraction.md) and is zero without a profile.
The stage adds the same terms for a run with `refraction_applied = 1`
(`src/propagation/cpp/include/obdeect/propagation/sky.hpp`, lines 297 to 327).
Columns `zem`, `w`, `wavelength`, `cx`, `cy` are kept. `z` is added. Telescope positions come from
the run keys `telescope_x_m`, `telescope_y_m`, `telescope_z_m`, the keys that the stage uses. Telescope
record keys `x_m`, `y_m`, `z_m` must agree with them if present. Telescope records get frame sky.
All other records pass through unchanged.

Rows that cannot be converted (`z_em <= h_t`, `c_z` not downwards, non-finite values) are removed.
The count is printed as `dropped` in the stderr summary. Remaining rows keep their photon index, so
the random draws of later stages do not change.

## Example

```bash
obdeect stream from-corsika run.iact \
  | obdeect propagation to-sky --profile atmprof.dat \
  | obdeect propagation run --scattering ...
```

`run --no-extinction` with the same `--profile` option reverses the conversion.

In Python (the same conversion, returns the counts):

```python
import obdeect.propagation as propagation

summary = propagation.to_sky("ground.pst", "sky.pst", profile="atmprof.dat")
# {'photons_in': ..., 'photons_out': ..., 'dropped': 0}
back = propagation.propagate(
    "sky.pst",
    "ground_again.pst",
    extinction=False,
    site_altitude_m=2000.0,
    seed=1,
    profile="atmprof.dat",
)
```

## Limits

- Round trip: exact to float32 rounding when the input has `refraction_applied = 1` (CORSIKA with
  refraction) or when no profile is used. With `--profile`, a run with `refraction_applied` missing
  or 0 is rejected, because bending is not inverted.
- The emission point is the straight back-projection along the arrival direction. Curvature of the
  refracted path is neglected: the error is of order metres at large zenith angles (see
  [refraction](refraction.md)).
- Rejected input: `curvature_applied = 1` (CORSIKA CURVED photons have no flat sky geometry),
  frames other than G, a missing `zem` column, a missing observation level, and a telescope record
  whose position differs from the run.
- Times are float32: the rounding of the converted time is about half the float32 spacing at the
  larger of $t$ and $t_\mathrm{sky}$ (about 0.002 ns near 3e4 ns).
