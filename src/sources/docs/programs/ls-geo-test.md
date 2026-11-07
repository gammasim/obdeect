# ls-geo-test

Geometry checks of `sim_telarray/LightEmission/ls-geo-test.cc`: point, line and plane
distances, rotations and intersections, with the original test points.

```bash
dlamb-ls-geo-test
```

## Output

- This program emits no photons, so it writes a text report to stdout (lengths in metres;
  the original uses cm with the same numbers). It writes no photon stream.
- Rotations use the right-hand rule. The original's `Rotate_*` sign convention was not
  checked against it.
- The report and the checks in `cpp/tests/test_ls_geo_test.cpp` are analytic.
