# ls-test

The octo figure-eight with `4e7` photons per position. Port of
`sim_telarray/LightEmission/ls-test.cc`.

```bash
dlamb-ls-test --seed 1 -o ls-test.ps
```

## Notes

- `ls-test.cc` passes `nb_per_pos = 4e7` as the photon number of `Emit()` (`ls-test.cc:53`,
  `:87`), so it is octo with half the photons. The name suggests bunches; the code uses
  photons. The port follows the code.
- The original does not build as written: `main()` (`ls-test.cc:45`) takes no arguments but
  uses `argc` and `argv` (line 47). The Makefile does not list it. The port is written from
  the intended logic (see the octo page for geometry and output).
- Options: `--photons` (default `4e7`), `--seed`, `-o`. Output is the same as octo.
