# Validation

Manual reference checks. They compare obdeect output with sim_telarray LightEmission programs and
check the stream codec in a container. They need the `ghcr.io/gammasim/simtools-dev` container
(podman) or external programs. They are not run in CI.

Scripts (run from the repository root after building):

- `sources/compare_flasher-dish.sh`: flasher `--type dish` vs `ff-1m`.
  Run: `WORK=/path/out CFG=/path/cfg validation/sources/compare_flasher-dish.sh`
- `sources/compare_flasher-secondary.sh`: flasher `--type secondary` vs `ff-gct`.
  Run: `WORK=/path/out CFG=/path/cfg validation/sources/compare_flasher-secondary.sh`
- `sources/compare_led.sh`: LED vs `pixled`, uses `sources/compare-led-camera.txt`.
  Run: `WORK=/path/out CFG=/path/cfg validation/sources/compare_led.sh`
- `sources/compare_stream.py`: comparison helper run inside the container by the scripts above.
  Run: `python3 validation/sources/compare_stream.py REFERENCE.dat.gz SOURCES.gst`
- `sources/summarize_iact.py`: prints a summary of an IACT EventIO file (inside the container).
  Run: `python3 validation/sources/summarize_iact.py showers.iact.gz`
- `stream/validate_container.sh`: CORSIKA, `stream from-corsika | stream to-simtel`, and sim_telarray
  end-to-end check. Run: `ROOT=/path/checkout validation/stream/validate_container.sh OUTDIR`

The automated end-to-end chain test is `validation/test_chain.py`. Run it with `pytest`.
