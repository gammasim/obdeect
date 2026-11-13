#!/bin/bash
# Cross-check obdeect sources flasher --type dish (default options, 20 events) against the reference ff-1m program
# in the simtools-dev container. Needs the built obdeect sources tool (build/debug/src/sources/obdeect-sources-flasher --type dish) and podman.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
: "${WORK:?set WORK to a scratch directory for the cross-check outputs}"
: "${CFG:?set CFG to the configuration directory of the reference programs}"
mkdir -p "$WORK"
"$ROOT/build/debug/src/sources/obdeect-sources-flasher" --type dish --events 20 -o "$WORK/sources.gst"
podman run --rm -v "$CFG:/cfg:ro" -v "$WORK:/out" -v "$ROOT/validation/sources:/scripts:ro" \
  ghcr.io/gammasim/simtools-dev:latest sh -c \
  "cd /out && /workdir/simulation_software/sim_telarray/LightEmission/ff-1m -I/cfg --events 20 -o /out/orig.dat.gz > /out/orig.log 2>&1 && python3 /scripts/compare_stream.py /out/orig.dat.gz /out/sources.gst"
