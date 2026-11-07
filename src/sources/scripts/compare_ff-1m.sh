#!/bin/bash
# Cross-check dlamb-ff-1m (default options, 20 events) against sim_telarray LightEmission ff-1m
# in the simtools-dev container. Needs the built dlamb tool (build/debug/dlamb-ff-1m) and podman.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${WORK:-/private/tmp/claude-502/xcheck/ff-1m}"
CFG="${CFG:-/Users/maierg/DESYCloud/WIP/cta/202611/full-simulation-chain/sim_telarray/cfg/common}"
mkdir -p "$WORK"
"$ROOT/build/debug/dlamb-ff-1m" --events 20 -o "$WORK/dlamb.gst"
podman run --rm -v "$CFG:/cfg:ro" -v "$WORK:/out" -v "$ROOT/scripts:/scripts:ro" \
  ghcr.io/gammasim/simtools-dev:latest sh -c \
  "cd /out && /workdir/simulation_software/sim_telarray/LightEmission/ff-1m -I/cfg --events 20 -o /out/orig.dat.gz > /out/orig.log 2>&1 && python3 /scripts/compare_stream.py /out/orig.dat.gz /out/dlamb.gst"
