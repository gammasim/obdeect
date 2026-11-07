#!/bin/bash
# Cross-check dlamb-pixled against sim_telarray LightEmission pixled (test camera, pixels
# 0 1 2 with two pixels per event, 10 photons per pixel). Both write pixel positions.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${WORK:-/private/tmp/claude-502/xcheck/pixled}"
CFG="${CFG:-/Users/maierg/DESYCloud/WIP/cta/202611/full-simulation-chain/sim_telarray/cfg/common}"
mkdir -p "$WORK"
cp "$ROOT/scripts/compare_pixled_camera.txt" "$WORK/camera.txt"
"$ROOT/build/debug/dlamb-pixled" -c "$WORK/camera.txt" --use "0 1 2" --require 2 --photons 10 -o "$WORK/dlamb.gst"
podman run --rm -v "$CFG:/cfg:ro" -v "$WORK:/out" -v "$ROOT/scripts:/scripts:ro" \
  ghcr.io/gammasim/simtools-dev:latest sh -c \
  "cd /out && /workdir/simulation_software/sim_telarray/LightEmission/pixled -I/cfg -c /out/camera.txt --use '0 1 2' --require 2 --photons 10 -o /out/orig.dat.gz > /out/orig.log 2>&1 && python3 /scripts/compare_stream.py /out/orig.dat.gz /out/dlamb.gst"
