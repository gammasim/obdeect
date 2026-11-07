#!/bin/bash
# Cross-check dlamb-flasher --type secondary against sim_telarray LightEmission/ff-gct.cc in the simtools-dev
# container. Photons per flasher reduced to 1e5 (the default 2.5e6 gives ~360 MB of output).
# The original's photon bunches carry cx, cy only (eventio), so cz is not comparable here.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
: "${WORK:?set WORK to a scratch directory for the cross-check outputs}"
: "${CFG:?set CFG to the sim_telarray cfg/common directory}"
mkdir -p "$WORK"
"$ROOT/build/debug/dlamb-flasher" --type secondary --photons 1e5 --events 1 -o "$WORK/dlamb.gst"
podman run --rm -v "$CFG:/cfg:ro" -v "$WORK:/out" -v "$ROOT/scripts:/scripts:ro" \
  ghcr.io/gammasim/simtools-dev:latest sh -c \
  "cd /out && /workdir/simulation_software/sim_telarray/LightEmission/ff-gct -I/cfg --photons 1e5 --events 1 -o /out/orig.dat.gz > /out/orig.log 2>&1 && python3 /scripts/compare_stream.py /out/orig.dat.gz /out/dlamb.gst"
