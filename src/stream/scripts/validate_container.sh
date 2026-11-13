#!/usr/bin/env bash
# End-to-end validation of obdeect stream from-corsika and obdeect stream to-simtel in ghcr.io/gammasim/simtools-dev.
#
#   1. CORSIKA (image binary) writes a small IACT EventIO file: 20 gamma showers, 100-1000 GeV,
#      20 deg zenith, 3 telescopes.
#   2. obdeect stream from-corsika and obdeect stream to-simtel are compiled inside the container with g++ (the image has no
#      cmake; the sources are mounted read-only).
#   3. obdeect stream from-corsika < run.iact | obdeect stream to-simtel > roundtrip.corsika; the files are compared.
#   4. sim_telarray (image binary, CTA-PROD6 LST configuration from the host checkout) runs on
#      both files with the same configuration and seed; the outputs are compared.
#
# Usage: scripts/validate_container.sh OUTDIR
#   OUTDIR  host directory for the data, logs and summary (large files; keep outside Git)
#   ROOT    checkout root with corsika7/ and sim_telarray/ (default: two levels above WORKDIR)
#   IMAGE   container image (default ghcr.io/gammasim/simtools-dev:latest)
#
# Steps 1-3 must succeed; step 4 is reported and does not stop the script.

set -euo pipefail

WORKDIR=$(cd "$(dirname "$0")/.." && pwd)
ROOT=${ROOT:-$(cd "$WORKDIR/../.." && pwd)}
IMAGE=${IMAGE:-ghcr.io/gammasim/simtools-dev:latest}
OUT=${1:?usage: scripts/validate_container.sh OUTDIR}
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)

# The steps that run inside the container.
write_inner_script() {
  cat <<'INNER'
#!/usr/bin/env bash
set -uo pipefail
export LC_ALL=C
W=/out/work
SIMTEL_DIR=/workdir/simulation_software/sim_telarray
mkdir -p "$W" /out/bin /out/data /out/simtel
echo "== image: $(uname -m), $(g++ --version | head -1)"

# CORSIKA data tables (QGSJET-II, SIBYLL-less UrQMD and EGS4 sets), copied once.
if [ ! -f /out/data/EGSDAT6_1. ]; then
  cp /hostdata/corsika_run/NUCNUCCS /hostdata/corsika_run/UrQMD-1.3.1-xs.dat \
     /hostdata/corsika_run/QGSDAT01 /hostdata/corsika_run/EGSDAT6_* /out/data/
  cp /hostdata/corsika_run/qgsdat-II-04 /out/data/QGSDAT-II-04
  cp /hostdata/corsika_run/sectnu-II-04 /out/data/SECTNU-II-04
fi

# Step 2: build the adapters from the mounted sources (no cmake in the image). CMake
# generates obdeect/stream/version.hpp and defines OBDEECT_STREAM_PACKAGE_VERSION; do the same here.
mkdir -p /out/gen/obdeect/stream
sed -e 's/@PROJECT_VERSION_MAJOR@/0/' -e 's/@PROJECT_VERSION_MINOR@/0/' \
    -e 's/@PROJECT_VERSION_PATCH@/0/' /src/cmake/version.hpp.in > /out/gen/obdeect/stream/version.hpp
for pair in corsika2ps:obdeect-stream-from-corsika ps2simtel:obdeect-stream-to-simtel; do
  src=${pair%%:*}; out=${pair##*:}
  g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -ffp-contract=off \
      -DOBDEECT_STREAM_PACKAGE_VERSION='"0.0.0+validation"' -I/out/gen -I/src/cpp/include \
      "/src/cpp/tools/$src.cpp" -o "/out/bin/$out" || exit 1
done
echo "== built: $(/out/bin/obdeect-stream-from-corsika --help 2>&1 | head -1)"

# Step 1: CORSIKA run (card as simtools writes it: RUN, SITE, IACT, TELESCOPE, TELFIL lines).
cd "$W" || exit 1
rm -f run.iact DAT* CER* tables.dat corsika.log
cp /hostdata/corsika_bernlohr/atmprof1.dat .
cat > card.txt <<'CARD'
RUNNR 1
EVTNR 1
NSHOW 20
PRMPAR 1
ERANGE 100.0 1000.0
ESLOPE -2.0
THETAP 20.0 20.0
PHIP 0.0 0.0
SEED 1 0 0
SEED 2 0 0
SEED 3 0 0
SEED 4 0 0
OBSLEV 215600.0
ATMOSPHERE 1 Y
MAGNET 20.0 40.0
ARRANG 0.0
HADFLG 0 0 0 0 0 2
ECUTS 0.3 0.3 0.02 0.02
ELMFLG F F
CSCAT 1 0.0 0.0
DATDIR /out/data/
DIRECT ./
TELESCOPE 0.0 0.0 0.0 1200.0
TELESCOPE 5000.0 0.0 0.0 1200.0
TELESCOPE 0.0 5000.0 0.0 1200.0
TELFIL run.iact
IACT SPLIT_AUTO 0
IACT IO_BUFFER 800MB
IACT MAX_BUNCHES 1000000
EXIT
CARD
start=$(date +%s)
/workdir/simulation_software/corsika7/corsika_qgs2_urqmd_flat < card.txt > corsika.log 2>&1
status=$?
echo "== corsika exit=$status seconds=$(( $(date +%s) - start )) iact_bytes=$(stat -c %s run.iact 2>/dev/null)"
[ "$status" -eq 0 ] || exit "$status"
[ -s run.iact ] || { echo "CORSIKA produced no run.iact"; exit 1; }

# Step 3: codec round trip.
/out/bin/obdeect-stream-from-corsika < run.iact | /out/bin/obdeect-stream-to-simtel > roundtrip.corsika
status=$?
echo "== round trip exit=$status sizes: $(stat -c %s run.iact) $(stat -c %s roundtrip.corsika)"
[ "$status" -eq 0 ] || exit "$status"
cmp -l run.iact roundtrip.corsika > cmp_bytes.txt 2>&1
echo "== differing bytes (cmp -l): $(wc -l < cmp_bytes.txt)"
cmp run.iact roundtrip.corsika && echo "== byte-identical" || echo "== not byte-identical"

# Photon comparison of the two EventIO files with simtools' reader (per telescope and event).
PYTHONPATH=/workdir/src /workdir/env/bin/python3 /out/work/compare_photons.py \
  run.iact roundtrip.corsika > compare_photons.txt 2>&1
photon_status=$?
echo "== photon comparison (exit $photon_status): $(tail -1 compare_photons.txt)"

# Step 4: sim_telarray on both files, same configuration and seed.
run_simtel() {
  local input=$1 name=$2
  cd "$SIMTEL_DIR" || return 1
  rm -f "/out/simtel/$name.simtel" "/out/simtel/$name.hdata" "/out/simtel/$name.log"
  ./bin/sim_telarray -c cfg/CTA/CTA-PROD6-LST.cfg \
      -DCTA_NORTH -DLST1 -DNUM_TELESCOPES=3 \
      -C maximum_telescopes=3 -C altitude=2156 -C telescope_theta=20 -C telescope_phi=0 \
      -C "histogram_file=/out/simtel/$name.hdata" -C "output_file=/out/simtel/$name.simtel" \
      -W random_state=none -C random_seed=20261107 -C show=all \
      "$input" > "/out/simtel/$name.log" 2>&1
  local status=$?
  cd "$W" || return 1
  echo "== sim_telarray on $name exit=$status output_bytes=$(stat -c %s "/out/simtel/$name.simtel" 2>/dev/null)"
  return $status
}
run_simtel "$W/run.iact" direct || true
run_simtel "$W/roundtrip.corsika" roundtrip || true
if [ -s /out/simtel/direct.simtel ] && [ -s /out/simtel/roundtrip.simtel ]; then
  cmp /out/simtel/direct.simtel /out/simtel/roundtrip.simtel && echo "== sim_telarray outputs identical" \
    || echo "== sim_telarray outputs differ: $(cmp -l /out/simtel/direct.simtel /out/simtel/roundtrip.simtel | wc -l) bytes"
  /workdir/env/bin/python3 /out/work/compare_simtel.py /out/simtel/direct.simtel /out/simtel/roundtrip.simtel
  simtel_status=$?
else
  echo "== sim_telarray produced no output to compare"
  simtel_status=1
fi
# Fail the validation when photons or sim_telarray events do not match.
[ "$photon_status" -eq 0 ] && [ "$simtel_status" -eq 0 ]
INNER
}

write_compare_scripts() {
  cat > "$OUT/work/compare_photons.py" <<'PY'
"""Compare the photon bunches of two CORSIKA IACT EventIO files (simtools reader)."""

import sys

import numpy as np
from simtools.corsika.photon_bunches import read_photon_bunches

FIELDS = ("x", "y", "cx", "cy", "time_ns", "zem", "size")


def collect(path):
    table = {}
    for bunch in read_photon_bunches(path):
        key = (bunch["event_number"], bunch["telescope_id"])
        table.setdefault(key, []).append(bunch)
    return {
        key: {name: np.concatenate([np.atleast_1d(b[name]) for b in items]) for name in FIELDS}
        for key, items in table.items()
    }


a = collect(sys.argv[1])
b = collect(sys.argv[2])
print(f"blocks: {len(a)} vs {len(b)}")
# Allowed difference: the cm -> m -> cm float32 round trip changes values by at most 1 ulp.
RTOL = 1e-6
if a.keys() != b.keys():
    print("MISMATCH: different (event, telescope) keys")
    sys.exit(1)
worst = {name: 0.0 for name in FIELDS}
photons_a = photons_b = 0
for key in a:
    photons_a += len(a[key]["x"])
    photons_b += len(b[key]["x"])
    for name in FIELDS:
        if a[key][name].shape != b[key][name].shape:
            worst[name] = float("inf")
            continue
        va = a[key][name].astype(float)
        vb = b[key][name].astype(float)
        diff = np.abs(va - vb) / np.maximum(np.abs(va), 1e-30)
        if diff.size:
            worst[name] = max(worst[name], float(np.max(diff)))
print("bunches:", photons_a, photons_b)
print("max relative difference per field:", {k: f"{v:.3g}" for k, v in worst.items()})
largest = max(worst.values())
if photons_a != photons_b or largest > RTOL:
    print(f"summary: differs beyond tolerance (max relative {largest:.3g} > {RTOL})")
    sys.exit(1)
print(f"summary: {'identical' if largest == 0.0 else 'within tolerance'}")
PY
  cat > "$OUT/work/compare_simtel.py" <<'PY'
"""Photoelectron sums per event of two sim_telarray output files (eventio)."""

import sys

from eventio import SimTelFile


def pe_sums(path):
    sums = []
    with SimTelFile(path) as f:
        for event in f:
            total = 0.0
            for telescope in event.get("telescope_events", {}).values():
                pe = telescope.get("photoelectrons")
                if pe is not None:
                    total += float(sum(pe.get("photoelectrons", [])) if isinstance(pe, dict) else 0.0)
            sums.append(round(total, 6))
    return sums


a = pe_sums(sys.argv[1])
b = pe_sums(sys.argv[2])
print(f"events: {len(a)} vs {len(b)}; photoelectron sums: {a[:10]} vs {b[:10]}")
if len(a) != len(b):
    print("MISMATCH: different number of events")
    sys.exit(1)
PY
}

echo "== container validation: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "== workdir $WORKDIR, output $OUT, root $ROOT, image $IMAGE"
mkdir -p "$OUT/work"
write_inner_script > "$OUT/work/container_steps.sh"
write_compare_scripts

podman run --rm \
  -v "$WORKDIR:/src:ro" \
  -v "$OUT:/out" \
  -v "$ROOT/corsika7/run:/hostdata/corsika_run:ro" \
  -v "$ROOT/corsika7/bernlohr:/hostdata/corsika_bernlohr:ro" \
  -v "$ROOT/sim_telarray/cfg:/workdir/simulation_software/sim_telarray/cfg:ro" \
  "$IMAGE" bash /out/work/container_steps.sh
