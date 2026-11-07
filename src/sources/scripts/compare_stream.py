"""Compare a sim_telarray LightEmission IACT file with a dlamb gstrickd stream.

Usage (inside the simtools-dev container, which has eventio and numpy):
    python3 compare_stream.py ORIGINAL.dat.gz DLAMB.gst

Per event: photon count (original: sum of bunch sizes; dlamb: sum of weights).
Over all photons: means of cx, cy, cz, time [ns], wavelength [nm] with standard errors.
Original cz = -sqrt(1 - cx^2 - cy^2) (IACT stores w = -cz with cz < 0 downwards).
Positions are not compared: the original writes them at the observation plane.
Times: the original adds the flight time to the observation plane; dlamb writes emission
times, so the t difference equals the path length / c (ff-1m: 16 m, 53.4 ns).
"""

import struct
import sys

import numpy as np
from eventio.iact import IACTFile

MAGIC = 0x52545347
RECORD_EVENT = 3
RECORD_PHOTONS = 5
COLUMNS = ["x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "wavelength", "p"]


def read_original(path):
    counts, cols = [], {"cx": [], "cy": [], "t": [], "wavelength": []}
    with IACTFile(path) as f:
        for event in f:
            bunches = event.photon_bunches
            n = 0.0
            for tel in bunches.values():
                n += float(np.sum(tel["photons"]))
                cols["cx"].append(tel["cx"])
                cols["cy"].append(tel["cy"])
                # Original times are arrival times at the observation plane, relative to
                # time_offset (the event's first interaction); add it back.
                cols["t"].append(tel["time"] + event.time_offset)
                cols["wavelength"].append(tel["wavelength"])
            counts.append(n)
    out = {k: np.concatenate(v) if v else np.zeros(0) for k, v in cols.items()}
    out["cz"] = -np.sqrt(np.clip(1.0 - out["cx"] ** 2 - out["cy"] ** 2, 0.0, None))
    return counts, out


def read_dlamb(path):
    counts, cols = [], {k: [] for k in ["cx", "cy", "cz", "t", "wavelength"]}
    weight = None
    with open(path, "rb") as f:
        data = f.read()
    pos = 0
    while pos + 16 <= len(data):
        magic, rtype, _flags, length = struct.unpack_from("<IHHQ", data, pos)
        if magic != MAGIC:
            raise ValueError("not a gstrickd stream")
        payload = data[pos + 16 : pos + 16 + length]
        pos += 16 + length
        if rtype == RECORD_EVENT:
            if weight is not None:
                counts.append(weight)
            weight = 0.0
        elif rtype == RECORD_PHOTONS:
            _tel, _frame, mask, rows = struct.unpack_from("<IIII", payload, 0)
            body = np.frombuffer(payload, dtype="<f4", offset=24)
            arrays, k = {}, 0
            for bit, name in enumerate(COLUMNS):
                if mask >> bit & 1:
                    arrays[name] = body[k * rows : (k + 1) * rows].astype(np.float64)
                    k += 1
            weight += float(np.sum(arrays["w"]))
            for name in cols:
                if name in arrays:
                    cols[name].append(arrays[name])
    if weight is not None:
        counts.append(weight)
    return counts, {k: np.concatenate(v) if v else np.zeros(0) for k, v in cols.items()}


def moments(a):
    n = len(a)
    if n == 0:
        return float("nan"), float("nan")
    return float(np.mean(a)), float(np.std(a) / np.sqrt(n)) if n > 1 else float("nan")


def main(original, dlamb):
    oc, oo = read_original(original)
    dc, dd = read_dlamb(dlamb)
    print(f"events: original {len(oc)}, dlamb {len(dc)}")
    if oc and dc:
        print(
            f"photons per event: original mean {np.mean(oc):.2f} (sum {np.sum(oc):.0f}), "
            f"dlamb mean {np.mean(dc):.2f} (sum {np.sum(dc):.0f})"
        )
    print(f"{'quantity':12s} {'original':>24s} {'dlamb':>24s} {'diff/sigma':>11s}")
    for name in ["cx", "cy", "cz", "t", "wavelength"]:
        mo, so = moments(oo[name])
        md, sd = moments(dd[name])
        sig = np.hypot(so, sd)
        ratio = (md - mo) / sig if sig > 0 else float("nan")
        print(f"{name:12s} {mo:12.6g} +- {so:8.2g} {md:12.6g} +- {sd:8.2g} {ratio:11.2f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
