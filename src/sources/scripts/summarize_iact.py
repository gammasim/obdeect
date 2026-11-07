"""Summarise photon bunches of an IACT (EventIO) file written by a LightEmission program.

Prints per telescope: bunches, sum of bunch sizes (photons), mean wavelength and the mean
and r.m.s. of the bunch impact positions at the observation level [m]. Run inside a
container with python eventio, e.g. `python3 summarize_iact.py figure-eight.iact.gz`.
"""

import sys

import numpy as np
from eventio import iact


def main(path):
    totals = {}
    for event in iact.IACTFile(path):
        for tel, bunches in event.photon_bunches.items():
            totals.setdefault(tel, []).append(bunches)
    for tel, arrays in sorted(totals.items()):
        b = np.concatenate(arrays) if arrays else np.zeros(0)
        if b.size == 0:
            print(f"telescope {tel}: 0 bunches")
            continue
        x, y = b["x"] / 100.0, b["y"] / 100.0
        print(
            f"telescope {tel}: bunches {b.size}, photons {float(np.sum(b['photons'])):.1f}, "
            f"wavelength mean {float(np.mean(b['wavelength'])):.3f} nm, "
            f"impact mean ({float(np.mean(x)):.3f}, {float(np.mean(y)):.3f}) m, "
            f"rms ({float(np.std(x)):.3f}, {float(np.std(y)):.3f}) m"
        )


if __name__ == "__main__":
    main(sys.argv[1])
