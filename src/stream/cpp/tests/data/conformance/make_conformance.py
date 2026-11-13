"""Writes the conformance streams and expected.json, using only the standard library.

An independent implementation of docs/format.md can read valid.ps and must reproduce the
content in expected.json, and must reject every file listed under "invalid".
Run: python make_conformance.py  (in this directory)
"""

import json
import struct
from pathlib import Path

MAGIC = 0x52545347
TYPES = {
    "stream": 1,
    "run": 2,
    "event": 3,
    "telescope": 4,
    "photons": 5,
    "telescope_end": 6,
    "event_end": 7,
    "run_end": 8,
    "blob": 9,
    "efficiency_status": 10,
    "stream_end": 15,
}
COLUMNS = ["x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "wavelength", "p"]


def record(kind, payload=b""):
    return struct.pack("<IHHQ", MAGIC, TYPES[kind], 0, len(payload)) + payload


def metadata(entries):
    out = struct.pack("<I", len(entries))
    for key, value in entries.items():
        out += struct.pack("<H", len(key)) + key.encode()
        if isinstance(value, int):
            out += struct.pack("<Bq", 1, value)
        elif isinstance(value, float):
            out += struct.pack("<Bd", 2, value)
        elif isinstance(value, str):
            out += struct.pack("<BI", 3, len(value.encode())) + value.encode()
        else:
            out += struct.pack("<BI", 4, len(value)) + struct.pack(f"<{len(value)}d", *value)
    return out


def photons(telescope, frame, first, columns):
    rows = len(next(iter(columns.values()))) if columns else 0
    mask = sum(1 << COLUMNS.index(name) for name in columns)
    payload = struct.pack("<IIIIQ", telescope, frame, mask, rows, first)
    for name in sorted(columns, key=COLUMNS.index):
        payload += struct.pack(f"<{rows}f", *columns[name])
    return record("photons", payload)


STREAM = {"format_major": 1, "format_minor": 0, "producer": "conformance"}
RUN = {"run_number": 7, "energy_min_gev": 0.25, "name": "site", "telescope_x_m": [0.0, 10.5]}
EFFICIENCY_STATUS = {"factors": "qe", "qe.state": "applied", "qe.source": "test"}
EVENT = {"event_number": 3}
TELESCOPE = {"telescope_id": 2, "frame": 1}
CHUNK = {
    "x": [0.0, 1.0, 2.0],
    "y": [0.5, 1.5, 2.5],
    "cx": [0.0, 0.25, 0.5],
    "cy": [0.0, -0.25, -0.5],
    "t": [10.0, 11.0, 12.0],
    "w": [1.0, 1.0, 0.5],
    "zem": [8000.0, 8100.0, 8200.0],
}
END = {"photons": 2.5}
BLOB = bytes([1, 2, 3, 255])

valid = b"".join([
    record("stream", metadata(STREAM)),
    record("run", metadata(RUN)),
    record("efficiency_status", metadata(EFFICIENCY_STATUS)),
    record("event", metadata(EVENT)),
    record("telescope", metadata(TELESCOPE)),
    photons(2, 1, 0, CHUNK),
    photons(2, 1, 3, {}),
    record("blob", struct.pack("<II", 1204, 9) + BLOB),
    record("telescope_end", metadata(END)),
    record("event_end", metadata({})),
    record("run_end", metadata({})),
    record("stream_end"),
])

# A record of unknown type (99) between two known ones must be skipped by readers.
unknown = struct.pack("<IHHQ", MAGIC, 99, 0, 3) + b"abc"
valid = valid.replace(
    record("event_end", metadata({})), unknown + record("event_end", metadata({}))
)

rows_too_many = struct.pack("<IIIIQ", 1, 1, 0b1_0101_1011, 65537, 0)
offset = valid.index(photons(2, 1, 0, CHUNK))
invalid = {
    "bad_magic.ps": b"XSTR" + valid[4:],
    "unknown_major.ps": record("stream", metadata({**STREAM, "format_major": 2}))
    + record("stream_end"),
    "truncated_chunk.ps": valid[: offset + 16 + 24 + 20],
    "too_many_rows.ps": record("stream", metadata(STREAM))
    + record("photons", rows_too_many)
    + record("stream_end"),
}

here = Path(__file__).parent
(here / "valid.ps").write_bytes(valid)
for name, data in invalid.items():
    (here / name).write_bytes(data)

expected = {
    "format": "records in file order; the unknown-type record is skipped; stream_end is not listed",
    "valid": {
        "file": "valid.ps",
        "records": [
            {"type": "stream", "metadata": STREAM},
            {"type": "run", "metadata": RUN},
            {"type": "efficiency_status", "metadata": EFFICIENCY_STATUS},
            {"type": "event", "metadata": EVENT},
            {"type": "telescope", "metadata": TELESCOPE},
            {"type": "photons", "telescope_id": 2, "frame": 1, "first_index": 0, "columns": CHUNK},
            {"type": "photons", "telescope_id": 2, "frame": 1, "first_index": 3, "columns": {}},
            {"type": "blob", "eventio_type": 1204, "eventio_id": 9, "data": list(BLOB)},
            {"type": "telescope_end", "metadata": END},
            {"type": "event_end", "metadata": {}},
            {"type": "run_end", "metadata": {}},
        ],
    },
    "invalid": {
        "bad_magic.ps": "bad magic",
        "unknown_major.ps": "unsupported format major version",
        "truncated_chunk.ps": "truncated",
        "too_many_rows.ps": "2^16 rows",
    },
}
(here / "expected.json").write_text(json.dumps(expected, indent=1) + "\n")
