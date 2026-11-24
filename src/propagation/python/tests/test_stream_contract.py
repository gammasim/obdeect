"""obdeect propagation and an installed obdeect stream agree on the photon-stream format version."""

import re
from pathlib import Path

STREAMS_HEADER = (
    Path(__file__).resolve().parents[2]
    / "cpp"
    / "include"
    / "obdeect"
    / "propagation"
    / "streams.hpp"
)


def _expected_format_major():
    match = re.search(
        r"expected_gstrickd_format_major\s*=\s*(\d+)", STREAMS_HEADER.read_text(encoding="utf-8")
    )
    assert match, "expected_gstrickd_format_major not found in streams.hpp"
    return int(match.group(1))


def test_stream_format_major_matches_dlufd():
    import obdeect.stream as obdeect_stream

    # No skip: an installed stream without format_major is a broken installation.
    assert hasattr(obdeect_stream, "format_major"), "obdeect stream has no format_major"
    assert obdeect_stream.format_major == _expected_format_major()
