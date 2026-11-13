"""dlufd and an installed gstrickd agree on the photon-stream format version."""

import re
from pathlib import Path

import pytest

STREAMS_HEADER = Path(__file__).resolve().parents[2] / "cpp" / "include" / "obdeect" / "propagation" / "streams.hpp"


def _expected_format_major():
    match = re.search(
        r"expected_gstrickd_format_major\s*=\s*(\d+)", STREAMS_HEADER.read_text(encoding="utf-8")
    )
    assert match, "expected_gstrickd_format_major not found in streams.hpp"
    return int(match.group(1))


def test_gstrickd_format_major_matches_dlufd():
    obdeect_stream = pytest.importorskip("obdeect.stream")
    if not hasattr(obdeect_stream, "format_major"):
        pytest.skip("installed gstrickd has no format_major")
    assert obdeect_stream.format_major == _expected_format_major()
