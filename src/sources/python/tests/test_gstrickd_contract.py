"""dlamb.source keeps Python copies of gstrickd constants. They must equal the values of
the gstrickd headers compiled into dlamb (always checked) and of an installed gstrickd
Python package (checked when it can be imported)."""

import pytest
from dlamb import _core, source

CONTRACT = _core.gstrickd_contract()


def test_frames_match_compiled_gstrickd():
    assert source.FRAMES == CONTRACT["frames"]


def test_chunk_rows_match_compiled_gstrickd():
    assert source.CHUNK_ROWS == CONTRACT["max_chunk_rows"]


def test_columns_match_compiled_gstrickd():
    assert source.COLUMNS == tuple(CONTRACT["emitted_columns"])


def test_format_major_matches_dlamb_expectation():
    assert CONTRACT["format_major"] == CONTRACT["dlamb_format_major"]


def test_matches_gstrickd_package():
    gstrickd = pytest.importorskip("gstrickd")
    assert source.FRAMES == {name: gstrickd.frame_ids[name] for name in source.FRAMES}
    assert source.CHUNK_ROWS == gstrickd.max_chunk_rows
    columns = [name for name in gstrickd.column_names() if name in source.COLUMNS]
    assert source.COLUMNS == tuple(columns)
    assert gstrickd.format_major == CONTRACT["format_major"]
