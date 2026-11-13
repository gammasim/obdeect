"""obdeect.sources.source keeps Python copies of obdeect stream constants.
They must equal the values of the obdeect stream headers compiled into
obdeect sources (always checked) and of an installed obdeect stream Python
package (checked when it can be imported)."""

import pytest
from obdeect.sources import _core, source

CONTRACT = _core.stream_contract()


def test_frames_match_compiled_gstrickd():
    assert source.FRAMES == CONTRACT["frames"]


def test_chunk_rows_match_compiled_gstrickd():
    assert source.CHUNK_ROWS == CONTRACT["max_chunk_rows"]


def test_columns_match_compiled_gstrickd():
    assert source.COLUMNS == tuple(CONTRACT["emitted_columns"])


def test_format_major_matches_dlamb_expectation():
    assert CONTRACT["format_major"] == CONTRACT["sources_format_major"]


def test_matches_gstrickd_package():
    obdeect_stream = pytest.importorskip("obdeect.stream")
    assert source.FRAMES == {name: obdeect_stream.frame_ids[name] for name in source.FRAMES}
    assert source.CHUNK_ROWS == obdeect_stream.max_chunk_rows
    columns = [name for name in obdeect_stream.column_names() if name in source.COLUMNS]
    assert source.COLUMNS == tuple(columns)
    assert obdeect_stream.format_major == CONTRACT["format_major"]
