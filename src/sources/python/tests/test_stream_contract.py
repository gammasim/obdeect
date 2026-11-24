"""The stream constants of obdeect sources come from the compiled stream contract.

The contract is read from the headers compiled into obdeect sources (always checked); the frame
ids and columns must also equal those of an installed obdeect stream package (checked when it can
be imported)."""

import pytest
from obdeect.sources import _core, source
from obdeect.sources.presets import _common

CONTRACT = _core.stream_contract()

# Values of the stream headers (obdeect/stream/stages.hpp): never change them.
STAGE_PROGRAMS_HEADER = 0x444C0002


def test_frames_match_the_compiled_header():
    assert source.FRAMES == CONTRACT["frames"]


def test_chunk_rows_match_the_compiled_header():
    assert source.CHUNK_ROWS == CONTRACT["max_chunk_rows"]


def test_columns_match_the_compiled_header():
    assert source.COLUMNS == tuple(CONTRACT["emitted_columns"])


def test_format_major_of_sources_is_the_stream_format_major():
    assert CONTRACT["format_major"] == CONTRACT["sources_format_major"]


def test_preset_stage_is_the_stream_programs_stage():
    assert CONTRACT["stages"]["programs"] == STAGE_PROGRAMS_HEADER
    assert _common.STAGE_PRESETS == STAGE_PROGRAMS_HEADER
    assert CONTRACT["stages"]["sources"] == 0x444C0001


def test_matches_installed_stream_package():
    obdeect_stream = pytest.importorskip("obdeect.stream")
    assert source.FRAMES == {name: obdeect_stream.frame_ids[name] for name in source.FRAMES}
    assert source.CHUNK_ROWS == obdeect_stream.max_chunk_rows
    columns = [name for name in obdeect_stream.column_names() if name in source.COLUMNS]
    assert source.COLUMNS == tuple(columns)
    assert obdeect_stream.format_major == CONTRACT["format_major"]
