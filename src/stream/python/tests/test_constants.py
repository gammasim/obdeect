"""Format constants come from the extension; the Python names are built from them."""

import obdeect.stream as obdeect_stream
from obdeect.stream import stream


def test_format_constants():
    assert obdeect_stream.format_major == 1
    assert obdeect_stream.format_minor == 0
    assert obdeect_stream.max_chunk_rows == 2**16


def test_frame_ids_and_columns_come_from_the_extension():
    assert obdeect_stream.frame_ids == {"sky": 0, "ground": 1, "telescope": 2, "camera": 3}
    assert obdeect_stream.column_names()[:3] == ["x", "y", "z"]
    assert len(obdeect_stream.column_names()) == 11


def test_python_names_match_the_extension():
    assert stream.FRAMES == obdeect_stream.frame_ids
    assert stream.COLUMNS == tuple(obdeect_stream.column_names())
    assert obdeect_stream.FRAMES is stream.FRAMES
    assert obdeect_stream.COLUMNS is stream.COLUMNS
