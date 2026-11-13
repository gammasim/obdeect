"""Format constants come from the extension; the Python names are built from them."""

import gstrickd
from gstrickd import stream


def test_format_constants():
    assert gstrickd.format_major == 1
    assert gstrickd.format_minor == 0
    assert gstrickd.max_chunk_rows == 2**16


def test_frame_ids_and_columns_come_from_the_extension():
    assert gstrickd.frame_ids == {"sky": 0, "ground": 1, "telescope": 2, "camera": 3}
    assert gstrickd.column_names()[:3] == ["x", "y", "z"]
    assert len(gstrickd.column_names()) == 11


def test_python_names_match_the_extension():
    assert stream.FRAMES == gstrickd.frame_ids
    assert stream.COLUMNS == tuple(gstrickd.column_names())
    assert gstrickd.FRAMES is stream.FRAMES
    assert gstrickd.COLUMNS is stream.COLUMNS
