"""The package and its native extension import."""

import obdeect.stream as obdeect_stream
from obdeect.stream import _core


def test_version():
    assert isinstance(obdeect_stream.__version__, str)
    assert len(_core.core_version()) == 3
