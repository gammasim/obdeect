"""The package and its native extension import."""

import obdeect.sources as obdeect_sources
from obdeect.sources import _core


def test_version():
    assert isinstance(obdeect_sources.__version__, str)
    assert len(_core.core_version()) == 3
