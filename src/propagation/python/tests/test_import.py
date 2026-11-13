"""The package and its native extension import."""

import dlufd
from dlufd import _core


def test_version():
    assert isinstance(dlufd.__version__, str)
    assert len(_core.core_version()) == 3
