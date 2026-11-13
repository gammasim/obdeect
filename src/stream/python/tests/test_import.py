"""The package and its native extension import."""

import gstrickd
from gstrickd import _core


def test_version():
    assert isinstance(gstrickd.__version__, str)
    assert len(_core.core_version()) == 3
