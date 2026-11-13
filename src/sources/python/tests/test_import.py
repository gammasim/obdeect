"""The package and its native extension import."""

import dlamb
from dlamb import _core


def test_version():
    assert isinstance(dlamb.__version__, str)
    assert len(_core.core_version()) == 3
