"""The package and its native extension import."""

import obdeect.propagation as obdeect_propagation
from obdeect.propagation import _core


def test_version():
    assert isinstance(obdeect_propagation.__version__, str)
    assert len(_core.core_version()) == 3
