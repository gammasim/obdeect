"""The native constructors enforce the parameter counts of each kind and validate on creation."""

import pytest
from obdeect.sources import _core


@pytest.mark.parametrize(
    ("make", "kind"),
    [
        (_core.make_angular, "isotropic"),
        (_core.make_angular, "lambertian"),
        (_core.make_angular, "rayleigh"),
        (_core.make_spectrum, "table"),
        (_core.make_pulse, "table"),
    ],
)
def test_zero_parameter_kinds_reject_parameters(make, kind):
    with pytest.raises(ValueError, match="wrong number of parameters"):
        make(kind, [1.0], [0.0, 1.0], [1.0, 1.0])


def test_source_is_validated_on_creation():
    emitter = _core.Emitter([0, 0, 1], [0, 0, -1], _core.make_angular("isotropic", [], [], []))
    spectrum = _core.make_spectrum("fixed", [400.0], [], [])
    pulse = _core.make_pulse("fixed", [0.0], [], [])
    with pytest.raises(ValueError):
        _core.Source(0, emitter, spectrum, pulse, -1.0)
