"""Profile plots (optional matplotlib)."""

import pytest
from obdeect.sources.plot import profile_figure

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

FLASHER = {
    "frame": "camera",
    "photons": 1e5,
    "emitter": {
        "position_m": [0.0, 0.0, 1.5],
        "pointing": [0.0, 0.0, -1.0],
        "angular": {"kind": "lambertian"},
    },
    "spectrum": {"kind": "uniform", "low_nm": 380.0, "high_nm": 420.0},
    "pulse": {"kind": "gauss", "sigma_ns": 1.0},
}


def test_profile_figure_has_angle_wavelength_and_time_panels():
    import matplotlib.pyplot as plt

    figure = profile_figure(FLASHER, seed=1, photons=20000)
    assert len(figure.axes) == 3
    assert figure.axes[0].get_xlabel().startswith("angle")
    assert figure.axes[1].get_xlabel().startswith("wavelength")
    assert figure.axes[2].get_xlabel().startswith("emission time")
    plt.close(figure)
