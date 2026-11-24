"""Fixtures shared by the optical-model tests: the segmented model and the stream module."""

import pytest
from obdeect.ray_tracing import OpticalModel
from telescopes import segmented


@pytest.fixture
def obdeect_stream():
    return pytest.importorskip("obdeect.stream")


@pytest.fixture(scope="module")
def model():
    records, files = segmented()
    return OpticalModel.from_records(
        records, files.__getitem__, configuration_seed=1, propagation_group_index=1.0
    )
