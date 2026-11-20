"""Fixed optical models compiled from the synthetic records in the former compiler tests."""

from pathlib import Path

from obdeect.ray_tracing import OpticalModel

_FIXTURES = Path(__file__).parent / "data" / "optical-models"


def load_model(name: str) -> OpticalModel:
    """Load a checked-in compiled optical model by fixture name."""
    return OpticalModel.load(_FIXTURES / f"{name}.json")
