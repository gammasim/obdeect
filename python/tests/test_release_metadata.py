"""Release metadata must agree with the version of the tag being built."""

import re
from pathlib import Path

import obdeect
import pytest

CITATION = Path(__file__).resolve().parents[2] / "CITATION.cff"


def test_citation_version_matches_release():
    if not CITATION.exists():
        pytest.skip("CITATION.cff is not part of this source tree")
    version = obdeect.__version__
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        pytest.skip(f"{version} is not a release version")
    cited = re.search(r"^version:\s*(\S+)", CITATION.read_text(encoding="utf-8"), re.M)
    assert cited is not None
    assert cited.group(1) == version, "CITATION.cff version must equal the release version"
