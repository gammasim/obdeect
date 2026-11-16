"""Release files at the repository root must agree with each other.

All components are released together, so these checks belong to the repository, not to one
component. They only read files, so no obdeect package needs to be installed.
"""

import re
from pathlib import Path

CITATION = Path(__file__).resolve().parents[1] / "CITATION.cff"
CHANGELOG = Path(__file__).resolve().parents[1] / "CHANGELOG.md"


def test_citation_version_matches_latest_changelog_release():
    latest = re.search(r"^## \[v?(\d+\.\d+\.\d+)\]", CHANGELOG.read_text(encoding="utf-8"), re.M)
    assert latest is not None, "CHANGELOG.md has no release heading"
    cited = re.search(r"^version:\s*(\S+)", CITATION.read_text(encoding="utf-8"), re.M)
    assert cited is not None, "CITATION.cff has no version"
    assert cited.group(1) == latest.group(1), (
        f"CITATION.cff version {cited.group(1)} must equal the latest CHANGELOG.md release "
        f"{latest.group(1)}"
    )


def test_citation_date_matches_latest_changelog_release():
    latest = re.search(
        r"^## \[v?\d+\.\d+\.\d+\] - (\d{4}-\d{2}-\d{2})",
        CHANGELOG.read_text(encoding="utf-8"),
        re.M,
    )
    assert latest is not None, "CHANGELOG.md has no dated release heading"
    released = re.search(r"^date-released:\s*(\S+)", CITATION.read_text(encoding="utf-8"), re.M)
    assert released is not None, "CITATION.cff has no date-released"
    assert released.group(1) == latest.group(1), (
        f"CITATION.cff date-released {released.group(1)} must equal the latest CHANGELOG.md "
        f"release date {latest.group(1)}"
    )
