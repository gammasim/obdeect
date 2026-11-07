"""Sphinx configuration of the dlamb documentation."""

import tomllib
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]

project = "dlamb"
copyright = "2026, dlamb developers"  # noqa: A001
with open(ROOT / "CITATION.cff", encoding="utf-8") as file:
    author = ", ".join(
        f"{a['given-names']} {a['family-names']}" for a in yaml.safe_load(file)["authors"]
    )
with open(ROOT / "pyproject.toml", "rb") as file:
    distribution = tomllib.load(file)["project"]["name"]
try:
    release = version(distribution)
except PackageNotFoundError:
    release = "unknown"
version = release

extensions = ["myst_parser", "sphinx.ext.autodoc", "sphinx.ext.mathjax", "sphinx.ext.githubpages"]
source_suffix = {".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "changes"]
myst_heading_anchors = 3
myst_enable_extensions = ["colon_fence", "dollarmath"]

html_theme = "pydata_sphinx_theme"
html_title = project
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/gammasim/dlamb",
            "icon": "fa-brands fa-github",
        },
    ],
}
