"""Release metadata stays consistent across pyproject, the package, CI and CITATION."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import fdatools

# tomllib is standard library from Python 3.11; on 3.10 these checks are skipped
# (the 3.11-3.13 CI jobs still run them).
tomllib = pytest.importorskip("tomllib")

REPO_ROOT = Path(__file__).resolve().parents[2]
PYPROJECT = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
SUPPORTED = ["3.10", "3.11", "3.12", "3.13"]


def test_version_matches_everywhere() -> None:
    version = PYPROJECT["project"]["version"]
    assert fdatools.__version__ == version
    citation = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert re.search(rf"^version: {re.escape(version)}$", citation, re.MULTILINE)
    changelog = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"## [{version}]" in changelog


def test_typed_marker_and_package_data_present() -> None:
    package = REPO_ROOT / "src" / "fdatools"
    assert (package / "py.typed").is_file()
    for name in ("growth", "gait", "pinch"):
        assert (package / "_data" / f"{name}.npz").is_file()
        assert (package / "_data" / f"{name}.json").is_file()


def test_python_versions_declared_and_tested() -> None:
    project = PYPROJECT["project"]
    assert project["requires-python"] == ">=3.10"
    classifiers = set(project["classifiers"])
    for version in SUPPORTED:
        assert f"Programming Language :: Python :: {version}" in classifiers
    ci = (REPO_ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    matrix = '["' + '", "'.join(SUPPORTED) + '"]'
    assert ci.count(f"python-version: {matrix}") >= 2  # unit tests and wheel smoke test


def test_license_is_bsd3_spdx_without_legacy_classifier() -> None:
    project = PYPROJECT["project"]
    assert project["license"] == "BSD-3-Clause"
    assert not any(c.startswith("License ::") for c in project["classifiers"])
    license_text = (REPO_ROOT / "LICENSE").read_text(encoding="utf-8")
    assert license_text.startswith("BSD 3-Clause License")
