"""Validate every ``tests/golden/*.json`` file against the schema in
``docs/dev/conventions.md``.

Structure only -- numeric parity against ``fdatools`` is exercised by
``tests/parity/test_<module>.py``, not here.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

GOLDEN_DIR = Path(__file__).resolve().parent.parent / "golden"
GOLDEN_FILES = sorted(GOLDEN_DIR.glob("*.json"))

REQUIRED_META_KEYS = {"module", "r_version", "fda_version", "generated", "seed", "rtol"}
REQUIRED_CASE_KEYS = {"name", "r_call", "input", "output"}


def _load(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as fh:
        doc: dict[str, Any] = json.load(fh)
        return doc


@pytest.fixture(params=GOLDEN_FILES, ids=[p.stem for p in GOLDEN_FILES])
def golden_path(request: pytest.FixtureRequest) -> Path:
    """One golden file path per parametrization."""
    path: Path = request.param
    return path


def test_golden_files_exist() -> None:
    assert GOLDEN_FILES, f"no golden files found in {GOLDEN_DIR}"


def test_top_level_keys(golden_path: Path) -> None:
    doc = _load(golden_path)
    assert set(doc.keys()) == {"meta", "cases"}


def test_meta_schema(golden_path: Path) -> None:
    doc = _load(golden_path)
    meta = doc["meta"]
    missing = REQUIRED_META_KEYS - meta.keys()
    assert not missing, f"{golden_path.name}: meta missing keys {missing}"
    assert meta["module"] == golden_path.stem
    assert meta["fda_version"] == "6.3.0"
    assert isinstance(meta["r_version"], str)
    assert meta["r_version"]
    assert isinstance(meta["generated"], str)
    assert meta["generated"]
    assert isinstance(meta["seed"], int)
    assert isinstance(meta["rtol"], float)
    assert meta["rtol"] > 0


def test_cases_schema(golden_path: Path) -> None:
    doc = _load(golden_path)
    cases = doc["cases"]
    assert isinstance(cases, list)
    assert len(cases) > 0, f"{golden_path.name}: no cases"

    names = []
    for case in cases:
        missing = REQUIRED_CASE_KEYS - case.keys()
        assert not missing, f"{golden_path.name}: case missing keys {missing}: {case.get('name')}"
        assert isinstance(case["name"], str)
        assert case["name"]
        assert isinstance(case["r_call"], str)
        assert case["r_call"]
        assert isinstance(case["input"], dict)
        assert isinstance(case["output"], dict)
        assert case["output"], f"{golden_path.name}: empty output in case {case['name']}"
        if "rtol" in case:
            assert isinstance(case["rtol"], float), (
                f"{golden_path.name}: case {case['name']} rtol must be a float"
            )
            assert case["rtol"] > 0, f"{golden_path.name}: case {case['name']} rtol must be > 0"
        names.append(case["name"])

    dupes = {n for n in names if names.count(n) > 1}
    assert not dupes, f"{golden_path.name}: duplicate case names {dupes}"
