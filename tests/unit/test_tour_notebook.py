"""Tests for ``tools/build_tour_notebook.py`` and the committed ``notebooks/tour.ipynb``."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import nbformat

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL = REPO_ROOT / "tools" / "build_tour_notebook.py"

_spec = importlib.util.spec_from_file_location("build_tour_notebook", TOOL)
assert _spec is not None
assert _spec.loader is not None
bt = importlib.util.module_from_spec(_spec)
sys.modules["build_tour_notebook"] = bt
_spec.loader.exec_module(bt)


def test_build_notebook_splits_markdown_and_code() -> None:
    text = "# %% [markdown]\n# # Title\n#\n# Text.\n\n# %%\nx = 1\n"
    notebook = bt.build_notebook(text)
    assert [cell.cell_type for cell in notebook.cells] == ["markdown", "code"]
    assert notebook.cells[0].source == "# Title\n\nText."
    assert notebook.cells[1].source == "x = 1"
    assert notebook.metadata["kernelspec"]["name"] == "python3"


def test_committed_notebook_matches_its_source() -> None:
    assert bt.main(["--check"]) == 0


def test_check_reports_a_stale_notebook(tmp_path: Path) -> None:
    source = tmp_path / "tour.py"
    output = tmp_path / "tour.ipynb"
    source.write_text("# %%\nx = 1\n", encoding="utf-8")
    assert bt.main(["--source", str(source), "--output", str(output)]) == 0
    assert bt.main(["--source", str(source), "--output", str(output), "--check"]) == 0
    source.write_text("# %%\nx = 2\n", encoding="utf-8")
    assert bt.main(["--source", str(source), "--output", str(output), "--check"]) == 1


def test_committed_notebook_has_no_outputs() -> None:
    path = str(REPO_ROOT / "notebooks" / "tour.ipynb")
    notebook = nbformat.read(path, as_version=4)  # type: ignore[no-untyped-call]
    assert all(not cell.get("outputs") for cell in notebook.cells if cell.cell_type == "code")
