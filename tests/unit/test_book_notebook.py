"""Tests for ``tools/build_book_notebook.py``, the book-figures notebook assembler."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any

import nbformat
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL = REPO_ROOT / "tools" / "build_book_notebook.py"

_spec = importlib.util.spec_from_file_location("build_book_notebook", TOOL)
assert _spec is not None
assert _spec.loader is not None
bb = importlib.util.module_from_spec(_spec)
sys.modules["build_book_notebook"] = bb
_spec.loader.exec_module(bb)

IMPORTS = """\
# %%
import matplotlib.pyplot as plt

plt.rcParams["figure.max_open_warning"] = 0
"""


def _figure(chapter: int, number: int, last: str = "fig") -> str:
    return f"""
# %% [markdown]
# ### Figure {chapter}.{number}
# A description.

# %%
fig, ax = plt.subplots()
{last}
"""


def _chapter(number: int, figures: list[int], title: str = "A title") -> str:
    head = f"# %% [markdown]\n# # Chapter {number}: {title}\n#\n# Figures listed.\n\n"
    return head + IMPORTS + "".join(_figure(number, m) for m in figures)


def _write(directory: Path, name: str, text: str) -> Path:
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# parse_percent
# --------------------------------------------------------------------------- #


def test_parse_percent_splits_markdown_and_code() -> None:
    cells = bb.parse_percent(_chapter(3, [1]))
    assert [c.kind for c in cells] == ["markdown", "code", "markdown", "code"]
    assert cells[0].source == "# Chapter 3: A title\n\nFigures listed."
    assert cells[2].source == "### Figure 3.1\nA description."
    assert cells[3].source == "fig, ax = plt.subplots()\nfig"


def test_parse_percent_keeps_blank_lines_inside_code() -> None:
    cells = bb.parse_percent("# %%\na = 1\n\nb = 2\n\n\n")
    assert cells == [bb.Cell("code", "a = 1\n\nb = 2")]


def test_parse_percent_rejects_text_before_first_marker() -> None:
    with pytest.raises(bb.BookError, match="before the first cell marker"):
        bb.parse_percent("import os\n# %%\nx = 1\n")


def test_parse_percent_allows_blank_lines_before_first_marker() -> None:
    assert bb.parse_percent("\n\n# %%\nx = 1\n") == [bb.Cell("code", "x = 1")]


def test_parse_percent_rejects_unknown_marker() -> None:
    with pytest.raises(bb.BookError, match="unknown cell marker"):
        bb.parse_percent("# %% [raw]\n# text\n")


def test_parse_percent_rejects_unprefixed_markdown() -> None:
    with pytest.raises(bb.BookError, match="must start with '# '"):
        bb.parse_percent("# %% [markdown]\nplain text\n")


def test_parse_percent_rejects_empty_cell() -> None:
    with pytest.raises(bb.BookError, match="empty code cell"):
        bb.parse_percent("# %%\n\n# %%\nx = 1\n")


def test_parse_percent_markdown_blank_lines() -> None:
    cells = bb.parse_percent("# %% [markdown]\n# one\n#\n\n# two\n")
    assert cells == [bb.Cell("markdown", "one\n\n\ntwo")]


# --------------------------------------------------------------------------- #
# validate_chapter
# --------------------------------------------------------------------------- #


def _validate(name: str, text: str) -> Any:
    return bb.validate_chapter(Path(name), bb.parse_percent(text, name))


def test_validate_chapter_collects_figures_and_title() -> None:
    chapter = _validate("ch05_smoothing.py", _chapter(5, [1, 2, 4], title="Smoothing"))
    assert chapter.number == 5
    assert chapter.title == "Smoothing"
    assert chapter.figures == ((5, 1), (5, 2), (5, 4))


def test_validate_chapter_without_figures() -> None:
    assert _validate("ch02_matlab_r.py", _chapter(2, [])).figures == ()


@pytest.mark.parametrize("name", ["chapter1.py", "ch1_intro.py", "ch01-intro.py", "ch01_Intro.py"])
def test_validate_chapter_rejects_bad_file_names(name: str) -> None:
    with pytest.raises(bb.BookError, match="chNN_<slug>"):
        _validate(name, _chapter(1, [1]))


def test_validate_chapter_rejects_wrong_heading() -> None:
    text = _chapter(1, [1]).replace("# # Chapter 1:", "# # Kapitel 1:")
    with pytest.raises(bb.BookError, match="first cell must be markdown"):
        _validate("ch01_intro.py", text)


def test_validate_chapter_rejects_code_first_cell() -> None:
    with pytest.raises(bb.BookError, match="first cell must be markdown"):
        _validate("ch01_intro.py", IMPORTS + IMPORTS)


def test_validate_chapter_rejects_heading_number_mismatch() -> None:
    with pytest.raises(bb.BookError, match="heading says chapter 2, file says 1"):
        _validate("ch01_intro.py", _chapter(2, []))


def test_validate_chapter_needs_two_cells() -> None:
    with pytest.raises(bb.BookError, match="heading cell and an import cell"):
        _validate("ch01_intro.py", "# %% [markdown]\n# # Chapter 1: T\n")


def test_validate_chapter_needs_max_open_warning() -> None:
    text = _chapter(1, [1]).replace('plt.rcParams["figure.max_open_warning"] = 0', "x = 0")
    with pytest.raises(bb.BookError, match="max_open_warning"):
        _validate("ch01_intro.py", text)


def test_validate_chapter_accepts_single_quoted_max_open_warning() -> None:
    text = _chapter(1, [1]).replace('"figure.max_open_warning"', "'figure.max_open_warning'")
    assert _validate("ch01_intro.py", text).figures == ((1, 1),)


def test_validate_chapter_rejects_figure_without_fig_at_end() -> None:
    text = _chapter(1, []) + _figure(1, 1, last="plt.close(fig)")
    with pytest.raises(bb.BookError, match=r"Figure 1\.1 must end with 'fig'"):
        _validate("ch01_intro.py", text)


def test_validate_chapter_rejects_figure_followed_by_markdown() -> None:
    text = _chapter(1, []) + "\n# %% [markdown]\n# ### Figure 1.1\n\n# %% [markdown]\n# more\n"
    with pytest.raises(bb.BookError, match="not followed by a code cell"):
        _validate("ch01_intro.py", text)


def test_validate_chapter_rejects_figure_at_end_of_file() -> None:
    text = _chapter(1, []) + "\n# %% [markdown]\n# ### Figure 1.1\n"
    with pytest.raises(bb.BookError, match="not followed by a code cell"):
        _validate("ch01_intro.py", text)


def test_validate_chapter_accepts_fig_with_trailing_comment() -> None:
    text = _chapter(1, []) + _figure(1, 1, last="fig  # noqa: B018")
    assert _validate("ch01_intro.py", text).figures == ((1, 1),)


def test_validate_chapter_rejects_fig_prefix() -> None:
    text = _chapter(1, []) + _figure(1, 1, last="figure")
    with pytest.raises(bb.BookError, match="must end with 'fig'"):
        _validate("ch01_intro.py", text)


def test_validate_chapter_rejects_figure_of_other_chapter() -> None:
    text = _chapter(1, []) + _figure(2, 1)
    with pytest.raises(bb.BookError, match=r"Figure 2\.1 does not belong to chapter 1"):
        _validate("ch01_intro.py", text)


def test_validate_chapter_rejects_malformed_figure_heading() -> None:
    text = _chapter(1, []) + _figure(1, 1).replace("### Figure 1.1", "### Figure 1.1 (growth)")
    with pytest.raises(bb.BookError, match="malformed figure heading"):
        _validate("ch01_intro.py", text)


def test_validate_chapter_rejects_buried_figure_heading() -> None:
    text = _chapter(1, []) + _figure(1, 1).replace("# ### Figure 1.1", "# Intro\n# ### Figure 1.1")
    with pytest.raises(bb.BookError, match="first line of its cell"):
        _validate("ch01_intro.py", text)


@pytest.mark.parametrize(
    ("code", "problem"),
    [("plt.show()", r"plt\.show"), ("from fdatools._internal import x", "_internal")],
)
def test_validate_chapter_rejects_forbidden_code(code: str, problem: str) -> None:
    text = _chapter(1, [1]).replace("fig, ax = plt.subplots()", code)
    with pytest.raises(bb.BookError, match=problem):
        _validate("ch01_intro.py", text)


# --------------------------------------------------------------------------- #
# validate_book, discover, load_chapters
# --------------------------------------------------------------------------- #


def test_validate_book_accepts_increasing_figures(tmp_path: Path) -> None:
    _write(tmp_path, "ch01_a.py", _chapter(1, [1, 2]))
    _write(tmp_path, "ch03_b.py", _chapter(3, [1]))
    chapters = bb.load_chapters(tmp_path)
    assert [c.number for c in chapters] == [1, 3]


def test_validate_book_rejects_duplicate_figure(tmp_path: Path) -> None:
    _write(tmp_path, "ch01_a.py", _chapter(1, [1, 1]))
    with pytest.raises(bb.BookError, match=r"Figure 1\.1 is duplicated or out of order"):
        bb.load_chapters(tmp_path)


def test_validate_book_rejects_decreasing_figures(tmp_path: Path) -> None:
    _write(tmp_path, "ch01_a.py", _chapter(1, [2, 1]))
    with pytest.raises(bb.BookError, match=r"Figure 1\.1 is duplicated or out of order"):
        bb.load_chapters(tmp_path)


def test_validate_book_rejects_duplicate_chapter(tmp_path: Path) -> None:
    _write(tmp_path, "ch01_a.py", _chapter(1, [1]))
    _write(tmp_path, "ch01_b.py", _chapter(1, [2]))
    with pytest.raises(bb.BookError, match="chapter 1 is duplicated"):
        bb.load_chapters(tmp_path)


def test_discover_orders_by_chapter_number(tmp_path: Path) -> None:
    for name in ["ch10_x.py", "ch02_y.py", "ch01_z.py"]:
        _write(tmp_path, name, "")
    assert [p.name for p in bb.discover(tmp_path)] == ["ch01_z.py", "ch02_y.py", "ch10_x.py"]


def test_discover_sorts_badly_named_files_last(tmp_path: Path) -> None:
    _write(tmp_path, "chapter.py", "")
    _write(tmp_path, "ch01_a.py", "")
    assert [p.name for p in bb.discover(tmp_path)] == ["ch01_a.py", "chapter.py"]


def test_discover_needs_chapter_files(tmp_path: Path) -> None:
    with pytest.raises(bb.BookError, match="no chapter files"):
        bb.discover(tmp_path)


# --------------------------------------------------------------------------- #
# build_notebook, figure_table, main
# --------------------------------------------------------------------------- #


def test_build_notebook_structure(tmp_path: Path) -> None:
    _write(tmp_path, "ch01_a.py", _chapter(1, [1, 2]))
    notebook = bb.build_notebook(bb.load_chapters(tmp_path))
    assert notebook.nbformat == 4
    assert notebook.metadata["kernelspec"]["name"] == "python3"
    kinds = [cell.cell_type for cell in notebook.cells]
    assert kinds == ["markdown", "markdown", "code", "markdown", "code", "markdown", "code"]
    assert "FDATOOLS_DATA_DIR" in notebook.cells[0].source
    assert "data-v1" in notebook.cells[0].source
    code = [cell for cell in notebook.cells if cell.cell_type == "code"]
    assert all(cell.outputs == [] and cell.execution_count is None for cell in code)
    nbformat.validate(notebook)


def test_figure_table_counts() -> None:
    chapters = [
        bb.validate_chapter(Path("ch01_a.py"), bb.parse_percent(_chapter(1, [1, 2]))),
        bb.validate_chapter(Path("ch03_b.py"), bb.parse_percent(_chapter(3, [1]))),
    ]
    table = bb.figure_table(chapters).splitlines()
    assert table[1].split() == ["1", "2", "ch01_a.py"]
    assert table[2].split() == ["3", "1", "ch03_b.py"]
    assert table[-1].split() == ["total", "3"]


def test_main_writes_notebook(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    source = tmp_path / "book"
    source.mkdir()
    _write(source, "ch01_a.py", _chapter(1, [1]))
    output = tmp_path / "out" / "book.ipynb"
    assert bb.main(["--source", str(source), "--output", str(output), "--expect", "1"]) == 0
    notebook = nbformat.read(str(output), as_version=4)  # type: ignore[no-untyped-call]
    assert notebook.cells[-1].source.endswith("fig")
    assert "wrote" in capsys.readouterr().out


def test_main_check_writes_nothing(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write(tmp_path, "ch01_a.py", _chapter(1, [1]))
    output = tmp_path / "book.ipynb"
    assert bb.main(["--source", str(tmp_path), "--output", str(output), "--check"]) == 0
    assert not output.exists()
    assert "check passed" in capsys.readouterr().out


def test_main_expect_mismatch(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write(tmp_path, "ch01_a.py", _chapter(1, [1, 2]))
    assert bb.main(["--source", str(tmp_path), "--check", "--expect", "76"]) == 1
    assert "expected 76 figures, found 2" in capsys.readouterr().err


def test_main_reports_contract_error(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write(tmp_path, "ch01_a.py", _chapter(1, [1, 1]))
    assert bb.main(["--source", str(tmp_path), "--check"]) == 1
    assert "duplicated or out of order" in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# the real chapter files
# --------------------------------------------------------------------------- #


def test_repository_chapters_pass_check(capsys: pytest.CaptureFixture[str]) -> None:
    assert bb.main(["--check"]) == 0
    out = capsys.readouterr().out
    assert "check passed" in out
    chapters = bb.load_chapters(bb.SOURCE_DIR)
    assert chapters[0].number == 1
    assert all(chapter.figures or chapter.number == 2 for chapter in chapters)
