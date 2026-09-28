"""Assemble ``notebooks/tour.ipynb`` from ``notebooks/tour.py``.

Usage::

    .venv/bin/python tools/build_tour_notebook.py [--output PATH] [--check]

``notebooks/tour.py`` is a plain Python file in the jupytext "percent" format
(see ``tools/build_book_notebook.py``, whose parser this tool reuses), so ruff
lints it like any other source file.  The notebook is written without outputs;
run it in Jupyter or with ``pytest --nbmake notebooks/tour.ipynb``.

``--check`` writes nothing and exits with ``1`` when the committed notebook is
not the one ``tour.py`` produces.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from collections.abc import Sequence
from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO_ROOT / "notebooks" / "tour.py"
OUTPUT = REPO_ROOT / "notebooks" / "tour.ipynb"


def _parse_percent() -> object:
    """Return ``parse_percent`` from the book notebook tool next to this file."""
    path = Path(__file__).with_name("build_book_notebook.py")
    spec = importlib.util.spec_from_file_location("build_book_notebook", path)
    if spec is None or spec.loader is None:  # pragma: no cover - the file ships with the repo
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("build_book_notebook", module)
    spec.loader.exec_module(module)
    return module.parse_percent


def build_notebook(text: str, name: str = "<string>") -> nbformat.NotebookNode:
    """Build the tour notebook from percent-format source.

    Parameters
    ----------
    text : str
        Contents of ``notebooks/tour.py``.
    name : str, optional
        File name used in error messages.

    Returns
    -------
    nbformat.NotebookNode
        A version 4 notebook with no outputs and a ``python3`` kernelspec.
    """
    parse = _parse_percent()
    cells = [
        new_markdown_cell(cell.source)  # type: ignore[no-untyped-call]
        if cell.kind == "markdown"
        else new_code_cell(cell.source)  # type: ignore[no-untyped-call]
        for cell in parse(text, name)  # type: ignore[operator]
    ]
    notebook: nbformat.NotebookNode = new_notebook(  # type: ignore[no-untyped-call]
        cells=cells,
        metadata={
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python"},
        },
    )
    nbformat.validate(notebook)
    return notebook


def _sources(notebook: nbformat.NotebookNode) -> list[tuple[str, str]]:
    """Return ``(cell_type, source)`` for every cell, ignoring ids and outputs."""
    return [(cell.cell_type, cell.source) for cell in notebook.cells]


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line tool.

    Parameters
    ----------
    argv : sequence of str, optional
        Arguments; defaults to ``sys.argv[1:]``.

    Returns
    -------
    int
        Exit status: ``0`` on success, ``1`` when ``--check`` finds a stale notebook.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, default=SOURCE, help="percent-format source")
    parser.add_argument("--output", type=Path, default=OUTPUT, help="notebook to write")
    parser.add_argument("--check", action="store_true", help="compare only, write nothing")
    args = parser.parse_args(argv)

    notebook = build_notebook(args.source.read_text(encoding="utf-8"), str(args.source))
    if args.check:
        current = nbformat.read(str(args.output), as_version=4)  # type: ignore[no-untyped-call]
        if _sources(current) != _sources(notebook):
            print(f"error: {args.output} is stale; rebuild it from {args.source}", file=sys.stderr)
            return 1
        print("check passed")
        return 0
    nbformat.write(notebook, str(args.output))  # type: ignore[no-untyped-call]
    print(f"wrote {args.output} ({len(notebook.cells)} cells)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
