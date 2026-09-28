# Every target runs tools from the project venv. A bare `uv run` would re-sync
# .venv without the extras and remove pytest, torch and rpy2, so it is not used.
# Run `make sync` once to create or refresh the venv with every extra.
# Override the interpreter with `make PY=/path/to/python <target>`.
PY ?= .venv/bin/python

# Datasets not shipped in the package (GATE 5 runs the book notebook and the
# dataset tests); default to the local release staging directory.
FDATOOLS_DATA_DIR ?= $(CURDIR)/data_release

.PHONY: sync lint type test parity golden bench docs build book gate5

sync:
	uv sync --all-extras

lint:
	$(PY) -m ruff check .
	$(PY) -m ruff format --check .

type:
	$(PY) -m mypy src tests tools benchmarks

test:
	$(PY) -m pytest tests/unit -m "not parity and not slow and not gpu" --cov=fdatools --cov-report=term-missing

parity:
	$(PY) -m pytest tests/unit/test_golden_schema.py tests/parity -m parity

golden:
	$(PY) tools/make_golden.py basis
	$(PY) tools/make_golden.py core

bench:
	$(PY) -m pytest benchmarks --benchmark-only

docs:
	$(PY) -m mkdocs build --strict

build:
	$(PY) -m build
	$(PY) -m twine check dist/*

# Rebuild notebooks/book_figures.ipynb from notebooks/book/ch*.py and execute it.
# Datasets not shipped in the package come from the data-v1 release (or FDATOOLS_DATA_DIR).
book:
	$(PY) tools/build_book_notebook.py
	$(PY) -m pytest --nbmake notebooks/book_figures.ipynb notebooks/tour.ipynb

# The exact GATE 5 (final) chain from WORKFLOW.md; stops at the first failure.
gate5: export FDATOOLS_DATA_DIR := $(FDATOOLS_DATA_DIR)
gate5:
	$(PY) -m pytest
	$(PY) -m pytest tests/parity
	$(PY) -m mypy --strict src/fdatools
	$(PY) -m ruff check .
	$(PY) -m mkdocs build --strict
	$(PY) -m build && $(PY) -m twine check dist/*
	$(PY) -m pytest --nbmake notebooks/book_figures.ipynb notebooks/tour.ipynb
