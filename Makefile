.PHONY: lint type test parity golden bench docs build book

lint:
	uv run ruff check .
	uv run ruff format --check .

type:
	uv run mypy src tests tools benchmarks

test:
	uv run pytest tests/unit -m "not parity and not slow and not gpu" --cov=fabel --cov-report=term-missing

parity:
	uv run pytest tests/unit/test_golden_schema.py tests/parity -m parity

golden:
	uv run python tools/make_golden.py basis
	uv run python tools/make_golden.py core

bench:
	uv run pytest benchmarks --benchmark-only

docs:
	uv run mkdocs build --strict

build:
	uv run python -m build
	uv run twine check dist/*

# Rebuild notebooks/book_figures.ipynb from notebooks/book/ch*.py and execute it.
# Datasets not shipped in the package come from the data-v1 release (or FABEL_DATA_DIR).
book:
	uv run python tools/build_book_notebook.py
	uv run pytest --nbmake notebooks/book_figures.ipynb
