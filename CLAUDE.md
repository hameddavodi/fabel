# CLAUDE.md — Fabel Project Contract

You are building **Fabel**: a production-grade Python rewrite of R's `fda` package (Ramsay, v6.3.0).
Authoritative docs: `SPEC.md` (API surface) and `WORKFLOW.md` (execution order). Read both before any work.

## Non-negotiables

1. **Autonomous.** Never ask the user questions. Decision order: `SPEC.md` → R `fda` 6.3.0 behavior → Ramsay/Silverman books → your best judgment (record it in `PROGRESS.md`).
2. **Production quality only.** No TODOs, no stubs, no `NotImplementedError` in shipped modules, no "beta" markers. Every merged module is final-quality.
3. **Parity is law.** Public API output must match R golden files (`rtol=1e-8`; `1e-5` for iterative solvers). Never weaken a tolerance to make a test pass — fix the math.
4. **Clean room.** Never translate R source code line-by-line. Implement from the math in SPEC.md and published references. License: BSD-3.
5. **API freeze.** Public surface = exactly the ~40 symbols in SPEC.md. New helpers go under `fabel._internal`.

## Commands

- Env: `docker compose run dev` (Python 3.11 + R 4.x + rpy2)
- Test: `pytest -x -q` | Parity only: `pytest tests/parity -x`
- Lint/type: `ruff check . && ruff format --check . && mypy --strict src/fabel`
- Golden files: `python tools/make_golden.py <module>`
- Bench: `pytest benchmarks --benchmark-only`
- Docs: `mkdocs build --strict`

## Definition of Done (every task)

- [ ] Tests written FIRST (golden-file parity + unit + property-based where applicable)
- [ ] `pytest` green, coverage of the touched module ≥ 90%
- [ ] `mypy --strict` and `ruff` clean
- [ ] Public functions: full docstrings (NumPy style) with an executable example
- [ ] Benchmark exists if module is on the hot path (basis eval, smoothing, inner products)
- [ ] `PROGRESS.md` updated (status, decisions, failed approaches + why)
- [ ] Conventional commit (`feat(basis): ...`), one commit per green task

## Session protocol

1. Read `PROGRESS.md` first — resume from last unchecked item in `WORKFLOW.md`.
2. Work ONE task at a time. Run the full gate before moving on.
3. If a gate fails 3 attempts in a row: record the failure + hypothesis in `PROGRESS.md`, move to the next independent task, revisit later.
4. Never delete golden files. Never edit files under `tests/golden/`.
5. End every session by updating `PROGRESS.md` and committing.

## Style

- `src/` layout, `pyproject.toml` (hatchling), Python ≥ 3.10
- Array API dispatch in `fabel/_backend.py` — core math must not import numpy directly
- Immutable dataclasses (`frozen=True`) for FData/Basis
- No pandas/torch imports in core modules — optional extras only (`fabel[torch]`, `fabel[pandas]`)
