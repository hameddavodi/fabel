# PROGRESS
## Status: Phase 0, task 1
## Done
## Decisions
- 2026-09-06 Python env managed with `uv` (AISMA §5.1: one package manager). Local dev on Python 3.12; Docker image on 3.11 per CLAUDE.md.
- 2026-09-06 Golden files: `tools/make_golden.py` drives R through rpy2 when importable, else falls back to `Rscript` subprocess with the same seeded R scripts. Output format identical.
- 2026-09-06 Commit style: project CLAUDE.md conventional commits (`feat(basis): ...`) win over AISMA prefix scheme (project file overrides global, AISMA rule 1: adopt repo conventions).
## Failed approaches (do not retry)
## Blocked
