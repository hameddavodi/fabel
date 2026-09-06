# PROGRESS
## Status: Phase 0 harness + Phase 1 golden files done (basis, core)
## Done
- `tools/make_golden.py` + `tools/golden_r/{common,basis,core}.R`: golden-file generator, rpy2-first with `Rscript` fallback.
- `tests/golden/basis.json` (256 cases), `tests/golden/core.json` (28 cases).
- `tests/parity/__init__.py`, `tests/parity/conftest.py` (`golden`/`golden_cases` helpers).
- `tests/unit/test_golden_schema.py`: schema validation for every `tests/golden/*.json`.
- `.github/workflows/ci.yml`, `.github/dependabot.yml`, `.github/PULL_REQUEST_TEMPLATE.md`, `Makefile`.
## Decisions
- 2026-09-06 Python env managed with `uv` (AISMA §5.1: one package manager). Local dev on Python 3.12; Docker image on 3.11 per CLAUDE.md.
- 2026-09-06 Golden files: `tools/make_golden.py` drives R through rpy2 when importable, else falls back to `Rscript` subprocess with the same seeded R scripts. Output format identical.
- 2026-09-06 Commit style: project CLAUDE.md conventional commits (`feat(basis): ...`) win over AISMA prefix scheme (project file overrides global, AISMA rule 1: adopt repo conventions).
- 2026-09-06 `getbasispenalty()` throws "Basis type not recognizable" for monomial and constant basis types in fda 6.3.0 (all nbasis tested); `eval.penalty()` works uniformly across all 7 basis types. `tools/golden_r/basis.R` uses `eval.penalty()` exclusively for consistency.
- 2026-09-06 B-spline penalty validity constraint empirically determined to be `Lfdobj <= norder - 2` (not the naive `Lfdobj < norder`); `getbasispenalty`/`eval.penalty` errors otherwise (e.g. order 3, Lfdobj 2 fails). Verified by systematic probing over norder 1-6 x Lfdobj 0-5.
- 2026-09-06 fda 6.3.0 exposes `sd.fd`, `stddev.fd`, `std.fd` as equivalent entry points (confirmed via `exists()`); `tools/golden_r/core.R` standardizes on `sd.fd`.
- 2026-09-06 `core.json` rtol set to 1e-6 (vs 1e-8 for `basis.json`): `inprod()` uses Richardson extrapolation, "good to about four to five significant digits" per fda's own docs.
- 2026-09-06 fda 6.3.0's `bifd()` is unconditionally broken for a 3-D coef array (`ndim == 3`, reps dim with no vars dim): its body sets `defaultnames` only in `if (ndim == 2)`/`if (ndim == 4)` branches (copy-paste bug omits `ndim == 3`) but then unconditionally runs `names(defaultnames) <- c(...)`, throwing "object 'defaultnames' not found" even when `fdnames` is passed explicitly (that line never references the `fdnames` arg). Worked around in `tools/golden_r/core.R` by using a 4-D array (`dim = c(nbasis_s, nbasis_t, nrep, nvar=1)`) instead, which hits the working `ndim == 4` branch and still exercises the reps dimension.
- 2026-09-06 GitHub Actions in `.github/workflows/ci.yml` pinned to full commit SHAs verified via `gh api repos/<owner>/<repo>/tags` (actions/checkout v7.0.1, astral-sh/setup-uv v10.0.1, actions/upload-artifact v7.0.1), each with a `# vX.Y.Z` comment.
- 2026-09-06 `parity` CI job runs on plain Python (uv-managed), no R/Docker: golden JSON files are pre-generated and committed, so parity tests only need `fabel` + the JSON, not R/rpy2. R is only needed to regenerate golden files via `make golden` / `tools/make_golden.py`.
## Failed approaches (do not retry)
- Naive B-spline penalty Lfdobj bound `Lfdobj < norder` (from initial task heuristic): fails for e.g. order 3, Lfdobj 2. Correct bound is `Lfdobj <= norder - 2`.
- Passing `fdnames` explicitly to `bifd()` to work around its 3-D `defaultnames` bug: does not help, since the crashing line (`names(defaultnames) <- ...`) never reads the `fdnames` argument. Must avoid ndim==3 arrays entirely (use ndim==2 or ndim==4).
## Blocked
