"""Parity of :mod:`fabel.io.read_rds` against values recorded straight from R.

Fixtures live at ``tests/fixtures/*.rds`` (written by R's ``saveRDS()``) with a
companion ``*_expected.json`` recording R's own ``eval.fd``/``eval.bifd``
output at ``digits=17``, generated the same way as every other golden file
(see ``docs/dev/conventions.md``); unlike ``tests/golden/*.json`` these are
one-off I/O fixtures, not replayed against ``tools/golden_r/``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import fabel as fb

pytest.importorskip("rdata", reason="io extra (rdata) not installed")

pytestmark = pytest.mark.parity

FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"
RTOL = 1e-8


def _expected(name: str) -> dict[str, Any]:
    with (FIXTURES / name).open(encoding="utf-8") as fh:
        loaded: dict[str, Any] = json.load(fh)
    return loaded


def test_fd_coefs_and_eval_match_r() -> None:
    fdobj = fb.read_rds(FIXTURES / "bspline_fd.rds")
    assert isinstance(fdobj, fb.FData)
    expected = _expected("bspline_fd_expected.json")

    np.testing.assert_allclose(fdobj.coefs, expected["coefs"], rtol=RTOL, atol=1e-12)
    evalt = np.array(expected["evalt"])
    np.testing.assert_allclose(fdobj(evalt), np.array(expected["evalvals"]), rtol=RTOL, atol=1e-12)
    assert tuple(fdobj.basis.domain) == tuple(expected["domain"])
    assert fdobj.basis.n_basis == np.asarray(expected["nbasis"]).item()


def test_bifd_coefs_and_eval_match_r() -> None:
    bifdobj = fb.read_rds(FIXTURES / "bspline_bifd.rds")
    assert isinstance(bifdobj, fb.BiFData)
    expected = _expected("bspline_bifd_expected.json")

    np.testing.assert_allclose(bifdobj.coefs, expected["coef"], rtol=RTOL, atol=1e-12)
    s = np.array(expected["evals"])
    t = np.array(expected["evalt"])
    np.testing.assert_allclose(bifdobj(s, t), np.array(expected["evalvals"]), rtol=RTOL, atol=1e-12)
    assert tuple(bifdobj.sbasis.domain) == tuple(expected["sdomain"])
    assert tuple(bifdobj.tbasis.domain) == tuple(expected["tdomain"])
