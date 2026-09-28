"""Tests for ``tools/parity_report.py``, the PARITY_REPORT.md generator."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL = REPO_ROOT / "tools" / "parity_report.py"

_spec = importlib.util.spec_from_file_location("parity_report", TOOL)
assert _spec is not None
assert _spec.loader is not None
pr = importlib.util.module_from_spec(_spec)
sys.modules["parity_report"] = pr
_spec.loader.exec_module(pr)


def _node(outcome: str, **kw: object) -> object:
    node = pr.NodeResult(
        nodeid=str(kw.get("nodeid", "tests/parity/test_x.py::t")),
        module=str(kw.get("module", "basis")),
        case=str(kw.get("case", "c")),
        r_call="",
        symbol=kw.get("symbol", "BSpline"),
    )
    node.outcome = outcome
    node.reason = str(kw.get("reason", ""))
    if "abs_err" in kw:
        node.record(float(kw["abs_err"]), float(kw["rel_err"]), 1e-8)  # type: ignore[arg-type]
    return node


# ------------------------------------------------------------------ mapping


def test_symbol_names_unique_and_rules_point_at_symbols() -> None:
    names = [s.name for s in pr.SYMBOLS]
    assert len(names) == len(set(names))
    assert {rule[2] for rule in pr.CASE_RULES} <= set(names)


def test_every_golden_case_is_attributed() -> None:
    for path in sorted((REPO_ROOT / "tests" / "golden").glob("*.json")):
        doc = json.loads(path.read_text(encoding="utf-8"))
        for case in doc["cases"]:
            text = case.get("r_call") or case["name"]
            assert pr.symbol_for(path.stem, text) is not None, (path.stem, case["name"])


@pytest.mark.parametrize(
    ("module", "text", "expected"),
    [
        ("basis", "eval.penalty(create.bspline.basis(c(0,1), 4, 2), int2Lfd(0))", "BSpline"),
        ("basis", "eval.basis(t, create.polygonal.basis(argvals=c(0,1)), 0)", "Polygonal"),
        ("core", "inprod.bspline(fdobj, fdobj, 0, 0)", "inprod"),
        ("core", "eval.fd(t, fdobj, vec2Lfd(c(0, 0.7), c(0,1)))", "LDO"),
        ("core", "bifd(coef=array(rnorm(24)), s, t)", "BiFData"),
        ("core", "mean.fd(fdobj)", "FData"),
        ("smoothing", "lambda2gcv(-2, 1:365, temp, fdParobj)", "gcv_curve"),
        ("smoothing", "df2lambda(1:365, basis, df=10)", "df_to_lambda"),
        ("smoothing", "smooth.monotone(age, y, fdPar)", "smooth"),
        ("decomposition", "cca.fd(a, b)", "FCCA"),
        ("decomposition", "varmx.pca.fd(pcafd)", "FPCA"),
        ("registration", "landmarkreg(f, m)", "landmark_register"),
        ("registration", "AmpPhaseDecomp(x, y, h)", "register"),
        ("stats", "set.seed(1); Fperm.fd(y, x, b)", "stats.f_test"),
        ("stats", "eval.bifd(grid, grid, var.fd(temp_fd))", "stats.cov"),
        ("io", "test_fd_coefs_and_eval_match_r", "read_rds"),
    ],
)
def test_symbol_for(module: str, text: str, expected: str) -> None:
    assert pr.symbol_for(module, text) == expected


def test_symbol_for_unknown_module() -> None:
    assert pr.symbol_for("nonexistent", "anything") is None


# ------------------------------------------------------------------- errors


def test_errors_basic() -> None:
    assert pr.errors([1.0, 2.0, 4.0], [1.0, 2.5, 4.0]) == (0.5, 0.125)


def test_errors_nan_zero_and_broadcast() -> None:
    abs_err, rel_err = pr.errors([np.nan, 1e-3, 3.0], [np.nan, 0.0, 2.0])
    assert abs_err == 1.0
    assert rel_err == 0.5
    assert pr.errors(3.0, [3.0, 3.0]) == (0.0, 0.0)


def test_errors_edge_cases() -> None:
    assert pr.errors([], []) == (0.0, 0.0)
    assert pr.errors([0.0], [0.0]) == (0.0, 0.0)
    assert pr.errors([1.0], [0.0]) == (1.0, float("inf"))
    assert pr.errors([1.0], [np.inf]) == (0.0, 0.0)
    assert pr.errors(["x"], ["y"]) is None
    assert pr.errors([1.0, 2.0], [1.0, 2.0, 3.0]) is None


def test_record_keeps_maxima() -> None:
    node = _node("pass", abs_err=1e-9, rel_err=1e-10)
    node.record(1e-12, 1e-8, None)  # type: ignore[attr-defined]
    node.record(1e-3, 1e-12, 1e-5)  # type: ignore[attr-defined]
    assert node.abs_err == 1e-3  # type: ignore[attr-defined]
    assert node.rel_err == 1e-8  # type: ignore[attr-defined]
    assert node.rtol == 1e-5  # type: ignore[attr-defined]
    assert node.comparisons == 3  # type: ignore[attr-defined]


# ------------------------------------------------------------------- status


def _symbol(name: str) -> object:
    return next(s for s in pr.SYMBOLS if s.name == name)


def test_symbol_status_branches() -> None:
    bs = _symbol("BSpline")
    assert pr.symbol_status(bs, [_node("pass")], 1) == "pass"
    assert pr.symbol_status(bs, [_node("pass"), _node("xfail"), _node("skip")], 3) == (
        "pass, 1 xfail (R defect), 1 skipped"
    )
    assert pr.symbol_status(bs, [_node("fail"), _node("pass")], 2).startswith("**FAIL**")
    assert pr.symbol_status(bs, [_node("skip")], 1) == "skipped (test data missing)"
    assert pr.symbol_status(bs, [_node("xfail")], 1) == "all 1 checks xfail (R defect)"
    assert pr.symbol_status(bs, [], 4) == "golden cases exist, no parity test yet"
    assert pr.symbol_status(bs, [], 0) == "no golden cases"
    ghost = pr.Symbol("ghost", "fdatools.no_such_module", "X", "none")
    assert pr.symbol_status(ghost, [], 3) == "not built on this branch"
    missing_attr = pr.Symbol("ghost", "fdatools.core", "NoSuchThing", "none")
    assert not pr.available(missing_attr)


def test_formatting_helpers() -> None:
    assert pr._fmt(None) == "—"
    assert pr._fmt(0.0) == "0"
    assert pr._fmt(1.5e-9) == "1.50e-09"
    assert pr._fmt(float("inf")) == "inf"
    assert pr._cell("a |\n b") == "a \\| b"


def test_render_sections() -> None:
    results = [
        _node("pass", case="ok", abs_err=1e-12, rel_err=1e-13),
        _node("xfail", case="rdefect", reason="R is | wrong\nhere", abs_err=1.0, rel_err=0.1),
        _node("fail", case="broken", nodeid="tests/parity/test_x.py::broken"),
        _node("skip", case="nodata", reason="dump missing", symbol=None),
    ]
    text = pr.render(results, ["tests/parity"])
    assert text.startswith("# Parity report")
    assert "| `BSpline` |" in text
    assert "1.00e-12" in text
    assert "`rdefect` (basis, rel err 1.00e-01): R is \\| wrong here" in text
    assert "`tests/parity/test_x.py::broken` -- fail" in text
    assert "(unattributed)" in text
    assert "**FAIL**" in text
    assert "not built on this branch" in text or "no parity test yet" in text


def test_render_without_xfails() -> None:
    text = pr.render([_node("pass", abs_err=0.0, rel_err=0.0)], ["x"])
    assert "## Strict xfails (R fda is the less accurate side)\n\nNone." in text
    assert "## Failing or skipped checks" not in text


# --------------------------------------------------------------- end to end


def _run(tmp_path: Path, *tests: str) -> tuple[int, str]:
    out = tmp_path / "REPORT.md"
    env = dict(os.environ, PYTHONPATH=str(REPO_ROOT / "src"))
    proc = subprocess.run(
        [sys.executable, str(TOOL), "--output", str(out), "--tests", *tests],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=600,
    )
    return proc.returncode, out.read_text(encoding="utf-8")


@pytest.mark.slow
def test_end_to_end_on_real_parity_tests(tmp_path: Path) -> None:
    pytest.importorskip("rdata")
    code, text = _run(tmp_path, "tests/parity/test_io.py")
    assert code == 0
    assert "| io | `test_fd_coefs_and_eval_match_r` | `read_rds` | 1.00e-08 |" in text
    assert "2 parity checks -- 2 pass" in text


@pytest.mark.slow
def test_end_to_end_reports_failures_and_xfails(tmp_path: Path) -> None:
    suite = tmp_path / "test_fake.py"
    suite.write_text(
        "import numpy as np\n"
        "import pytest\n\n"
        "def test_bad():\n"
        "    np.testing.assert_allclose([1.0], [2.0], rtol=1e-8)\n\n"
        "@pytest.mark.xfail(strict=True, reason='R is off by 50%')\n"
        "def test_known():\n"
        "    np.testing.assert_allclose([3.0], [2.0], rtol=1e-8)\n\n"
        "def test_good():\n"
        "    np.testing.assert_allclose([2.0], [2.0])\n",
        encoding="utf-8",
    )
    code, text = _run(tmp_path, str(suite))
    assert code == 1
    assert "3 parity checks -- 1 pass, 1 strict xfail (R defect), 1 fail" in text
    assert "R is off by 50%" in text
    assert "::test_bad` -- fail" in text


# ------------------------------------------------------- recorder, in-process


class _CallSpec:
    def __init__(self, params: dict[str, object]) -> None:
        self.params = params


class _Item:
    def __init__(self, nodeid: str, name: str, params: dict[str, object] | None = None) -> None:
        self.nodeid = nodeid
        self.name = name
        self.fspath = nodeid.split("::")[0]
        if params is not None:
            self.callspec = _CallSpec(params)


class _Report:
    def __init__(self, nodeid: str, when: str, outcome: str, **extra: object) -> None:
        self.nodeid = nodeid
        self.when = when
        self.passed = outcome == "passed"
        self.failed = outcome == "failed"
        self.skipped = outcome == "skipped"
        self.longrepr = extra.pop("longrepr", None)
        for key, value in extra.items():
            setattr(self, key, value)


SMOOTH_CASE = {"name": "smooth_basis_x", "r_call": "smooth.basis(t, y, fdPar(b, 2, 1))"}
LAMBDA_CASE = {"name": "lambda2df_x", "r_call": "lambda2df(t, basis, lambda=1)"}


def _drive(recorder: Any, items: list[_Item], bodies: Mapping[str, Callable[[], object]]) -> None:
    """Replay pytest's hook sequence on the recorder with fake items/reports."""
    recorder.pytest_collection_modifyitems(items)
    for item in items:
        gen = recorder.pytest_runtest_call(item)
        next(gen)
        outcome = "passed"
        try:
            bodies[item.nodeid]()
        except AssertionError:
            outcome = "failed"
        with pytest.raises(StopIteration):
            next(gen)
        report = _Report(item.nodeid, "call", outcome)
        recorder.pytest_runtest_logreport(report)


def test_recorder_measures_and_attributes() -> None:
    recorder = pr.Recorder()
    items = [
        _Item(
            "tests/parity/test_smoothing.py::t[a-coefs]",
            "t[smooth_basis_x-coefs]",
            {"case": SMOOTH_CASE, "field": "coefs"},
        ),
        _Item("tests/parity/test_smoothing.py::t[b]", "t[lambda2df_x]", {"case": LAMBDA_CASE}),
        _Item("tests/parity/test_io.py::test_fd", "test_fd"),
        _Item("tests/parity/test_datasets.py::test_growth", "test_growth"),
    ]
    original = np.testing.assert_allclose
    bodies = {
        items[0].nodeid: lambda: np.testing.assert_allclose([1.0, 2.0], [1.0, 2.0], 1e-8),
        items[1].nodeid: lambda: np.testing.assert_allclose([1.0], [1.5], rtol=1e-8),
        items[2].nodeid: lambda: np.testing.assert_allclose([1.0, 2.0], [1.0, 2.0, 3.0]),
        items[3].nodeid: lambda: None,
    }
    _drive(recorder, items, bodies)
    assert np.testing.assert_allclose is original
    res = recorder.results
    first = res[items[0].nodeid]
    assert (first.case, first.symbol, first.outcome, first.rtol) == (
        "smooth_basis_x-coefs",
        "smooth",
        "pass",
        1e-8,
    )
    second = res[items[1].nodeid]
    assert (second.case, second.symbol, second.outcome) == ("lambda2df_x", "lambda_to_df", "fail")
    assert second.abs_err == 0.5
    third = res[items[2].nodeid]
    assert (third.symbol, third.comparisons, third.outcome) == ("read_rds", 0, "fail")
    assert res[items[3].nodeid].symbol == "datasets.load_*"


def test_recorder_outcomes_from_reports() -> None:
    recorder = pr.Recorder()
    ids = [f"tests/parity/test_core.py::t{i}" for i in range(4)]
    recorder.pytest_collection_modifyitems([_Item(i, i.split("::")[1]) for i in ids])
    log = recorder.pytest_runtest_logreport
    log(_Report(ids[0], "call", "skipped", wasxfail="reason: R off"))
    log(_Report(ids[1], "setup", "skipped", longrepr=("f.py", 1, "Skipped: no dump")))
    log(_Report(ids[2], "setup", "skipped", longrepr="plain"))
    log(_Report(ids[3], "setup", "passed"))
    log(_Report("not/collected.py::x", "call", "passed"))
    res = recorder.results
    assert (res[ids[0]].outcome, res[ids[0]].reason) == ("xfail", "R off")
    assert res[ids[1]].reason == "Skipped: no dump"
    assert res[ids[2]].reason == "plain"
    assert res[ids[3]].outcome == "unknown"
    assert "not/collected.py::x" not in res


def test_main_writes_report_and_sets_exit_status(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[list[str]] = []

    def fake_main(args: list[str], plugins: list[object]) -> int:
        seen.append(args)
        item = _Item("tests/parity/test_io.py::test_x", "test_x")
        _drive(plugins[0], [item], {item.nodeid: lambda: np.testing.assert_allclose(1.0, 1.0)})
        return 0

    monkeypatch.setattr(pr.pytest, "main", fake_main)
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "R.md"
    assert pr.main(["--output", str(out), "--tests", "tests/parity/test_io.py"]) == 0
    assert seen[0][-1] == "tests/parity/test_io.py"
    assert "1 parity checks -- 1 pass" in out.read_text(encoding="utf-8")

    monkeypatch.setattr(pr.pytest, "main", lambda args, plugins: 5)
    assert pr.main(["--output", str(out)]) == 1
    assert "0 parity checks" in out.read_text(encoding="utf-8")
