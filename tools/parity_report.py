"""Generate ``PARITY_REPORT.md``: Fabel's numeric agreement with R ``fda`` 6.3.0.

Usage::

    PYTHONPATH=src .venv/bin/python tools/parity_report.py [--output PARITY_REPORT.md]
                                                          [--tests tests/parity]

How it works
------------
The report is measured, not recomputed by hand.  The tool runs the golden-file
parity suite (``tests/parity`` by default) in-process through ``pytest.main``
with a small plugin attached.  While a test body runs, the plugin wraps
``numpy.testing.assert_allclose`` -- the one comparison every parity test goes
through -- and records, for each comparison, the largest absolute error
``max |actual - expected|`` and the normwise relative error
``max |actual - expected| / max |expected|`` (both over the finite expected
entries).  The normwise form is used on purpose: an entrywise ratio divides by
matrix entries that are zero up to rounding (e.g. ``1e-17`` in a penalty matrix
whose largest entry is ``1e6``) and reports meaningless numbers, while the
parity tests themselves absorb those entries with an ``atol`` scaled by
``max |expected|``.  It also records the tolerance the test asked for and the
test outcome (pass, strict xfail with its written reason, fail, skip).

Each test is then attributed to one public Fabel symbol: parametrised parity
tests carry their golden case (``callspec.params["case"]``), whose ``r_call``
names the R function under test; unparametrised tests are attributed by their
test module.  Golden files with no parity test yet on the current branch are
still listed, with their case count, so the report shows every public symbol
whether or not its module exists.

Because the numbers come from the real test run, the report stays correct as
new parity tests are added: nothing in this file needs to change unless a new
public symbol or a new R function family appears (see :data:`SYMBOLS` and
:data:`CASE_RULES`).

Exit status is ``0`` when every collected parity test passed or failed as a
documented strict xfail, ``1`` otherwise (the report is written either way).
"""

from __future__ import annotations

import argparse
import importlib
import json
import math
import os
import re
import sys
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import numpy.testing
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
GOLDEN_DIR = REPO_ROOT / "tests" / "golden"


@dataclass(frozen=True)
class Symbol:
    """One public Fabel symbol and the R ``fda`` functions it replaces."""

    name: str
    module: str
    attr: str
    r_counterpart: str


#: Every public symbol of SPEC.md, in the order the report lists them.
SYMBOLS: tuple[Symbol, ...] = (
    Symbol(
        "FData", "fabel.core", "FData", "fd, eval.fd, mean.fd, sd.fd, center.fd, deriv.fd, +.fd"
    ),
    Symbol("BiFData", "fabel.core", "BiFData", "bifd, eval.bifd"),
    Symbol("LDO", "fabel.core", "LDO", "Lfd, int2Lfd, vec2Lfd"),
    Symbol("inprod", "fabel.core", "inprod", "inprod, inprod.bspline"),
    Symbol("Basis", "fabel.basis", "Basis", "basisfd"),
    Symbol("BSpline", "fabel.basis", "BSpline", "create.bspline.basis, bsplineS, bsplinepen"),
    Symbol("Fourier", "fabel.basis", "Fourier", "create.fourier.basis, fourier, fourierpen"),
    Symbol("Monomial", "fabel.basis", "Monomial", "create.monomial.basis, monomial"),
    Symbol("Exponential", "fabel.basis", "Exponential", "create.exponential.basis, expon"),
    Symbol("Power", "fabel.basis", "Power", "create.power.basis, powerbasis, powerpen"),
    Symbol("Constant", "fabel.basis", "Constant", "create.constant.basis"),
    Symbol("Polygonal", "fabel.basis", "Polygonal", "create.polygonal.basis, polyg, polygpen"),
    Symbol(
        "smooth",
        "fabel.smoothing",
        "smooth",
        "smooth.basis, smooth.basisPar, Data2fd, smooth.monotone, smooth.pos",
    ),
    Symbol("SmoothResult", "fabel.smoothing", "SmoothResult", "smooth.basis return list"),
    Symbol("Smoother", "fabel.smoothing", "Smoother", "(new: sklearn estimator)"),
    Symbol("gcv_curve", "fabel.smoothing", "gcv_curve", "lambda2gcv"),
    Symbol("lambda_to_df", "fabel.smoothing", "lambda_to_df", "lambda2df"),
    Symbol("df_to_lambda", "fabel.smoothing", "df_to_lambda", "df2lambda"),
    Symbol("FPCA", "fabel.decomposition", "FPCA", "pca.fd, varmx.pca.fd"),
    Symbol("FCCA", "fabel.decomposition", "FCCA", "cca.fd"),
    Symbol(
        "fregress",
        "fabel.regression",
        "fregress",
        "fRegress, predict.fRegress, fRegress.stderr, fRegress.CV",
    ),
    Symbol("FRegress", "fabel.regression", "FRegress", "fRegress (estimator form)"),
    Symbol("register", "fabel.registration", "register", "register.fd, AmpPhaseDecomp"),
    Symbol("landmark_register", "fabel.registration", "landmark_register", "landmarkreg"),
    Symbol("Registrator", "fabel.registration", "Registrator", "register.fd (estimator form)"),
    Symbol("PDA", "fabel.dynamics", "PDA", "pda.fd, pda.overlay"),
    Symbol("phase_plane", "fabel.dynamics", "phase_plane", "phaseplanePlot"),
    Symbol("stats.cov", "fabel.stats", "cov", "var.fd"),
    Symbol("stats.cor", "fabel.stats", "cor", "cor.fd"),
    Symbol("stats.depth", "fabel.stats", "depth", "fdepth"),
    Symbol("stats.boxplot", "fabel.stats", "boxplot", "fbplot, boxplot.fd"),
    Symbol("stats.f_test", "fabel.stats", "f_test", "Fperm.fd"),
    Symbol("stats.t_test", "fabel.stats", "t_test", "tperm.fd"),
    Symbol("datasets.load_*", "fabel.datasets", "load_growth", "data(package='fda')"),
    Symbol("nn.BasisLayer", "fabel.nn", "BasisLayer", "(new: PyTorch layer)"),
    Symbol("nn.FDataDataset", "fabel.nn", "FDataDataset", "(new: PyTorch dataset)"),
    Symbol("from_pandas", "fabel.io", "from_pandas", "(new)"),
    Symbol("to_pandas", "fabel.io", "to_pandas", "(new)"),
    Symbol("to_xarray", "fabel.io", "to_xarray", "(new)"),
    Symbol("read_rds", "fabel.io", "read_rds", "readRDS on fd / bifd / basisfd objects"),
)

#: ``(golden module, regex searched in the case's r_call or name, symbol name)``.
#: The first matching rule wins, so specific patterns come before general ones.
CASE_RULES: tuple[tuple[str, str, str], ...] = (
    ("basis", r"create\.bspline\.basis", "BSpline"),
    ("basis", r"create\.fourier\.basis", "Fourier"),
    ("basis", r"create\.monomial\.basis", "Monomial"),
    ("basis", r"create\.exponential\.basis", "Exponential"),
    ("basis", r"create\.power\.basis", "Power"),
    ("basis", r"create\.constant\.basis", "Constant"),
    ("basis", r"create\.polygonal\.basis", "Polygonal"),
    ("core", r"\binprod", "inprod"),
    ("core", r"\bbifd\(", "BiFData"),
    ("core", r"eval\.fd\(t, fdobj, (int2Lfd|vec2Lfd)", "LDO"),
    ("core", r".", "FData"),
    ("smoothing", r"lambda2gcv", "gcv_curve"),
    ("smoothing", r"lambda2df", "lambda_to_df"),
    ("smoothing", r"df2lambda", "df_to_lambda"),
    ("smoothing", r".", "smooth"),
    ("decomposition", r"cca\.fd", "FCCA"),
    ("decomposition", r".", "FPCA"),
    ("regression", r".", "fregress"),
    ("registration", r"landmarkreg", "landmark_register"),
    ("registration", r".", "register"),
    ("dynamics", r".", "PDA"),
    ("stats", r"cor\.fd", "stats.cor"),
    ("stats", r"fbplot", "stats.boxplot"),
    ("stats", r"fdepth", "stats.depth"),
    ("stats", r"tperm\.fd", "stats.t_test"),
    ("stats", r"Fperm\.fd", "stats.f_test"),
    ("stats", r"var\.fd", "stats.cov"),
    ("stats", r".", "FData"),
    ("datasets", r".", "datasets.load_*"),
    ("io", r".", "read_rds"),
)


def symbol_for(module: str, text: str) -> str | None:
    """Return the public symbol a golden case exercises.

    Parameters
    ----------
    module : str
        Golden module name (``tests/golden/<module>.json``) or parity test module
        stem without the ``test_`` prefix.
    text : str
        The case's ``r_call`` (preferred) or its name.

    Returns
    -------
    str or None
        The symbol name from :data:`SYMBOLS`, or ``None`` when no rule matches.

    Examples
    --------
    >>> symbol_for("basis", "eval.basis(t, create.fourier.basis(c(0,1), 3, 1), 0)")
    'Fourier'
    >>> symbol_for("stats", "fdepth(list(argvals=grid, y=Y))")
    'stats.depth'
    """
    for rule_module, pattern, name in CASE_RULES:
        if rule_module == module and re.search(pattern, text):
            return name
    return None


def errors(actual: Any, expected: Any) -> tuple[float, float] | None:
    """Return ``(max absolute error, normwise relative error)`` of one comparison.

    Both are taken over the entries whose expected value is finite; entries that
    are NaN on both sides are ignored.  The relative error is
    ``max |actual - expected| / max |expected|``; it is ``inf`` when the expected
    values are all zero but the actual ones are not.  Returns ``None`` when the
    inputs are not numeric or their shapes cannot be broadcast.

    Examples
    --------
    >>> errors([1.0, 2.0], [1.0, 2.5])
    (0.5, 0.2)
    >>> errors(["a"], ["b"]) is None
    True
    """
    try:
        a = np.asarray(actual, dtype=float)
        e = np.asarray(expected, dtype=float)
        a, e = np.broadcast_arrays(a, e)
    except (TypeError, ValueError):
        return None
    finite = np.isfinite(e)
    if not finite.any():
        return 0.0, 0.0
    diff = np.abs(a[finite] - e[finite])
    abs_err = float(np.max(diff))
    scale = float(np.max(np.abs(e[finite])))
    if scale == 0.0:
        return abs_err, (math.inf if abs_err else 0.0)
    return abs_err, abs_err / scale


@dataclass
class NodeResult:
    """Measured outcome of one parity test node."""

    nodeid: str
    module: str
    case: str
    r_call: str
    symbol: str | None
    outcome: str = "unknown"
    reason: str = ""
    abs_err: float | None = None
    rel_err: float | None = None
    rtol: float | None = None
    comparisons: int = 0

    def record(self, abs_err: float, rel_err: float, rtol: float | None) -> None:
        """Fold one comparison into the node's running maxima."""
        self.comparisons += 1
        self.abs_err = abs_err if self.abs_err is None else max(self.abs_err, abs_err)
        self.rel_err = rel_err if self.rel_err is None else max(self.rel_err, rel_err)
        if rtol is not None:
            self.rtol = rtol if self.rtol is None else max(self.rtol, rtol)


@dataclass
class Recorder:
    """pytest plugin: measures every ``assert_allclose`` made by a parity test."""

    results: dict[str, NodeResult] = field(default_factory=dict)
    _current: NodeResult | None = None

    @staticmethod
    def _describe(item: pytest.Item) -> NodeResult:
        module = Path(str(item.fspath)).stem.removeprefix("test_")
        callspec = getattr(item, "callspec", None)
        case = callspec.params.get("case") if callspec is not None else None
        if isinstance(case, dict):
            name = str(case.get("name", item.name))
            r_call = str(case.get("r_call", ""))
            suffix = item.name.split("[", 1)[1].rstrip("]") if "[" in item.name else name
            label = suffix if suffix.startswith(name) else name
        else:
            label, r_call = item.name, ""
        symbol = symbol_for(module, r_call or label)
        return NodeResult(item.nodeid, module, label, r_call, symbol)

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_call(self, item: pytest.Item) -> Iterator[None]:
        """Wrap ``numpy.testing.assert_allclose`` while the test body runs."""
        node = self.results.setdefault(item.nodeid, self._describe(item))
        original: Callable[..., None] = numpy.testing.assert_allclose

        def measured(actual: Any, desired: Any, *args: Any, **kwargs: Any) -> None:
            measured_err = errors(actual, desired)
            rtol = kwargs.get("rtol", args[0] if args else None)
            if measured_err is not None:
                node.record(*measured_err, float(rtol) if rtol is not None else None)
            original(actual, desired, *args, **kwargs)

        numpy.testing.assert_allclose = measured
        try:
            yield
        finally:
            numpy.testing.assert_allclose = original

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        """Register every collected node, so skipped ones are reported too."""
        for item in items:
            self.results.setdefault(item.nodeid, self._describe(item))

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        """Record the final outcome of each node."""
        node = self.results.get(report.nodeid)
        if node is None:
            return
        if hasattr(report, "wasxfail"):
            node.outcome, node.reason = "xfail", str(report.wasxfail).removeprefix("reason: ")
        elif report.failed:
            node.outcome = "fail"
        elif report.skipped:
            longrepr = report.longrepr
            node.outcome = "skip"
            node.reason = str(longrepr[2]) if isinstance(longrepr, tuple) else str(longrepr)
        elif report.when == "call" and report.passed:
            node.outcome = "pass"


def golden_inventory() -> dict[str, dict[str, Any]]:
    """Return ``{module: {"meta": ..., "per_symbol": {symbol: n_cases}}}`` for every golden file.

    Examples
    --------
    >>> inv = golden_inventory()
    >>> inv["basis"]["per_symbol"]["BSpline"] > 0
    True
    """
    inventory: dict[str, dict[str, Any]] = {}
    for path in sorted(GOLDEN_DIR.glob("*.json")):
        doc = json.loads(path.read_text(encoding="utf-8"))
        module = path.stem
        per_symbol: dict[str, int] = {}
        for case in doc["cases"]:
            symbol = symbol_for(module, str(case.get("r_call") or case["name"]))
            if symbol is not None:
                per_symbol[symbol] = per_symbol.get(symbol, 0) + 1
        inventory[module] = {"meta": doc["meta"], "per_symbol": per_symbol}
    return inventory


def available(symbol: Symbol) -> bool:
    """Return whether ``symbol`` can be imported on the current branch."""
    try:
        module = importlib.import_module(symbol.module)
    except ImportError:
        return False
    return hasattr(module, symbol.attr)


def _fmt(value: float | None) -> str:
    if value is None:
        return "—"
    if value == 0.0:
        return "0"
    if math.isnan(value) or math.isinf(value):
        return str(value)
    return f"{value:.2e}"


def _cell(text: str) -> str:
    return " ".join(text.split()).replace("|", "\\|")


def symbol_status(symbol: Symbol, nodes: Sequence[NodeResult], golden_cases: int) -> str:
    """Summarise one symbol's parity status in a few words.

    Examples
    --------
    >>> s = Symbol("X", "fabel.core", "FData", "x")
    >>> symbol_status(s, [], 0)
    'no golden cases'
    """
    counts = {o: sum(n.outcome == o for n in nodes) for o in ("pass", "xfail", "fail", "skip")}
    if counts["fail"]:
        return f"**FAIL** ({counts['fail']} failing)"
    if nodes and counts["skip"] == len(nodes):
        return "skipped (test data missing)"
    if counts["xfail"] and not counts["pass"]:
        return f"all {counts['xfail']} checks xfail (R defect)"
    if counts["pass"]:
        extra = f", {counts['xfail']} xfail (R defect)" if counts["xfail"] else ""
        skipped = f", {counts['skip']} skipped" if counts["skip"] else ""
        return f"pass{extra}{skipped}"
    if not available(symbol):
        return "not built on this branch"
    if golden_cases:
        return "golden cases exist, no parity test yet"
    return "no golden cases"


def render(results: Sequence[NodeResult], test_paths: Sequence[str]) -> str:
    """Render the Markdown report from measured node results."""
    inventory = golden_inventory()
    golden_per_symbol: dict[str, int] = {}
    for entry in inventory.values():
        for name, n in entry["per_symbol"].items():
            golden_per_symbol[name] = golden_per_symbol.get(name, 0) + n
    by_symbol: dict[str, list[NodeResult]] = {}
    for node in results:
        by_symbol.setdefault(node.symbol or "(unattributed)", []).append(node)

    fabel = importlib.import_module("fabel")
    metas = [entry["meta"] for entry in inventory.values()]
    r_version = metas[0]["r_version"] if metas else "unknown"
    fda_version = metas[0]["fda_version"] if metas else "unknown"
    totals = {o: sum(n.outcome == o for n in results) for o in ("pass", "xfail", "fail", "skip")}

    lines = [
        "# Parity report",
        "",
        f"Fabel {fabel.__version__} against R `fda` {fda_version} ({r_version}).",
        f"Generated {date.today().isoformat()} by `tools/parity_report.py` "
        f"from a live run of `{' '.join(test_paths)}`.",
        "",
        "Errors are measured on every `assert_allclose` a parity test makes:",
        "absolute error is `max |fabel - R|`, relative error is the normwise",
        "`max |fabel - R| / max |R|` (entrywise ratios are meaningless on matrix entries",
        "that are zero up to rounding). The summary columns show",
        "the worst value over **passing** cases only; strict xfails (cases where R is",
        "demonstrably the less accurate side) are listed with their measured reason below.",
        "",
        f"**Totals:** {len(results)} parity checks -- {totals['pass']} pass, "
        f"{totals['xfail']} strict xfail (R defect), {totals['fail']} fail, "
        f"{totals['skip']} skipped.",
        "",
        "## By public symbol",
        "",
        "| Symbol | R counterpart | Golden cases | Checks | Max abs err | Max rel err | Status |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for symbol in SYMBOLS:
        nodes = by_symbol.get(symbol.name, [])
        passing = [n for n in nodes if n.outcome == "pass"]
        abs_errs = [n.abs_err for n in passing if n.abs_err is not None]
        rel_errs = [n.rel_err for n in passing if n.rel_err is not None]
        lines.append(
            f"| `{symbol.name}` | {_cell(symbol.r_counterpart)} "
            f"| {golden_per_symbol.get(symbol.name, 0)} | {len(nodes)} "
            f"| {_fmt(max(abs_errs) if abs_errs else None)} "
            f"| {_fmt(max(rel_errs) if rel_errs else None)} "
            f"| {symbol_status(symbol, nodes, golden_per_symbol.get(symbol.name, 0))} |"
        )

    xfails = [n for n in results if n.outcome == "xfail"]
    lines += ["", "## Strict xfails (R fda is the less accurate side)", ""]
    if not xfails:
        lines.append("None.")
    for node in xfails:
        lines.append(
            f"- `{node.case}` ({node.module}, rel err {_fmt(node.rel_err)}): {_cell(node.reason)}"
        )

    problems = [n for n in results if n.outcome in ("fail", "skip", "unknown")]
    if problems:
        lines += ["", "## Failing or skipped checks", ""]
        for node in problems:
            reason = f": {_cell(node.reason)}" if node.reason else ""
            lines.append(f"- `{node.nodeid}` -- {node.outcome}{reason}")

    lines += [
        "",
        "## Every check",
        "",
        "| Module | Case | Symbol | rtol | Max abs err | Max rel err | Status |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for node in results:
        lines.append(
            f"| {node.module} | `{node.case}` | `{node.symbol or '(unattributed)'}` "
            f"| {_fmt(node.rtol)} | {_fmt(node.abs_err)} | {_fmt(node.rel_err)} | {node.outcome} |"
        )
    return "\n".join(lines) + "\n"


def run(test_paths: Sequence[str], output: Path) -> int:
    """Run the parity tests, write the report and return the exit status.

    Parameters
    ----------
    test_paths : sequence of str
        Test files or directories handed to pytest.
    output : Path
        Where the Markdown report is written.

    Returns
    -------
    int
        ``0`` when every check passed or is a documented strict xfail, else ``1``.
    """
    recorder = Recorder()
    os.chdir(REPO_ROOT)
    pytest.main(["-q", "-p", "no:cacheprovider", "--no-header", *test_paths], plugins=[recorder])
    results = list(recorder.results.values())
    output.write_text(render(results, test_paths), encoding="utf-8")
    bad = [n for n in results if n.outcome in ("fail", "unknown")]
    return 1 if bad or not results else 0


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--output", type=Path, default=REPO_ROOT / "PARITY_REPORT.md")
    parser.add_argument("--tests", nargs="+", default=["tests/parity"])
    args = parser.parse_args(argv)
    return run(args.tests, args.output.resolve())


if __name__ == "__main__":
    sys.exit(main())
