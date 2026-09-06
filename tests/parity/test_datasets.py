"""Parity of :mod:`fabel.datasets` against golden output from R ``fda`` 6.3.0.

Every golden case in ``tests/golden/datasets.json`` covers one of the ten
datasets exported by the original ``tools/golden_r/datasets.R`` (see
``PROGRESS.md`` for the decision to scope strict golden parity to exactly
those ten; the four extra dataset families added later --
MontrealTemp/daily/infantGrowth/nondurables/lip -- only have unit tests).

Compared at ``rtol=1e-12``: the golden values were exported from R at full
(``digits=17``) precision, and Fabel's loaders repackage that same export
without recomputation, so any real discrepancy is a bug, not accumulated
floating-point drift.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import fabel.datasets as ds

from .conftest import golden_cases

pytestmark = pytest.mark.parity

MODULE = "datasets"
RTOL = 1e-12
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_RELEASE_DIR = REPO_ROOT / "data_release"


def _case(name: str) -> dict[str, Any]:
    return next(c for c in golden_cases(MODULE) if c["name"] == name)


def _spot(arr: Any) -> tuple[np.ndarray, np.ndarray, float]:
    """Return ``(head3, tail3, sum)`` matching R's ``as.numeric()`` column-major flatten."""
    flat = np.asarray(arr, dtype=np.float64).flatten(order="F")
    return flat[:3], flat[-3:], float(np.nansum(flat))


def _assert_numeric(arr: Any, component: dict[str, Any], rtol: float = RTOL) -> None:
    head, tail, total = _spot(arr)
    exp_head = np.array([np.nan if v is None else v for v in component["head3"]])
    exp_tail = np.array([np.nan if v is None else v for v in component["tail3"]])
    np.testing.assert_allclose(head, exp_head, rtol=rtol, atol=1e-12, equal_nan=True)
    np.testing.assert_allclose(tail, exp_tail, rtol=rtol, atol=1e-12, equal_nan=True)
    np.testing.assert_allclose(total, component["sum"], rtol=rtol, atol=1e-12)


def _assert_strings(values: list[str], component: dict[str, Any]) -> None:
    assert len(values) == component["length"]
    assert values[:3] == component["head3"]
    assert values[-3:] == component["tail3"]


@pytest.fixture(autouse=True)
def _local_downloads(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Redirect every dataset "download" to the local ``data_release/`` staging files.

    No network: this is the fixture required by WORKFLOW.md so parity tests
    never hit the real GitHub release.
    """

    def _copy_local(url: str, dest: Path) -> None:
        name_ext = url.rsplit("/", 1)[-1]
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(DATA_RELEASE_DIR / name_ext, dest)

    monkeypatch.setattr(ds, "_download_file", _copy_local)
    monkeypatch.setenv("FABEL_DATA_DIR", str(tmp_path))


def test_canadian_weather() -> None:
    case = _case("dataset_CanadianWeather")
    comps = case["output"]["components"]
    cw = ds.load_canadian_weather()

    daily_av = np.stack([cw.temp, cw.precip, cw.log10precip], axis=2)
    _assert_numeric(daily_av, comps["dailyAv"])
    _assert_numeric(cw.coordinates, comps["coordinates"])
    _assert_numeric(cw.monthly_temp, comps["monthlyTemp"])
    _assert_numeric(cw.monthly_precip, comps["monthlyPrecip"])
    _assert_numeric(cw.geogindex, comps["geogindex"])
    _assert_strings(cw.stations, comps["place"])
    _assert_strings(cw.province, comps["province"])
    _assert_strings(cw.region, comps["region"])


def test_growth() -> None:
    case = _case("dataset_growth")
    comps = case["output"]["components"]
    g = ds.load_growth()

    _assert_numeric(g.hgtm, comps["hgtm"])
    _assert_numeric(g.hgtf, comps["hgtf"])
    _assert_numeric(g.age, comps["age"])


def test_gait() -> None:
    case = _case("dataset_gait")
    comps = case["output"]["components"]
    gait = ds.load_gait()

    _assert_numeric(gait.value, comps["value"])


def test_handwrit() -> None:
    case = _case("dataset_handwrit")
    comps = case["output"]["components"]
    hw = ds.load_handwriting()

    _assert_numeric(hw.value, comps["value"])


def test_handwrit_time() -> None:
    case = _case("dataset_handwritTime")
    comps = case["output"]["components"]
    hw = ds.load_handwriting()

    _assert_numeric(hw.t, comps["value"])


def test_pinch() -> None:
    case = _case("dataset_pinch")
    comps = case["output"]["components"]
    p = ds.load_pinch()

    _assert_numeric(p.pinch, comps["value"])


def test_melanoma() -> None:
    case = _case("dataset_melanoma")
    comps = case["output"]["components"]
    m = ds.load_melanoma()

    _assert_numeric(m.value, comps["value"])


def test_refinery() -> None:
    case = _case("dataset_refinery")
    comps = case["output"]["components"]
    r = ds.load_refinery()

    _assert_numeric(r.time, comps["Time"])
    _assert_numeric(r.reflux, comps["Reflux"])
    _assert_numeric(r.tray47, comps["Tray47"])


def test_seabird() -> None:
    case = _case("dataset_seabird")
    comps = case["output"]["components"]
    sb = ds.load_seabird()

    for species in (
        "BAGO",
        "BLSC",
        "COME",
        "COMU",
        "CORM",
        "HADU",
        "HOGR",
        "LOON",
        "MAMU",
        "OLDS",
        "PIGU",
        "RBME",
        "RNGR",
        "SUSC",
        "WWSC",
    ):
        _assert_numeric(sb.counts[species], comps[species])
    _assert_numeric(sb.year, comps["Year"])
    _assert_numeric(sb.site, comps["Site"])
    _assert_numeric(sb.transect, comps["Transect"])
    _assert_numeric(sb.temp, comps["Temp"])
    assert sb.observ_cond[:3] == comps["ObservCond"]["head3"]
    assert sb.observ_cond[-3:] == comps["ObservCond"]["tail3"]
    assert sb.bay[:3] == comps["Bay"]["head3"]
    assert sb.bay[-3:] == comps["Bay"]["tail3"]
    assert sb.observ_cond_factor3[:3] == comps["ObservCondFactor3"]["head3"]
    assert sb.observ_cond_factor3[-3:] == comps["ObservCondFactor3"]["tail3"]


def test_regina_precip() -> None:
    case = _case("dataset_ReginaPrecip")
    comps = case["output"]["components"]
    rp = ds.load_regina_precip()

    _assert_numeric(rp.value, comps["value"])


def test_catalog_scope() -> None:
    """Every non-"dateAccessories" catalog entry is covered by a Fabel loader."""
    case = _case("dataset_catalog")
    catalog = case["output"]["catalog"]

    covered = {
        "CanadianWeather",
        "MontrealTemp",
        "ReginaPrecip",
        "daily",
        "gait",
        "growth",
        "handwrit",
        "handwritTime",
        "infantGrowth",
        "lip",
        "lipmarks",
        "liptime",
        "melanoma",
        "nondurables",
        "pinch",
        "pinchraw",
        "pinchtime",
        "refinery",
        "seabird",
    }
    for entry in catalog:
        base = entry.split(" (")[0]
        if entry.endswith("(dateAccessories)"):
            continue
        assert base in covered, f"catalog entry {entry!r} has no Fabel loader"
