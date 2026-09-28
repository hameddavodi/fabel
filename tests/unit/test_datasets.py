"""Unit tests for :mod:`fdatools.datasets`: shape/dtype sanity, storage/cache mechanics."""

from __future__ import annotations

import inspect
import os
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any, TypeVar

import numpy as np
import pytest

import fdatools.datasets as ds

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_RELEASE_DIR = REPO_ROOT / "data_release"

_F = TypeVar("_F", bound=Callable[..., Any])


def offline(func: _F) -> _F:
    """Mark a test that needs no ``data_release/`` files, so it always runs."""
    func._fdatools_offline = True  # type: ignore[attr-defined]
    return func


@pytest.fixture(autouse=True)
def _local_downloads(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, request: pytest.FixtureRequest
) -> None:
    """Redirect every dataset "download" to the local ``data_release/`` staging files.

    Tests that need those files are skipped when ``data_release/`` is missing;
    tests marked :func:`offline` still run.
    """
    if not DATA_RELEASE_DIR.exists() and not getattr(request.function, "_fdatools_offline", False):
        pytest.skip("data_release/ fixtures missing; run tools/build_data_release.py first")

    def _copy_local(url: str, dest: Path) -> None:
        name_ext = url.rsplit("/", 1)[-1]
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(DATA_RELEASE_DIR / name_ext, dest)

    monkeypatch.setattr(ds, "_download_file", _copy_local)
    monkeypatch.setenv("FDATOOLS_DATA_DIR", str(tmp_path))


# --------------------------------------------------------------------------- #
# in-package loaders (no network, no monkeypatch needed)
# --------------------------------------------------------------------------- #


def test_load_growth_shapes() -> None:
    g = ds.load_growth()
    assert g.hgtm.shape == (31, 39)
    assert g.hgtf.shape == (31, 54)
    assert g.age.shape == (31,)
    assert g.hgtm.dtype == np.float64


def test_load_gait_shapes_and_views() -> None:
    gait = ds.load_gait()
    assert gait.value.shape == (20, 39, 2)
    assert gait.t.shape == (20,)
    assert len(gait.subjects) == 39
    assert gait.variables == ["Hip Angle", "Knee Angle"]
    np.testing.assert_array_equal(gait.hip_angle, gait.value[..., 0])
    np.testing.assert_array_equal(gait.knee_angle, gait.value[..., 1])


def test_load_pinch_shapes() -> None:
    p = ds.load_pinch()
    assert p.pinch.shape == (151, 20)
    assert p.pinchraw.shape == (151, 20)
    assert p.t.shape == (151,)


def test_in_package_loaders_never_touch_the_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def _boom(url: str, dest: Path) -> None:
        raise AssertionError(f"unexpected download attempt: {url}")

    monkeypatch.setattr(ds, "_download_file", _boom)
    ds.load_growth()
    ds.load_gait()
    ds.load_pinch()


# --------------------------------------------------------------------------- #
# downloaded + cached loaders
# --------------------------------------------------------------------------- #


def test_load_canadian_weather_shapes() -> None:
    cw = ds.load_canadian_weather()
    assert cw.temp.shape == (365, 35)
    assert cw.precip.shape == (365, 35)
    assert cw.log10precip.shape == (365, 35)
    assert cw.t.shape == (365,)
    assert len(cw.stations) == 35
    assert len(cw.province) == 35
    assert len(cw.region) == 35
    assert cw.coordinates.shape == (35, 2)
    assert cw.monthly_temp.shape == (12, 35)
    assert cw.monthly_precip.shape == (12, 35)
    assert cw.geogindex.shape == (35,)


def test_load_handwriting_shapes_and_views() -> None:
    hw = ds.load_handwriting()
    assert hw.value.shape == (1401, 20, 2)
    assert hw.t.shape == (1401,)
    assert len(hw.subjects) == 20
    np.testing.assert_array_equal(hw.x, hw.value[..., 0])
    np.testing.assert_array_equal(hw.y, hw.value[..., 1])


def test_load_melanoma_shape() -> None:
    m = ds.load_melanoma()
    assert m.value.shape == (37, 3)
    assert m.columns == ["index", "year", "incidence"]


def test_load_refinery_shapes() -> None:
    r = ds.load_refinery()
    assert r.time.shape == (194,)
    assert r.reflux.shape == (194,)
    assert r.tray47.shape == (194,)


def test_load_seabird_shapes() -> None:
    sb = ds.load_seabird()
    assert sb.year.shape == (3793,)
    assert len(sb.counts) == 15
    assert all(arr.shape == (3793,) for arr in sb.counts.values())
    assert len(sb.observ_cond) == 3793
    assert len(sb.bay) == 3793


def test_load_regina_precip_shape() -> None:
    rp = ds.load_regina_precip()
    assert rp.value.shape == (1006,)


def test_load_montreal_temp_shape() -> None:
    mt = ds.load_montreal_temp()
    assert mt.value.shape == (34, 365)
    assert len(mt.years) == 34


def test_load_daily_shapes() -> None:
    d = ds.load_daily()
    assert d.tempav.shape == (365, 35)
    assert d.precav.shape == (365, 35)
    assert len(d.place) == 35


def test_load_infant_growth_shapes() -> None:
    ig = ds.load_infant_growth()
    assert ig.day.shape == (40,)
    assert ig.tibia_length.shape == (40,)
    assert ig.sd_length.shape == (40,)


def test_load_nondurables_shape() -> None:
    nd = ds.load_nondurables()
    assert nd.value.shape == (1377,)
    assert nd.start == "1919-01"
    assert nd.frequency == 12


def test_load_lip_shapes() -> None:
    lip = ds.load_lip()
    assert lip.value.shape == (51, 20)
    assert lip.t.shape == (51,)
    assert lip.left_elbow.shape == (20,)
    assert lip.right_elbow.shape == (20,)


# --------------------------------------------------------------------------- #
# cache / checksum mechanics
# --------------------------------------------------------------------------- #


def test_cache_dir_honours_env_var(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("FDATOOLS_DATA_DIR", str(tmp_path / "custom"))
    assert ds._cache_dir() == tmp_path / "custom"


def test_cache_dir_defaults_to_home_cache(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("FDATOOLS_DATA_DIR", raising=False)
    assert ds._cache_dir() == Path.home() / ".cache" / "fdatools"


def test_second_load_reuses_cache_without_downloading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FDATOOLS_DATA_DIR", str(tmp_path))
    calls = []

    def _copy_local(url: str, dest: Path) -> None:
        calls.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(DATA_RELEASE_DIR / url.rsplit("/", 1)[-1], dest)

    monkeypatch.setattr(ds, "_download_file", _copy_local)
    ds.load_regina_precip()
    assert len(calls) == 2  # npz + json
    ds.load_regina_precip()
    assert len(calls) == 2  # cached, no new download


def test_checksum_mismatch_raises(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("FDATOOLS_DATA_DIR", str(tmp_path))

    def _corrupt(url: str, dest: Path) -> None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"not the expected bytes")

    monkeypatch.setattr(ds, "_download_file", _corrupt)
    with pytest.raises(ValueError, match="checksum mismatch"):
        ds.load_melanoma()


def test_download_file_rejects_non_https(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Undo the autouse _local_downloads monkeypatch: this test exercises the
    # real _download_file, not the fixture-copying stand-in.
    monkeypatch.undo()
    with pytest.raises(ValueError, match="https"):
        ds._download_file("http://example.com/x.npz", tmp_path / "x.npz")


@pytest.mark.network
@pytest.mark.skipif(
    os.environ.get("FDATOOLS_RUN_NETWORK_TESTS") != "1",
    reason="hits the real network; set FDATOOLS_RUN_NETWORK_TESTS=1 to run",
)
def test_real_download_from_github_release(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Hits the real GitHub release URL. Skipped by default.

    Run explicitly with ``FDATOOLS_RUN_NETWORK_TESTS=1 pytest -m network`` once
    ``data-v1`` has been published (see ``docs/dev/data-release.md``).
    """
    monkeypatch.undo()  # use the real _download_file, not the local-fixture stand-in
    monkeypatch.setenv("FDATOOLS_DATA_DIR", str(tmp_path))
    rp = ds.load_regina_precip()
    assert rp.value.shape == (1006,)


# --------------------------------------------------------------------------- #
# docstrings: every value/time field states its physical unit (R help pages)
# --------------------------------------------------------------------------- #

# (class, attribute as written in the Attributes section, text naming its unit)
_DOCUMENTED_UNITS = [
    (ds.CanadianWeather, "temp, precip", "deg C"),
    (ds.CanadianWeather, "temp, precip", "mm"),
    (ds.CanadianWeather, "t", "Day of year"),
    (ds.CanadianWeather, "monthly_temp, monthly_precip", "deg C"),
    (ds.CanadianWeather, "monthly_temp, monthly_precip", "mm"),
    (ds.Growth, "hgtm", "cm"),
    (ds.Growth, "hgtf", "cm"),
    (ds.Growth, "age", "years"),
    (ds.Gait, "value", "degrees"),
    (ds.Gait, "t", "proportion of the cycle"),
    (ds.Handwriting, "value", "metres"),
    (ds.Handwriting, "t", "milliseconds"),
    (ds.Pinch, "pinch", "newtons"),
    (ds.Pinch, "pinchraw", "newtons"),
    (ds.Pinch, "t", "seconds"),
    (ds.Melanoma, "value", "per 100,000"),
    (ds.Melanoma, "value", "year"),
    (ds.Refinery, "time", "no unit"),
    (ds.Refinery, "reflux", "no unit"),
    (ds.Refinery, "tray47", "no unit"),
    (ds.Seabird, "counts", "sightings"),
    (ds.Seabird, "year, site, transect, temp", "no unit"),
    (ds.ReginaPrecip, "value", "mm"),
    (ds.MontrealTemp, "value", "degrees Celsius"),
    (ds.Daily, "tempav, precav", "deg C"),
    (ds.Daily, "tempav, precav", "mm"),
    (ds.InfantGrowth, "day", "days"),
    (ds.InfantGrowth, "tibia_length", "mm"),
    (ds.InfantGrowth, "sd_length", "mm"),
    (ds.Nondurables, "value", "no unit"),
    (ds.Lip, "value", "mm"),
    (ds.Lip, "t", "seconds"),
]


def _attribute_text(cls: type, name: str) -> str:
    """Return the description block of attribute ``name`` in ``cls``'s docstring."""
    lines = inspect.cleandoc(cls.__doc__ or "").splitlines()
    start = lines.index("Attributes") + 2
    block: list[str] = []
    inside = False
    for line in lines[start:]:
        if line and not line.startswith(" "):
            if inside:
                break
            inside = line.split(" : ")[0] == name
        elif inside:
            block.append(line.strip())
    assert block, f"{cls.__name__} documents no attribute {name!r}"
    return " ".join(block)


@offline
@pytest.mark.parametrize(
    ("cls", "attribute", "unit"),
    _DOCUMENTED_UNITS,
    ids=[f"{c.__name__}.{a}.{u}" for c, a, u in _DOCUMENTED_UNITS],
)
def test_dataset_docstring_states_unit(cls: type, attribute: str, unit: str) -> None:
    assert unit in _attribute_text(cls, attribute)


@offline
def test_every_dataset_class_is_audited_for_units() -> None:
    classes = {
        obj
        for obj in vars(ds).values()
        if inspect.isclass(obj)
        and obj.__module__ == ds.__name__
        and "Attributes" in (obj.__doc__ or "")
    }
    assert classes == {cls for cls, _, _ in _DOCUMENTED_UNITS}
