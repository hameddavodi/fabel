"""Dataset access for parity tests.

Golden cases built on a real R ``fda`` dataset record only the design (argvals,
basis, lambda), never the response matrix -- the data would dominate the golden
file.  The arrays are read from the ``data_export/`` dump produced alongside the
golden files; when that dump is absent the parity test is skipped rather than
failed, so a checkout without it still runs the synthetic cases.
"""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[2] / "data_export"

#: Third axis of ``CanadianWeather$dailyAv``.
WEATHER_VARIABLES = ("Temperature.C", "Precipitation.mm", "log10precip")


@cache
def _load(name: str) -> Any:
    path = DATA_DIR / f"{name}.json"
    if not path.exists():
        pytest.skip(f"dataset dump {path} is not available")
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def weather_daily(variable: str = "Temperature.C") -> np.ndarray:
    """Return ``CanadianWeather$dailyAv[, , variable]`` as a ``(365, 35)`` array."""
    raw = np.asarray(_load("CanadianWeather")["dailyAv"], dtype=float)
    return np.ascontiguousarray(raw[:, :, WEATHER_VARIABLES.index(variable)])


def weather_region() -> np.ndarray:
    """Return the 35 station regions."""
    return np.asarray(_load("CanadianWeather")["region"], dtype=object)


def weather_place() -> np.ndarray:
    """Return the 35 station names."""
    return np.asarray(_load("CanadianWeather")["place"], dtype=object)


def growth(field: str) -> np.ndarray:
    """Return one of ``growth$hgtm``, ``growth$hgtf`` or ``growth$age``."""
    return np.asarray(_load("growth")[field], dtype=float)


def handwriting() -> tuple[np.ndarray, np.ndarray]:
    """Return ``(handwritTime, handwrit)`` shaped ``(1401,)`` and ``(1401, 20, 2)``."""
    return (
        np.asarray(_load("handwritTime"), dtype=float),
        np.asarray(_load("handwrit"), dtype=float),
    )
