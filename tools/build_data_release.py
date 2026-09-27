"""Build ``data_release/<name>.npz`` (+ ``.json`` sidecar) for every fabel dataset.

Reads the full-precision R dumps in ``data_export/<name>.json`` (produced by
``Rscript tools/export_datasets.R``, gitignored, not committed) and writes, for
each dataset:

- ``data_release/<name>.npz`` -- every numeric array, cast to ``float64``,
  compressed with :func:`numpy.savez_compressed`.
- ``data_release/<name>.json`` -- string/label metadata (station names,
  subject ids, factor levels, column order) that does not belong in an npz.

``data_release/`` is gitignored: it is regenerated from ``data_export/`` (in
turn regenerated from R) and is the staging area for the human release step
(see ``docs/dev/data-release.md``). This script also prints the SHA-256 of
every file it writes -- copy that table into the ``_CHECKSUMS`` constant in
``src/fabel/datasets.py``.

Usage::

    .venv/bin/python tools/build_data_release.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_EXPORT_DIR = REPO_ROOT / "data_export"
DATA_RELEASE_DIR = REPO_ROOT / "data_release"


def _load_export(name: str) -> Any:
    """Return the parsed ``data_export/<name>.json`` document."""
    path = DATA_EXPORT_DIR / f"{name}.json"
    with path.open(encoding="utf-8") as fh:
        result: Any = json.load(fh)
        return result


def _arr(data: Any) -> np.ndarray:
    """Convert a (possibly nested) JSON list to a ``float64`` array."""
    return np.asarray(data, dtype=np.float64)


# Each builder returns (arrays, meta): arrays is float64-only (-> .npz),
# meta is JSON-safe string/int metadata (-> .json sidecar).


def _build_canadian_weather() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    d = _load_export("CanadianWeather")
    daily_av = _arr(d["dailyAv"])  # (365, 35, 3): Temperature.C, Precipitation.mm, log10precip
    arrays = {
        "temp": daily_av[:, :, 0],
        "precip": daily_av[:, :, 1],
        "log10precip": daily_av[:, :, 2],
        "coordinates": _arr(d["coordinates"]),
        "monthly_temp": _arr(d["monthlyTemp"]),
        "monthly_precip": _arr(d["monthlyPrecip"]),
        "geogindex": _arr(d["geogindex"]),
        "t": np.arange(1, 366, dtype=np.float64),
    }
    meta = {
        "stations": d["place"],
        "province": d["province"],
        "region": d["region"],
    }
    return arrays, meta


def _build_growth() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    d = _load_export("growth")
    arrays = {"hgtm": _arr(d["hgtm"]), "hgtf": _arr(d["hgtf"]), "age": _arr(d["age"])}
    return arrays, {}


def _build_gait() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    value = _arr(_load_export("gait"))  # (20, 39, 2): Hip Angle, Knee Angle
    t = np.arange(0.025, 1.0, 0.05, dtype=np.float64)[: value.shape[0]]
    arrays = {"value": value, "t": t}
    meta = {
        "subjects": [f"boy{i}" for i in range(1, value.shape[1] + 1)],
        "variables": ["Hip Angle", "Knee Angle"],
    }
    return arrays, meta


def _build_handwrit() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    value = _arr(_load_export("handwrit"))  # (1401, 20, 2): X, Y
    t = _arr(_load_export("handwritTime"))
    arrays = {"value": value, "t": t}
    meta = {
        "subjects": [f"rep{i:02d}" for i in range(1, value.shape[1] + 1)],
        "variables": ["X", "Y"],
    }
    return arrays, meta


def _build_pinch() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    arrays = {
        "pinch": _arr(_load_export("pinch")),
        "pinchraw": _arr(_load_export("pinchraw")),
        "t": _arr(_load_export("pinchtime")),
    }
    return arrays, {}


def _build_melanoma() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    arrays = {"value": _arr(_load_export("melanoma"))}  # (37, 3): index, year, incidence
    meta = {"columns": ["index", "year", "incidence"]}
    return arrays, meta


def _build_refinery() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    d = _load_export("refinery")  # list of row dicts (jsonlite data.frame-as-columns... see below)
    # export_datasets.R serializes a data.frame as a JSON array of row objects
    # (jsonlite's default for a list of equal-length atomic columns is
    # column-major, but a plain `toJSON(df)` call -- no dataframe="columns"
    # override -- emits one object per row); pivot back to columns here.
    columns = {key: [row[key] for row in d] for key in d[0]}
    arrays = {name: _arr(values) for name, values in columns.items()}
    meta = {"columns": list(columns.keys())}
    return arrays, meta


def _build_seabird() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    d = _load_export("seabird")
    columns = {key: [row[key] for row in d] for key in d[0]}
    factor_cols = ["ObservCond", "Bay", "ObservCondFactor3"]
    arrays = {
        name: _arr([np.nan if v is None else v for v in values])
        for name, values in columns.items()
        if name not in factor_cols
    }
    meta = {
        "columns": list(columns.keys()),
        **{name: columns[name] for name in factor_cols},
    }
    return arrays, meta


def _build_regina_precip() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    return {"value": _arr(_load_export("ReginaPrecip"))}, {}


def _build_montreal_temp() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    value = _arr(_load_export("MontrealTemp"))  # (34, 365)
    meta = {
        "years": [str(y) for y in range(1961, 1961 + value.shape[0])],
    }
    return {"value": value}, meta


def _build_daily() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    d = _load_export("daily")
    arrays = {"tempav": _arr(d["tempav"]), "precav": _arr(d["precav"])}
    meta = {"place": d["place"]}
    return arrays, meta


def _build_infant_growth() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    value = _arr(_load_export("infantGrowth"))  # (40, 3): day, tibiaLength, sd.length
    arrays = {"day": value[:, 0], "tibia_length": value[:, 1], "sd_length": value[:, 2]}
    return arrays, {}


def _build_nondurables() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    return {"value": _arr(_load_export("nondurables"))}, {"start": "1919-01", "frequency": 12}


def _build_lip() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    value = _arr(_load_export("lip"))  # (51, 20)
    t = _arr(_load_export("liptime"))
    marks = _load_export("lipmarks")  # list of {"leftElbow":..., "rightElbow":...}
    arrays = {
        "value": value,
        "t": t,
        "left_elbow": _arr([row["leftElbow"] for row in marks]),
        "right_elbow": _arr([row["rightElbow"] for row in marks]),
    }
    return arrays, {}


_BUILDERS = {
    "CanadianWeather": _build_canadian_weather,
    "growth": _build_growth,
    "gait": _build_gait,
    "handwrit": _build_handwrit,
    "pinch": _build_pinch,
    "melanoma": _build_melanoma,
    "refinery": _build_refinery,
    "seabird": _build_seabird,
    "ReginaPrecip": _build_regina_precip,
    "MontrealTemp": _build_montreal_temp,
    "daily": _build_daily,
    "infantGrowth": _build_infant_growth,
    "nondurables": _build_nondurables,
    "lip": _build_lip,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    """Write every release asset and its SHA-256 checksum into ``data_release/``."""
    DATA_RELEASE_DIR.mkdir(parents=True, exist_ok=True)
    checksums: dict[str, dict[str, str]] = {}
    for name, builder in _BUILDERS.items():
        arrays, meta = builder()
        npz_path = DATA_RELEASE_DIR / f"{name}.npz"
        json_path = DATA_RELEASE_DIR / f"{name}.json"
        # allow_pickle is numpy's default, spelled out so the stub cannot match a
        # ``**arrays`` entry against that bool keyword.
        np.savez_compressed(npz_path, allow_pickle=True, **arrays)
        with json_path.open("w", encoding="utf-8") as fh:
            json.dump(meta, fh)
        checksums[name] = {"npz": _sha256(npz_path), "json": _sha256(json_path)}
        print(f"wrote {npz_path} ({npz_path.stat().st_size} bytes)")
        print(f"wrote {json_path} ({json_path.stat().st_size} bytes)")

    print("\n# Copy into src/fabel/datasets.py::_CHECKSUMS\n")
    print("_CHECKSUMS: dict[str, dict[str, str]] = {")
    for name, sums in checksums.items():
        print(f'    "{name}": {{"npz": "{sums["npz"]}", "json": "{sums["json"]}"}},')
    print("}")


if __name__ == "__main__":
    main()
