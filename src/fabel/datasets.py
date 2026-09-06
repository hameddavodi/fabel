"""Loaders for the real datasets shipped with R's ``fda`` package (v6.3.0).

Every ``load_*()`` function returns a frozen dataclass of raw NumPy arrays plus
the string/label metadata R attaches to the object (station names, subject
ids, factor levels, ...), reproduced exactly as ``fda`` stores them. None of
these datasets is already a smooth functional object in R (``fda`` stores them
as plain matrices/arrays/data frames), so none of the dataclasses below carries
an ``.fd`` property -- fit one with :func:`fabel.smooth` using the raw arrays
and the accompanying ``t`` grid.

Storage: ``growth``, ``gait`` and ``pinch`` ship inside the wheel
(``src/fabel/_data/``, read via :mod:`importlib.resources`). Every other
dataset is downloaded on first use from the project's GitHub release
``data-v1``, SHA-256 verified, and cached under ``$FABEL_DATA_DIR`` (default
``~/.cache/fabel``). See ``docs/dev/data-release.md`` for how the release
archives are built and published.

``fda`` 6.3.0's dataset catalog has no ``CSTR`` dataset and no separate
"Chinese script" dataset (see ``tests/golden/datasets.json``'s
``dataset_catalog`` case) -- ``load_cstr``/``load_chinese_script`` from
SPEC.md §8 are not implemented; see ``PROGRESS.md`` for the record. The
catalog's "dateAccessories" group (``day.5``, ``monthLetters``, ...) are
plotting-axis constants, not datasets, and are likewise out of scope.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import urllib.request
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "CanadianWeather",
    "Daily",
    "Gait",
    "Growth",
    "Handwriting",
    "InfantGrowth",
    "Lip",
    "Melanoma",
    "MontrealTemp",
    "Nondurables",
    "Pinch",
    "Refinery",
    "ReginaPrecip",
    "Seabird",
    "load_canadian_weather",
    "load_daily",
    "load_gait",
    "load_growth",
    "load_handwriting",
    "load_infant_growth",
    "load_lip",
    "load_melanoma",
    "load_montreal_temp",
    "load_nondurables",
    "load_pinch",
    "load_refinery",
    "load_regina_precip",
    "load_seabird",
]

#: fda 6.3.0 datasets shipped inside the wheel (small; no network needed).
_IN_PACKAGE = frozenset({"growth", "gait", "pinch"})

#: GitHub release asset URL template for every other dataset.
_RELEASE_URL = "https://github.com/AISMAsrl/fabel/releases/download/data-v1/{name}.{ext}"

#: SHA-256 of every ``data_release/<name>.{npz,json}`` file (see
#: ``tools/build_data_release.py``, which prints this table). Verified after
#: every download; a mismatch raises instead of silently loading bad data.
_CHECKSUMS: dict[str, dict[str, str]] = {
    "CanadianWeather": {
        "npz": "2a24cbeaccb57ab9ec956a7fc16827fac39e48add95b9011e72006e07bb004a3",
        "json": "550d585700e9a7d4ff28758c8b7240f0e878ea059ea37dc41dfe807eebe0006b",
    },
    "growth": {
        "npz": "49c39d9a818a32f6ac66334b5b58304e1e1da3e5667e6f03b1fbd6fe7612d5b8",
        "json": "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a",
    },
    "gait": {
        "npz": "9012b748ae9eafb923cd81ca2de387cda14792701ef82c625a507505ba2ba882",
        "json": "ffb9f8972ceaa441702c932067e7a43d46c073e65e7952474a8bbc5c0b38f284",
    },
    "handwrit": {
        "npz": "2e46b68ae70fc1a0c8b1a25ade9fee12289cd69bc58e22a940bfa57f3fbafc66",
        "json": "86034e3ef6f913158b0599c3870a8ac472b09db5ff733dfd6e3a2bd86a57cb75",
    },
    "pinch": {
        "npz": "4ac050d8d4a3d636776f24d783ddeb0e42ba9496ef7fa6096fb8da6bc3d9d9ac",
        "json": "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a",
    },
    "melanoma": {
        "npz": "61fa23a586ac9de3e1a4e595f6d6480178fffcdf0210bde0947ec0c50081d1bf",
        "json": "9644293a416d940f17720027b515381e8583c48a7e43edf0c9ad1e250ba24af2",
    },
    "refinery": {
        "npz": "5fa33ce9c0ee13e3c400a6788a946f0a4e49b353524a8fd7fee25492c268e760",
        "json": "2a60ea9c9462882c9cee20f639fa337d03e8145e18b1228abd76b74c4ba5b0d2",
    },
    "seabird": {
        "npz": "a2768a6516631f9beeb87fc2424a0f2ebf29bb6715840309d23a9935df505388",
        "json": "28ae841d5c3f307ec233ee8dc57ebf1a3ec94c2fb88dafb9a5794dc1784821c1",
    },
    "ReginaPrecip": {
        "npz": "97ff9a0d9eaa2d2651a1eacb1f9a886c23c0dd7efc57a0d06a6a4f42fe35c3a8",
        "json": "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a",
    },
    "MontrealTemp": {
        "npz": "dbb00516f5ebc400fc209e57cb2749d3bef502de2eee2426c858ced029a1c482",
        "json": "e236071dd3d140706911ab3d21b405445a6a09259bb72248e3f2dd3a7c509d39",
    },
    "daily": {
        "npz": "e2354ea84784ffd288e0918d21f6039ab0c02d3728c495c0ca2e088e988dc9df",
        "json": "45c00eb711439f350ae5ca8ce99fe5b35abcbc26dee4e4621f6a56fd9cf8ad83",
    },
    "infantGrowth": {
        "npz": "b6ba10866603c75910c4bbee09da7ae3c1b46b05fc1fdbc06dc1c8b60a7d5be7",
        "json": "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a",
    },
    "nondurables": {
        "npz": "ad269753677f45a152a9449d808e504273cc258c2afb9a1ff291439876dc21ba",
        "json": "ed6a435416dda383714ba7cf28e2b9cff757c7bf8a5b61b42d8986f294cf7ba7",
    },
    "lip": {
        "npz": "a89784c69a6bcbd15025d0694290e0822dd545f44855a4e42b74c5a56da4a9ee",
        "json": "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a",
    },
}


def _cache_dir() -> Path:
    """Return the local cache directory, honouring ``$FABEL_DATA_DIR``."""
    env = os.environ.get("FABEL_DATA_DIR")
    return Path(env) if env else Path.home() / ".cache" / "fabel"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_file(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest`` (atomic: write to a temp file, then rename).

    A separate, monkeypatchable function so tests can redirect it to a local
    fixture instead of the network (see ``tests/unit/test_datasets.py``).
    """
    if not url.startswith("https://"):
        raise ValueError(f"refusing to download from a non-https URL: {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response:  # noqa: S310 -- scheme checked above
        data = response.read()
    tmp = dest.with_name(dest.name + ".part")
    tmp.write_bytes(data)
    tmp.replace(dest)


def _ensure_cached(name: str) -> tuple[Path, Path]:
    """Return local ``(npz_path, json_path)`` for ``name``, downloading+verifying if needed."""
    cache = _cache_dir()
    checksums = _CHECKSUMS[name]
    npz_path = cache / f"{name}.npz"
    json_path = cache / f"{name}.json"
    for path, ext in ((npz_path, "npz"), (json_path, "json")):
        if path.exists() and _sha256(path) == checksums[ext]:
            continue
        _download_file(_RELEASE_URL.format(name=name, ext=ext), path)
        actual = _sha256(path)
        if actual != checksums[ext]:
            path.unlink(missing_ok=True)
            raise ValueError(
                f"checksum mismatch downloading {name}.{ext}: "
                f"expected {checksums[ext]}, got {actual}"
            )
    return npz_path, json_path


def _load_in_package(name: str) -> tuple[dict[str, NDArray[np.float64]], dict[str, Any]]:
    pkg = resources.files("fabel._data")
    npz_bytes = pkg.joinpath(f"{name}.npz").read_bytes()
    meta_text = pkg.joinpath(f"{name}.json").read_text(encoding="utf-8")
    with np.load(io.BytesIO(npz_bytes)) as npz:
        arrays = {key: npz[key] for key in npz.files}
    meta: dict[str, Any] = json.loads(meta_text)
    return arrays, meta


def _load_downloaded(name: str) -> tuple[dict[str, NDArray[np.float64]], dict[str, Any]]:
    npz_path, json_path = _ensure_cached(name)
    with np.load(npz_path) as npz:
        arrays = {key: npz[key] for key in npz.files}
    meta: dict[str, Any] = json.loads(json_path.read_text(encoding="utf-8"))
    return arrays, meta


def _load(name: str) -> tuple[dict[str, NDArray[np.float64]], dict[str, Any]]:
    """Return ``(arrays, meta)`` for dataset ``name``, in-package or cached-download."""
    if name in _IN_PACKAGE:
        return _load_in_package(name)
    return _load_downloaded(name)


@dataclass(frozen=True, eq=False)
class CanadianWeather:
    """Daily temperature/precipitation at 35 Canadian weather stations, 1960-1994.

    Reproduces R ``fda::CanadianWeather`` (dailyAv's three layers split into
    separate fields) plus ``CanadianWeather$monthlyTemp``/``monthlyPrecip``.

    Attributes
    ----------
    temp, precip : ndarray, shape (365, 35)
        Average daily temperature (deg C) and precipitation (mm).
    log10precip : ndarray, shape (365, 35)
        ``log10(precip)`` with the 27 zero readings replaced by 0.05 mm first.
    t : ndarray, shape (365,)
        Day of year, 1-365.
    stations : list of str
        35 station names (R: ``CanadianWeather$place``).
    province, region : list of str
        Province and one of the 4 climate zones (Atlantic/Pacific/Continental/Arctic)
        for each station.
    coordinates : ndarray, shape (35, 2)
        North latitude, west longitude.
    monthly_temp, monthly_precip : ndarray, shape (12, 35)
        Monthly averages.
    geogindex : ndarray, shape (35,)
        Station order east-to-west-then-north (R: ``CanadianWeather$geogindex``).
    """

    temp: NDArray[np.float64]
    precip: NDArray[np.float64]
    log10precip: NDArray[np.float64]
    t: NDArray[np.float64]
    stations: list[str]
    province: list[str]
    region: list[str]
    coordinates: NDArray[np.float64]
    monthly_temp: NDArray[np.float64]
    monthly_precip: NDArray[np.float64]
    geogindex: NDArray[np.float64]


def load_canadian_weather() -> CanadianWeather:
    """Load :class:`CanadianWeather` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> cw = fb.datasets.load_canadian_weather()  # doctest: +SKIP
    >>> cw.temp.shape  # doctest: +SKIP
    (365, 35)
    """
    arrays, meta = _load("CanadianWeather")
    return CanadianWeather(
        temp=arrays["temp"],
        precip=arrays["precip"],
        log10precip=arrays["log10precip"],
        t=arrays["t"],
        stations=list(meta["stations"]),
        province=list(meta["province"]),
        region=list(meta["region"]),
        coordinates=arrays["coordinates"],
        monthly_temp=arrays["monthly_temp"],
        monthly_precip=arrays["monthly_precip"],
        geogindex=arrays["geogindex"],
    )


@dataclass(frozen=True, eq=False)
class Growth:
    """Berkeley Growth Study: heights of 39 boys and 54 girls, ages 1-18.

    Reproduces R ``fda::growth`` (``hgtm``, ``hgtf``, ``age``). Ages are not
    equally spaced.

    Attributes
    ----------
    hgtm : ndarray, shape (31, 39)
        Boys' heights in cm at each of 31 ages.
    hgtf : ndarray, shape (31, 54)
        Girls' heights in cm at each of 31 ages.
    age : ndarray, shape (31,)
        Ages in years at which heights were measured.
    """

    hgtm: NDArray[np.float64]
    hgtf: NDArray[np.float64]
    age: NDArray[np.float64]


def load_growth() -> Growth:
    """Load :class:`Growth` (ships in-package, no network).

    Examples
    --------
    >>> import fabel as fb
    >>> g = fb.datasets.load_growth()
    >>> g.hgtm.shape, g.hgtf.shape, g.age.shape
    ((31, 39), (31, 54), (31,))
    """
    arrays, _meta = _load("growth")
    return Growth(hgtm=arrays["hgtm"], hgtf=arrays["hgtf"], age=arrays["age"])


@dataclass(frozen=True, eq=False)
class Gait:
    """Hip and knee angle (degrees) through a 20-point gait cycle, 39 boys.

    Reproduces R ``fda::gait``, an array of dim ``(20, 39, 2)``.

    Attributes
    ----------
    value : ndarray, shape (20, 39, 2)
        ``value[..., 0]`` is hip angle, ``value[..., 1]`` is knee angle
        (see :attr:`hip_angle`/:attr:`knee_angle`).
    t : ndarray, shape (20,)
        Standardized gait time, proportion of the cycle (0.025 to 0.975).
    subjects : list of str
        39 subject ids, "boy1".."boy39".
    variables : list of str
        ``["Hip Angle", "Knee Angle"]``, the meaning of ``value``'s last axis.
    """

    value: NDArray[np.float64]
    t: NDArray[np.float64]
    subjects: list[str]
    variables: list[str]

    @property
    def hip_angle(self) -> NDArray[np.float64]:
        """Return ``value[..., 0]``, shape (20, 39)."""
        return self.value[..., 0]

    @property
    def knee_angle(self) -> NDArray[np.float64]:
        """Return ``value[..., 1]``, shape (20, 39)."""
        return self.value[..., 1]


def load_gait() -> Gait:
    """Load :class:`Gait` (ships in-package, no network).

    Examples
    --------
    >>> import fabel as fb
    >>> gait = fb.datasets.load_gait()
    >>> gait.value.shape
    (20, 39, 2)
    """
    arrays, meta = _load("gait")
    return Gait(
        value=arrays["value"],
        t=arrays["t"],
        subjects=list(meta["subjects"]),
        variables=list(meta["variables"]),
    )


@dataclass(frozen=True, eq=False)
class Handwriting:
    """20 cursive samples of 1401 (x, y) coordinates writing "fda".

    Reproduces R ``fda::handwrit`` + ``fda::handwritTime``. The subject was
    Jim Ramsay; each replicate has been resampled to a common 2300 ms length
    and registered so key features align across replicates.

    Attributes
    ----------
    value : ndarray, shape (1401, 20, 2)
        ``value[..., 0]`` is X, ``value[..., 1]`` is Y (see :attr:`x`/:attr:`y`).
    t : ndarray, shape (1401,)
        Sampling times in milliseconds, 0 to 2300.
    subjects : list of str
        20 replicate ids, "rep01".."rep20".
    variables : list of str
        ``["X", "Y"]``, the meaning of ``value``'s last axis.
    """

    value: NDArray[np.float64]
    t: NDArray[np.float64]
    subjects: list[str]
    variables: list[str]

    @property
    def x(self) -> NDArray[np.float64]:
        """Return ``value[..., 0]``, shape (1401, 20)."""
        return self.value[..., 0]

    @property
    def y(self) -> NDArray[np.float64]:
        """Return ``value[..., 1]``, shape (1401, 20)."""
        return self.value[..., 1]


def load_handwriting() -> Handwriting:
    """Load :class:`Handwriting` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> hw = fb.datasets.load_handwriting()  # doctest: +SKIP
    >>> hw.value.shape  # doctest: +SKIP
    (1401, 20, 2)
    """
    arrays, meta = _load("handwrit")
    return Handwriting(
        value=arrays["value"],
        t=arrays["t"],
        subjects=list(meta["subjects"]),
        variables=list(meta["variables"]),
    )


@dataclass(frozen=True, eq=False)
class Pinch:
    """151 measurements of pinch force across 20 replications.

    Reproduces R ``fda::pinch`` + ``pinchraw`` + ``pinchtime``. The original
    recording had 300 observations every 2 ms; ``pinchraw`` keeps the first
    151, while ``pinch`` instead selects the 151 observations so that every
    curve's maximum falls at 0.076 s.

    Attributes
    ----------
    pinch : ndarray, shape (151, 20)
        Pinch force (N), aligned so each curve peaks at 0.076 s.
    pinchraw : ndarray, shape (151, 20)
        Pinch force (N), the first 151 of the original 300 samples, unaligned.
    t : ndarray, shape (151,)
        Time in seconds from the start, every 2 ms (R: ``pinchtime``).
    """

    pinch: NDArray[np.float64]
    pinchraw: NDArray[np.float64]
    t: NDArray[np.float64]


def load_pinch() -> Pinch:
    """Load :class:`Pinch` (ships in-package, no network).

    Examples
    --------
    >>> import fabel as fb
    >>> p = fb.datasets.load_pinch()
    >>> p.pinch.shape, p.pinchraw.shape, p.t.shape
    ((151, 20), (151, 20), (151,))
    """
    arrays, _meta = _load("pinch")
    return Pinch(pinch=arrays["pinch"], pinchraw=arrays["pinchraw"], t=arrays["t"])


@dataclass(frozen=True, eq=False)
class Melanoma:
    """Age-adjusted melanoma incidence per 100,000 in Connecticut, 1936-1972.

    Reproduces R ``fda::melanoma``. The R object is a bare 37x3 matrix; its
    first column is a plain 1..37 row index (the documented format lists
    only ``year``/``incidence``, but the shipped object carries a leading
    index column -- reproduced as-is, see ``columns`` for the mapping).

    Attributes
    ----------
    value : ndarray, shape (37, 3)
        Columns as listed in ``columns``.
    columns : list of str
        ``["index", "year", "incidence"]``.
    """

    value: NDArray[np.float64]
    columns: list[str]


def load_melanoma() -> Melanoma:
    """Load :class:`Melanoma` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> m = fb.datasets.load_melanoma()  # doctest: +SKIP
    >>> m.value.shape  # doctest: +SKIP
    (37, 3)
    """
    arrays, meta = _load("melanoma")
    return Melanoma(value=arrays["value"], columns=list(meta["columns"]))


@dataclass(frozen=True, eq=False)
class Refinery:
    """194 observations of reflux and "tray 47 level" in an oil refinery column.

    Reproduces R ``fda::refinery`` (a data frame with 3 numeric columns).

    Attributes
    ----------
    time : ndarray, shape (194,)
        Observation time, 0-193.
    reflux : ndarray, shape (194,)
        Reflux flow, centered on the mean of the first 60 observations.
    tray47 : ndarray, shape (194,)
        Tray 47 level, centered on the mean of the first 60 observations.
    """

    time: NDArray[np.float64]
    reflux: NDArray[np.float64]
    tray47: NDArray[np.float64]


def load_refinery() -> Refinery:
    """Load :class:`Refinery` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> r = fb.datasets.load_refinery()  # doctest: +SKIP
    >>> r.time.shape  # doctest: +SKIP
    (194,)
    """
    arrays, _meta = _load("refinery")
    return Refinery(time=arrays["Time"], reflux=arrays["Reflux"], tray47=arrays["Tray47"])


@dataclass(frozen=True, eq=False)
class Seabird:
    """Seabird sighting counts by transect, 4 Kodiak Island bays, 1986-2005.

    Reproduces R ``fda::seabird`` (a data frame with 22 columns): 15 species
    count columns, 3 integer id columns, one numeric column, and 3 factor
    columns.

    Attributes
    ----------
    counts : dict of str to ndarray, each shape (3793,)
        One entry per species code (``BAGO``, ``BLSC``, ``COME``, ``COMU``,
        ``CORM``, ``HADU``, ``HOGR``, ``LOON``, ``MAMU``, ``OLDS``, ``PIGU``,
        ``RBME``, ``RNGR``, ``SUSC``, ``WWSC``); ``NaN`` marks R's ``NA``.
    year, site, transect, temp : ndarray, shape (3793,)
        Survey year, site code, transect code, and temperature.
    observ_cond : list of str
        Observing conditions factor (5 levels).
    bay : list of str
        Bay name factor (4 levels).
    observ_cond_factor3 : list of str
        ``observ_cond`` collapsed to 3 levels.
    """

    counts: dict[str, NDArray[np.float64]]
    year: NDArray[np.float64]
    site: NDArray[np.float64]
    transect: NDArray[np.float64]
    temp: NDArray[np.float64]
    observ_cond: list[str]
    bay: list[str]
    observ_cond_factor3: list[str] = field(default_factory=list)


_SEABIRD_SPECIES = (
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
)


def load_seabird() -> Seabird:
    """Load :class:`Seabird` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> sb = fb.datasets.load_seabird()  # doctest: +SKIP
    >>> sb.year.shape  # doctest: +SKIP
    (3793,)
    """
    arrays, meta = _load("seabird")
    return Seabird(
        counts={species: arrays[species] for species in _SEABIRD_SPECIES},
        year=arrays["Year"],
        site=arrays["Site"],
        transect=arrays["Transect"],
        temp=arrays["Temp"],
        observ_cond=list(meta["ObservCond"]),
        bay=list(meta["Bay"]),
        observ_cond_factor3=list(meta["ObservCondFactor3"]),
    )


@dataclass(frozen=True, eq=False)
class ReginaPrecip:
    """Daily June precipitation (mm) in Regina, Saskatchewan, 1960-1993.

    Reproduces R ``fda::ReginaPrecip``: 16 missing values already omitted.

    Attributes
    ----------
    value : ndarray, shape (1006,)
        Precipitation in mm.
    """

    value: NDArray[np.float64]


def load_regina_precip() -> ReginaPrecip:
    """Load :class:`ReginaPrecip` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> rp = fb.datasets.load_regina_precip()  # doctest: +SKIP
    >>> rp.value.shape  # doctest: +SKIP
    (1006,)
    """
    arrays, _meta = _load("ReginaPrecip")
    return ReginaPrecip(value=arrays["value"])


@dataclass(frozen=True, eq=False)
class MontrealTemp:
    """Daily temperature (deg C) in Montreal, 1961-1994.

    Reproduces R ``fda::MontrealTemp``, a 34 (years) x 365 (days) matrix.

    Attributes
    ----------
    value : ndarray, shape (34, 365)
        Daily average temperature in degrees Celsius.
    years : list of str
        The 34 years, "1961".."1994".
    """

    value: NDArray[np.float64]
    years: list[str]


def load_montreal_temp() -> MontrealTemp:
    """Load :class:`MontrealTemp` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> mt = fb.datasets.load_montreal_temp()  # doctest: +SKIP
    >>> mt.value.shape  # doctest: +SKIP
    (34, 365)
    """
    arrays, meta = _load("MontrealTemp")
    return MontrealTemp(value=arrays["value"], years=list(meta["years"]))


@dataclass(frozen=True, eq=False)
class Daily:
    """Daily temperature/precipitation at 35 Canadian stations (alternate layout).

    Reproduces R ``fda::daily``: the same data as :class:`CanadianWeather`,
    kept by ``fda`` "primarily for compatibility with scripts written before
    the other format became available" and with the Matlab ``fda`` code --
    ``place`` here is padded to a fixed 11 characters (unlike
    :attr:`CanadianWeather.stations`, which is trimmed).

    Attributes
    ----------
    tempav, precav : ndarray, shape (365, 35)
        Average daily temperature (deg C) and precipitation (mm).
    place : list of str
        35 station names, each padded to 11 characters with trailing blanks.
    """

    tempav: NDArray[np.float64]
    precav: NDArray[np.float64]
    place: list[str]


def load_daily() -> Daily:
    """Load :class:`Daily` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> d = fb.datasets.load_daily()  # doctest: +SKIP
    >>> d.tempav.shape  # doctest: +SKIP
    (365, 35)
    """
    arrays, meta = _load("daily")
    return Daily(tempav=arrays["tempav"], precav=arrays["precav"], place=list(meta["place"]))


@dataclass(frozen=True, eq=False)
class InfantGrowth:
    """Tibia length for one infant's first 40 days of life.

    Reproduces R ``fda::infantGrowth`` (a 3-column matrix): each length is the
    average of 5 measurements.

    Attributes
    ----------
    day : ndarray, shape (40,)
        Age in days.
    tibia_length : ndarray, shape (40,)
        Average tibia length in mm.
    sd_length : ndarray, shape (40,)
        Standard deviation of the 5 measurements, in mm.
    """

    day: NDArray[np.float64]
    tibia_length: NDArray[np.float64]
    sd_length: NDArray[np.float64]


def load_infant_growth() -> InfantGrowth:
    """Load :class:`InfantGrowth` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> ig = fb.datasets.load_infant_growth()  # doctest: +SKIP
    >>> ig.day.shape  # doctest: +SKIP
    (40,)
    """
    arrays, _meta = _load("infantGrowth")
    return InfantGrowth(
        day=arrays["day"], tibia_length=arrays["tibia_length"], sd_length=arrays["sd_length"]
    )


@dataclass(frozen=True, eq=False)
class Nondurables:
    """US nondurable goods index, monthly, January 1919 onward.

    Reproduces R ``fda::nondurables``.

    Attributes
    ----------
    value : ndarray, shape (1377,)
        The index value each month.
    start : str
        First month covered, ``"1919-01"``.
    frequency : int
        Observations per year (12, monthly).
    """

    value: NDArray[np.float64]
    start: str
    frequency: int


def load_nondurables() -> Nondurables:
    """Load :class:`Nondurables` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> nd = fb.datasets.load_nondurables()  # doctest: +SKIP
    >>> nd.value.shape  # doctest: +SKIP
    (1377,)
    """
    arrays, meta = _load("nondurables")
    return Nondurables(
        value=arrays["value"], start=str(meta["start"]), frequency=int(meta["frequency"])
    )


@dataclass(frozen=True, eq=False)
class Lip:
    """51 measurements of lower lip position, 20 repetitions of "bob".

    Reproduces R ``fda::lip`` + ``liptime`` + ``lipmarks``.

    Attributes
    ----------
    value : ndarray, shape (51, 20)
        Lower lip position, sampled every 7 ms for 350 ms.
    t : ndarray, shape (51,)
        Time in seconds from the start (R: ``liptime``).
    left_elbow, right_elbow : ndarray, shape (20,)
        Landmark times (R: ``lipmarks$leftElbow``/``lipmarks$rightElbow``)
        for each of the 20 repetitions.
    """

    value: NDArray[np.float64]
    t: NDArray[np.float64]
    left_elbow: NDArray[np.float64]
    right_elbow: NDArray[np.float64]


def load_lip() -> Lip:
    """Load :class:`Lip` (downloaded + cached on first use).

    Examples
    --------
    >>> import fabel as fb
    >>> lip = fb.datasets.load_lip()  # doctest: +SKIP
    >>> lip.value.shape  # doctest: +SKIP
    (51, 20)
    """
    arrays, _meta = _load("lip")
    return Lip(
        value=arrays["value"],
        t=arrays["t"],
        left_elbow=arrays["left_elbow"],
        right_elbow=arrays["right_elbow"],
    )
