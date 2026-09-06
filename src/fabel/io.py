"""Interoperability with pandas, xarray and R's ``saveRDS`` format.

``pandas``, ``xarray`` and ``rdata`` are optional dependencies (extras
``pandas``/``io``): every function here imports them lazily and raises a clear
:class:`ImportError` if the extra is missing, so importing :mod:`fabel` itself
never requires them.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from fabel._backend import to_numpy
from fabel.basis import (
    Basis,
    BSpline,
    Constant,
    Exponential,
    Fourier,
    Monomial,
    Polygonal,
    Power,
)
from fabel.core import BiFData, FData

if TYPE_CHECKING:  # pragma: no cover - typing only
    import pandas as pd
    import xarray as xr

__all__ = ["LongData", "from_pandas", "read_rds", "to_pandas", "to_xarray"]


@dataclass(frozen=True, eq=False)
class LongData:
    """Per-curve arrays recovered from a long-format table.

    Attributes
    ----------
    t : dict of str to ndarray, or ndarray
        Argument values for each curve.  A ``dict`` keyed by curve id when
        curves do not share a grid (the irregular path); a single 1-D array,
        shared by every curve, when they do.
    y : dict of str to ndarray, or ndarray, shape (n_t, n_curves)
        Observed values, laid out the same way as ``t``: one array per curve
        id for the irregular path, or a single 2-D array (argument x curve)
        when the grid is shared.
    ids : list of str
        Curve ids, in first-seen order.
    """

    t: dict[str, NDArray[np.float64]] | NDArray[np.float64]
    y: dict[str, NDArray[np.float64]] | NDArray[np.float64]
    ids: list[str]


def from_pandas(df: pd.DataFrame, id_col: str, t_col: str, y_col: str) -> LongData:
    """Reshape a long-format :class:`pandas.DataFrame` into per-curve arrays.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per (curve, argument value) observation.
    id_col : str
        Column identifying which curve a row belongs to.
    t_col : str
        Column holding the argument value.
    y_col : str
        Column holding the observed value.

    Returns
    -------
    LongData
        ``t``/``y`` are a single shared 2-D array when every curve has the
        same, identically-ordered ``t_col`` values; otherwise a dict of
        per-curve 1-D arrays (curves may then have different lengths).

    Raises
    ------
    ValueError
        If any of ``id_col``, ``t_col``, ``y_col`` is missing from ``df``, or
        ``df`` has no rows.

    Examples
    --------
    >>> import pandas as pd
    >>> import fabel as fb
    >>> df = pd.DataFrame(
    ...     {"id": ["a", "a", "b", "b"], "t": [0.0, 1.0, 0.0, 1.0], "y": [1.0, 2.0, 3.0, 4.0]}
    ... )
    >>> long = fb.from_pandas(df, "id", "t", "y")
    >>> long.ids
    ['a', 'b']
    >>> long.t.shape
    (2, 2)
    """
    missing = [col for col in (id_col, t_col, y_col) if col not in df.columns]
    if missing:
        raise ValueError(f"columns {missing} not found in the data frame")
    if df.empty:
        raise ValueError("from_pandas requires at least one row")

    ids = list(dict.fromkeys(df[id_col]))
    per_curve_t: dict[str, NDArray[np.float64]] = {}
    per_curve_y: dict[str, NDArray[np.float64]] = {}
    for curve_id in ids:
        rows = df[df[id_col] == curve_id]
        per_curve_t[str(curve_id)] = rows[t_col].to_numpy(dtype=np.float64)
        per_curve_y[str(curve_id)] = rows[y_col].to_numpy(dtype=np.float64)

    str_ids = [str(curve_id) for curve_id in ids]
    lengths = {len(arr) for arr in per_curve_t.values()}
    shared_grid = len(lengths) == 1 and all(
        np.array_equal(per_curve_t[str_ids[0]], arr) for arr in per_curve_t.values()
    )
    if shared_grid:
        t_shared = per_curve_t[str_ids[0]]
        y_shared = np.column_stack([per_curve_y[curve_id] for curve_id in str_ids])
        return LongData(t=t_shared, y=y_shared, ids=str_ids)
    return LongData(t=per_curve_t, y=per_curve_y, ids=str_ids)


def to_pandas(fd: FData, t: Any) -> pd.DataFrame:
    """Evaluate ``fd`` on ``t`` and return a long-format :class:`pandas.DataFrame`.

    Parameters
    ----------
    fd : FData
        The curves to evaluate.
    t : array_like
        Argument values to evaluate at.

    Returns
    -------
    pandas.DataFrame
        One row per (t, curve[, var]) combination, with columns ``t``,
        ``curve`` (integer curve index), ``value`` and, when ``fd.n_vars > 1``,
        ``var``.

    Examples
    --------
    >>> import numpy as np
    >>> import fabel as fb
    >>> fd = fb.FData(np.eye(4), fb.BSpline(domain=(0.0, 1.0), n_basis=4))
    >>> df = fb.to_pandas(fd, np.array([0.0, 1.0]))
    >>> list(df.columns)
    ['t', 'curve', 'value']
    """
    import pandas as pd

    t_arr = np.asarray(t, dtype=np.float64)
    values = to_numpy(fd(t_arr))
    n_t = t_arr.shape[0]
    if fd.n_vars == 1:
        t_grid, curve_grid = np.meshgrid(t_arr, np.arange(fd.n_curves), indexing="ij")
        return pd.DataFrame(
            {
                "t": t_grid.ravel(),
                "curve": curve_grid.ravel(),
                "value": values.reshape(n_t, fd.n_curves).ravel(),
            }
        )
    t_grid3, curve_grid3, var_grid3 = np.meshgrid(
        t_arr, np.arange(fd.n_curves), np.arange(fd.n_vars), indexing="ij"
    )
    return pd.DataFrame(
        {
            "t": t_grid3.ravel(),
            "curve": curve_grid3.ravel(),
            "var": var_grid3.ravel(),
            "value": values.reshape(n_t, fd.n_curves, fd.n_vars).ravel(),
        }
    )


def to_xarray(fd: FData, t: Any) -> xr.DataArray:
    """Evaluate ``fd`` on ``t`` and return an :class:`xarray.DataArray`.

    Parameters
    ----------
    fd : FData
        The curves to evaluate.
    t : array_like
        Argument values to evaluate at.

    Returns
    -------
    xarray.DataArray
        Dims ``("t", "curve")``, or ``("t", "curve", "var")`` when
        ``fd.n_vars > 1``, coordinated on ``t``.

    Examples
    --------
    >>> import numpy as np
    >>> import fabel as fb
    >>> fd = fb.FData(np.eye(4), fb.BSpline(domain=(0.0, 1.0), n_basis=4))
    >>> da = fb.to_xarray(fd, np.array([0.0, 1.0]))
    >>> da.dims
    ('t', 'curve')
    """
    import xarray as xr

    t_arr = np.asarray(t, dtype=np.float64)
    values = to_numpy(fd(t_arr))
    if fd.n_vars == 1:
        return xr.DataArray(values, dims=("t", "curve"), coords={"t": t_arr})
    return xr.DataArray(values, dims=("t", "curve", "var"), coords={"t": t_arr})


def _rds_basis(obj: dict[str, Any]) -> Basis:
    """Build a :class:`~fabel.basis.Basis` from a converted R ``basisfd`` object."""
    dropind = obj.get("dropind")
    if dropind is not None and np.asarray(dropind).size:
        raise ValueError("a basis with dropind (dropped basis functions) is not supported")

    kind = str(np.asarray(obj["type"]).ravel()[0])
    domain = tuple(float(v) for v in np.asarray(obj["rangeval"]).ravel())
    if len(domain) != 2:
        raise ValueError(f"expected a 2-element rangeval, got {domain}")
    n_basis = int(np.asarray(obj["nbasis"]).ravel()[0])
    params_raw = obj.get("params")
    params = None if params_raw is None else np.asarray(params_raw, dtype=np.float64).ravel()

    if kind == "bspline":
        n_interior = 0 if params is None else params.size
        order = n_basis - n_interior
        interior = () if params is None else tuple(float(v) for v in params)
        breaks = (domain[0], *interior, domain[1])
        return BSpline(domain=domain, order=order, breaks=breaks)
    if kind == "fourier":
        period = float(params[0]) if params is not None and params.size else domain[1] - domain[0]
        return Fourier(domain=domain, n_basis=n_basis, period=period)
    if kind == "monomial":
        exponents = tuple(round(e) for e in params) if params is not None else tuple(range(n_basis))
        return Monomial(domain=domain, exponents=exponents)
    if kind == "power":
        exponents = tuple(float(e) for e in params) if params is not None else (0.0, 1.0)
        return Power(domain=domain, exponents=exponents)
    if kind == "expon":
        rates = tuple(float(r) for r in params) if params is not None else (0.0, 1.0)
        return Exponential(domain=domain, rates=rates)
    if kind == "const":
        return Constant(domain=domain)
    if kind == "polygonal":
        argvals = tuple(float(v) for v in params) if params is not None else domain
        return Polygonal(argvals=argvals)
    raise ValueError(f"unsupported R basis type {kind!r}")


def _basisfd_constructor(obj: dict[str, Any], attrs: Mapping[str, Any]) -> Basis:
    del attrs
    return _rds_basis(obj)


def _coefs_array(value: Any) -> NDArray[np.float64]:
    """Convert an ``fd``/``bifd`` coefficient field (ndarray or labelled DataArray)."""
    as_numpy = getattr(value, "to_numpy", None)
    array = as_numpy() if callable(as_numpy) else value
    return np.asarray(array, dtype=np.float64)


def _fd_constructor(obj: dict[str, Any], attrs: Mapping[str, Any]) -> FData:
    del attrs
    return FData(_coefs_array(obj["coefs"]), obj["basis"])


def _bifd_constructor(obj: dict[str, Any], attrs: Mapping[str, Any]) -> BiFData:
    del attrs
    return BiFData(_coefs_array(obj["coefs"]), obj["sbasis"], obj["tbasis"])


def read_rds(path: Any) -> Basis | FData | BiFData:
    """Read an R ``.rds`` file holding a ``basisfd``, ``fd`` or ``bifd`` object.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the ``.rds`` file, as written by R's ``saveRDS()``.

    Returns
    -------
    Basis or FData or BiFData
        Depending on which of the three R classes was stored.

    Raises
    ------
    ValueError
        If the stored R object's class is not one of ``basisfd``, ``fd`` or
        ``bifd``, or it uses a basis feature Fabel does not support
        (``dropind``, or a basis type outside
        bspline/fourier/monomial/power/expon/const/polygonal).

    Examples
    --------
    >>> import fabel as fb
    >>> basis = fb.read_rds("tests/fixtures/bspline_basis.rds")  # doctest: +SKIP
    """
    try:
        import rdata
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError("read_rds requires the 'rdata' package: pip install fabel[io]") from exc

    from rdata.conversion import DEFAULT_CLASS_MAP

    constructor_dict = {
        **DEFAULT_CLASS_MAP,
        "basisfd": _basisfd_constructor,
        "fd": _fd_constructor,
        "bifd": _bifd_constructor,
    }
    result = rdata.read_rds(path, constructor_dict=constructor_dict)
    if not isinstance(result, Basis | FData | BiFData):
        raise ValueError(
            f"expected an R basisfd/fd/bifd object, got a converted {type(result).__name__}"
        )
    return result
