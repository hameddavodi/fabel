"""Unit tests for :mod:`fabel.io`: pandas/xarray round-trips and R ``.rds`` reading."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import fabel as fb
from fabel.basis import BSpline, Constant, Exponential, Fourier, Monomial, Polygonal, Power

pd = pytest.importorskip("pandas", reason="pandas extra not installed")
pytest.importorskip("xarray", reason="pandas extra (xarray) not installed")
pytest.importorskip("rdata", reason="io extra (rdata) not installed")

FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"


# --------------------------------------------------------------------------- #
# from_pandas
# --------------------------------------------------------------------------- #


def test_from_pandas_shared_grid_returns_arrays() -> None:
    df = pd.DataFrame(
        {"id": ["a", "a", "b", "b"], "t": [0.0, 1.0, 0.0, 1.0], "y": [1.0, 2.0, 3.0, 4.0]}
    )
    long = fb.from_pandas(df, "id", "t", "y")
    assert long.ids == ["a", "b"]
    assert isinstance(long.t, np.ndarray)
    assert long.t.shape == (2,)
    assert long.y.shape == (2, 2)
    np.testing.assert_array_equal(long.y[:, 0], [1.0, 2.0])
    np.testing.assert_array_equal(long.y[:, 1], [3.0, 4.0])


def test_from_pandas_irregular_grid_returns_dicts() -> None:
    df = pd.DataFrame({"id": ["a", "a", "b"], "t": [0.0, 1.0, 0.0], "y": [1.0, 2.0, 3.0]})
    long = fb.from_pandas(df, "id", "t", "y")
    assert isinstance(long.t, dict)
    assert isinstance(long.y, dict)
    np.testing.assert_array_equal(long.t["a"], [0.0, 1.0])
    np.testing.assert_array_equal(long.t["b"], [0.0])
    np.testing.assert_array_equal(long.y["b"], [3.0])


def test_from_pandas_missing_column_raises() -> None:
    df = pd.DataFrame({"id": ["a"], "t": [0.0], "y": [1.0]})
    with pytest.raises(ValueError, match="not found"):
        fb.from_pandas(df, "id", "missing", "y")


def test_from_pandas_empty_raises() -> None:
    df = pd.DataFrame({"id": [], "t": [], "y": []})
    with pytest.raises(ValueError, match="at least one row"):
        fb.from_pandas(df, "id", "t", "y")


# --------------------------------------------------------------------------- #
# to_pandas / to_xarray
# --------------------------------------------------------------------------- #


def test_to_pandas_univariate_columns_and_shape() -> None:
    fd = fb.FData(np.eye(4), BSpline(domain=(0.0, 1.0), n_basis=4))
    t = np.array([0.0, 0.5, 1.0])
    df = fd.to_pandas(t)
    assert list(df.columns) == ["t", "curve", "value"]
    assert len(df) == 3 * 4
    np.testing.assert_allclose(sorted(df["t"].unique()), t)


def test_to_pandas_multivariate_adds_var_column() -> None:
    coefs = np.random.default_rng(0).normal(size=(4, 3, 2))
    fd = fb.FData(coefs, BSpline(domain=(0.0, 1.0), n_basis=4))
    df = fd.to_pandas(np.array([0.0, 1.0]))
    assert list(df.columns) == ["t", "curve", "var", "value"]
    assert len(df) == 2 * 3 * 2


def test_to_xarray_univariate_dims() -> None:
    fd = fb.FData(np.eye(4), BSpline(domain=(0.0, 1.0), n_basis=4))
    t = np.array([0.0, 0.5, 1.0])
    da = fd.to_xarray(t)
    assert da.dims == ("t", "curve")
    assert da.shape == (3, 4)
    np.testing.assert_allclose(da.coords["t"].to_numpy(), t)


def test_to_xarray_multivariate_dims() -> None:
    coefs = np.random.default_rng(1).normal(size=(4, 3, 2))
    fd = fb.FData(coefs, BSpline(domain=(0.0, 1.0), n_basis=4))
    da = fd.to_xarray(np.array([0.0, 1.0]))
    assert da.dims == ("t", "curve", "var")
    assert da.shape == (2, 3, 2)


def test_to_pandas_and_to_xarray_agree_on_values() -> None:
    fd = fb.FData(np.eye(4), BSpline(domain=(0.0, 1.0), n_basis=4))
    t = np.array([0.0, 0.3, 1.0])
    df = fb.to_pandas(fd, t)
    da = fb.to_xarray(fd, t)
    for _, row in df.iterrows():
        expected = da.to_numpy()[int(np.where(t == row["t"])[0][0]), int(row["curve"])]
        assert row["value"] == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# read_rds
# --------------------------------------------------------------------------- #


def test_read_rds_bspline_basis() -> None:
    basis = fb.read_rds(FIXTURES / "bspline_basis.rds")
    assert isinstance(basis, BSpline)
    assert basis.domain == (0.0, 10.0)
    assert basis.n_basis == 7
    assert basis.order == 4


def test_read_rds_fourier_basis() -> None:
    basis = fb.read_rds(FIXTURES / "fourier_basis.rds")
    assert isinstance(basis, Fourier)
    assert basis.domain == (0.0, 12.0)
    assert basis.n_basis == 5
    assert basis.period == 12.0


def test_read_rds_fd() -> None:
    fd = fb.read_rds(FIXTURES / "bspline_fd.rds")
    assert isinstance(fd, fb.FData)
    assert fd.coefs.shape == (7, 3)
    assert isinstance(fd.basis, BSpline)


def test_read_rds_bifd() -> None:
    bifd = fb.read_rds(FIXTURES / "bspline_bifd.rds")
    assert isinstance(bifd, fb.BiFData)
    assert bifd.coefs.shape == (4, 4, 2, 1)


def test_read_rds_rejects_dropind() -> None:
    from fabel.io import _rds_basis

    obj = {
        "type": np.array(["bspline"]),
        "rangeval": np.array([0.0, 1.0]),
        "nbasis": np.array([4.0]),
        "params": np.array([0.5]),
        "dropind": np.array([1.0]),
    }
    with pytest.raises(ValueError, match="dropind"):
        _rds_basis(obj)


def test_read_rds_rejects_unsupported_basis_type() -> None:
    from fabel.io import _rds_basis

    obj = {
        "type": np.array(["nonexistent"]),
        "rangeval": np.array([0.0, 1.0]),
        "nbasis": np.array([4.0]),
        "params": None,
        "dropind": None,
    }
    with pytest.raises(ValueError, match="unsupported"):
        _rds_basis(obj)


def test_basis_roundtrip_matches_direct_construction() -> None:
    from fabel.io import _rds_basis

    for basis, obj in [
        (
            Monomial(domain=(0.0, 2.0), exponents=(0, 1, 2)),
            {
                "type": np.array(["monomial"]),
                "rangeval": np.array([0.0, 2.0]),
                "nbasis": np.array([3.0]),
                "params": np.array([0.0, 1.0, 2.0]),
                "dropind": None,
            },
        ),
        (
            Power(domain=(1.0, 2.0), exponents=(0.0, 0.5)),
            {
                "type": np.array(["power"]),
                "rangeval": np.array([1.0, 2.0]),
                "nbasis": np.array([2.0]),
                "params": np.array([0.0, 0.5]),
                "dropind": None,
            },
        ),
        (
            Exponential(domain=(0.0, 1.0), rates=(0.0, 1.0)),
            {
                "type": np.array(["expon"]),
                "rangeval": np.array([0.0, 1.0]),
                "nbasis": np.array([2.0]),
                "params": np.array([0.0, 1.0]),
                "dropind": None,
            },
        ),
        (
            Constant(domain=(0.0, 3.0)),
            {
                "type": np.array(["const"]),
                "rangeval": np.array([0.0, 3.0]),
                "nbasis": np.array([1.0]),
                "params": None,
                "dropind": None,
            },
        ),
        (
            Polygonal([0.0, 0.5, 1.0]),
            {
                "type": np.array(["polygonal"]),
                "rangeval": np.array([0.0, 1.0]),
                "nbasis": np.array([3.0]),
                "params": np.array([0.0, 0.5, 1.0]),
                "dropind": None,
            },
        ),
    ]:
        rebuilt = _rds_basis(obj)
        assert type(rebuilt) is type(basis)
        assert rebuilt.domain == basis.domain
        assert rebuilt.n_basis == basis.n_basis
