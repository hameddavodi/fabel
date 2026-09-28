"""Unit tests for fdatools._backend (array-API dispatch)."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fdatools import _backend as be

torch = pytest.importorskip("torch", reason="torch extra not installed")


# --------------------------------------------------------------------------- #
# default_namespace / array_namespace
# --------------------------------------------------------------------------- #


def test_default_namespace_is_numpy() -> None:
    xp = be.default_namespace()
    assert xp.__name__.endswith("numpy")


def test_array_namespace_with_no_arrays_falls_back_to_default() -> None:
    assert be.array_namespace() is be.default_namespace()
    assert be.array_namespace(1.0, 2, None) is be.default_namespace()


def test_array_namespace_from_list_is_default() -> None:
    assert be.array_namespace([1.0, 2.0]) is be.default_namespace()


def test_array_namespace_numpy() -> None:
    xp = be.array_namespace(np.zeros(3))
    assert xp.__name__.endswith("numpy")


def test_array_namespace_torch() -> None:
    xp = be.array_namespace(torch.zeros(3))
    assert "torch" in xp.__name__


def test_array_namespace_torch_wins_over_scalars() -> None:
    xp = be.array_namespace(1.0, torch.zeros(3), None)
    assert "torch" in xp.__name__


def test_array_namespace_mixed_backends_raises() -> None:
    with pytest.raises(TypeError):
        be.array_namespace(np.zeros(3), torch.zeros(3))


# --------------------------------------------------------------------------- #
# asarray
# --------------------------------------------------------------------------- #


def test_asarray_defaults_to_float64_numpy() -> None:
    a = be.asarray([1, 2, 3])
    assert isinstance(a, np.ndarray)
    assert a.dtype == np.float64


def test_asarray_preserves_torch() -> None:
    t = torch.tensor([1, 2, 3])
    a = be.asarray(t)
    assert isinstance(a, torch.Tensor)
    assert a.dtype == torch.float64


def test_asarray_explicit_namespace() -> None:
    xp = be.array_namespace(torch.zeros(1))
    a = be.asarray([1.0, 2.0], xp=xp)
    assert isinstance(a, torch.Tensor)


def test_asarray_dtype_none_keeps_dtype() -> None:
    t = torch.tensor([1.0, 2.0], dtype=torch.float32)
    assert be.asarray(t, dtype=None).dtype == torch.float32


def test_asarray_does_not_break_autograd() -> None:
    t = torch.tensor([1.0, 2.0], dtype=torch.float64, requires_grad=True)
    out = be.asarray(t)
    assert out.requires_grad
    out.sum().backward()
    assert t.grad is not None


# --------------------------------------------------------------------------- #
# to_numpy / is_torch
# --------------------------------------------------------------------------- #


def test_to_numpy_from_numpy() -> None:
    a = np.arange(3.0)
    assert be.to_numpy(a) is not None
    np.testing.assert_array_equal(be.to_numpy(a), a)


def test_to_numpy_from_torch() -> None:
    t = torch.tensor([1.0, 2.0], dtype=torch.float64)
    np.testing.assert_allclose(be.to_numpy(t), np.array([1.0, 2.0]))


def test_to_numpy_detaches_grad_tensor() -> None:
    t = torch.tensor([1.0, 2.0], dtype=torch.float64, requires_grad=True)
    out = be.to_numpy(t)
    assert isinstance(out, np.ndarray)


def test_to_numpy_from_list() -> None:
    np.testing.assert_allclose(be.to_numpy([1.0, 2.0]), np.array([1.0, 2.0]))


def test_is_torch() -> None:
    assert be.is_torch(torch.zeros(2))
    assert not be.is_torch(np.zeros(2))
    assert not be.is_torch(3.0)


# --------------------------------------------------------------------------- #
# same result numpy vs torch (WORKFLOW Phase 1 requirement, rtol 1e-10)
# --------------------------------------------------------------------------- #


@given(
    values=st.lists(
        st.floats(min_value=-50.0, max_value=50.0, allow_nan=False, allow_infinity=False),
        min_size=1,
        max_size=25,
    )
)
@settings(max_examples=50, deadline=None)
def test_numpy_torch_roundtrip_equal(values: list[float]) -> None:
    a = be.asarray(values)
    t = be.asarray(torch.tensor(values, dtype=torch.float64))
    np.testing.assert_allclose(be.to_numpy(a), be.to_numpy(t), rtol=1e-10, atol=0.0)


@given(
    values=st.lists(
        st.floats(min_value=-5.0, max_value=5.0, allow_nan=False, allow_infinity=False),
        min_size=1,
        max_size=20,
    )
)
@settings(max_examples=50, deadline=None)
def test_namespace_ops_agree_across_backends(values: list[float]) -> None:
    xn = be.array_namespace(np.zeros(1))
    xt = be.array_namespace(torch.zeros(1))
    a = be.asarray(values, xp=xn)
    b = be.asarray(values, xp=xt)
    got_n = be.to_numpy(xn.sum(xn.exp(a) * 2.0))
    got_t = be.to_numpy(xt.sum(xt.exp(b) * 2.0))
    np.testing.assert_allclose(got_n, got_t, rtol=1e-10, atol=0.0)
