"""Unit tests for fdatools._operator (the linear differential operator LDO)."""

from __future__ import annotations

from math import pi

import numpy as np
import pytest

from fdatools import LDO
from fdatools._backend import default_namespace


def test_order_shorthand_is_a_plain_derivative() -> None:
    op = LDO(2)
    assert op.order == 2
    assert op.weights == (0.0, 0.0)
    assert op.is_derivative


def test_default_is_the_identity_operator() -> None:
    assert LDO().order == 0
    assert LDO().weights == ()
    assert LDO().is_derivative


def test_weights_define_the_operator() -> None:
    op = LDO(weights=[1.0, 2.0])
    assert op.order == 2
    assert not op.is_derivative


def test_harmonic_accelerator_weights() -> None:
    op = LDO.harmonic(period=365.0)
    omega = 2.0 * pi / 365.0
    assert op.order == 3
    assert op.weights[0] == 0.0
    assert op.weights[1] == pytest.approx(omega**2)
    assert op.weights[2] == 0.0


def test_rejects_order_and_weights_together() -> None:
    with pytest.raises(ValueError, match="not both"):
        LDO(2, weights=[1.0])


def test_rejects_negative_order() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        LDO(-1)


def test_is_frozen_and_hashable() -> None:
    op = LDO(2)
    assert hash(op) == hash(LDO(2))
    assert op == LDO(2)
    assert op != LDO(3)
    with pytest.raises(AttributeError):
        op.weights = ()  # type: ignore[misc]


def test_apply_sums_weighted_derivatives() -> None:
    xp = default_namespace()
    t = np.linspace(0.0, 1.0, 4)

    def evaluate(j: int) -> np.ndarray:
        return np.full((4, 2), float(j))

    # L = 3 * D^0 + D^1  ->  3 * 0 + 1
    got = LDO(weights=[3.0]).apply(evaluate, t, xp)
    np.testing.assert_allclose(got, np.ones((4, 2)))


def test_apply_skips_zero_weights() -> None:
    xp = default_namespace()
    t = np.linspace(0.0, 1.0, 3)
    got = LDO(2).apply(lambda j: np.full((3, 1), float(j)), t, xp)
    np.testing.assert_allclose(got, np.full((3, 1), 2.0))


def test_apply_accepts_a_functional_weight() -> None:
    xp = default_namespace()
    t = np.linspace(0.0, 1.0, 5)
    op = LDO(weights=[lambda s: np.reshape(s, (-1, 1))])
    # L x = t * x + D x, with evaluate(j) == j
    got = op.apply(lambda j: np.full((5, 2), float(j)), t, xp)
    np.testing.assert_allclose(got, np.ones((5, 2)))


def test_functional_weight_must_be_a_single_curve() -> None:
    xp = default_namespace()
    t = np.linspace(0.0, 1.0, 4)
    op = LDO(weights=[lambda s: np.ones((len(s), 3))])
    with pytest.raises(ValueError, match="exactly one curve"):
        op.apply(lambda j: np.zeros((4, 1)), t, xp)
