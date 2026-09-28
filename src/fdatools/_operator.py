"""The linear differential operator :class:`LDO`.

Lives in a private module rather than in :mod:`fdatools.core` because
:mod:`fdatools.basis` needs it too and ``core`` imports ``basis``.  It is
re-exported from :mod:`fdatools.core` and :mod:`fdatools`, which is where users see it.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from math import pi
from types import ModuleType
from typing import Any

from fdatools._backend import asarray

__all__ = ["LDO"]

Weight = float | Any  # a scalar, or any callable functional weight (an FData)


@dataclass(frozen=True, init=False)
class LDO:
    r"""Linear differential operator ``L x = Σ_j w_j D^j x + D^m x``.

    Replaces R's ``Lfd`` / ``int2Lfd`` / ``vec2Lfd``.  The operator is defined by
    its weight functions ``w_0 .. w_{m-1}``; the highest derivative ``D^m`` always
    enters with unit weight, so ``m = len(weights)``.

    Parameters
    ----------
    order : int, optional
        Shorthand for a plain derivative: ``LDO(2)`` is ``D²``, i.e. weights
        ``(0, 0)``.  Mutually exclusive with ``weights``.  Defaults to ``0``
        (the identity operator) when both arguments are omitted.
    weights : sequence, optional
        Weights ``w_0 .. w_{m-1}``.  Each entry is either a real scalar or a
        callable functional weight (an :class:`~fdatools.core.FData` holding a
        single curve).

    Attributes
    ----------
    weights : tuple
        The stored weights.
    order : int
        The highest derivative ``m``.

    Raises
    ------
    ValueError
        If both ``order`` and ``weights`` are given, or ``order`` is negative.

    Examples
    --------
    >>> from fdatools import LDO
    >>> LDO(2).order
    2
    >>> LDO(2).is_derivative
    True
    >>> LDO.harmonic(period=365.0).order
    3
    """

    weights: tuple[Weight, ...]

    def __init__(self, order: int | None = None, weights: Sequence[Weight] | None = None) -> None:
        if order is not None and weights is not None:
            raise ValueError("give either order or weights, not both")
        if weights is None:
            resolved = 0 if order is None else int(order)
            if resolved < 0:
                raise ValueError(f"order must be non-negative, got {resolved}")
            object.__setattr__(self, "weights", (0.0,) * resolved)
        else:
            object.__setattr__(self, "weights", tuple(weights))

    @classmethod
    def harmonic(cls, period: float) -> LDO:
        r"""Return the harmonic accelerator ``D³ x + ω² D x`` with ``ω = 2π/period``.

        Parameters
        ----------
        period : float
            Period of the sinusoid annihilated by the operator.

        Returns
        -------
        LDO
            Operator with weights ``(0, ω², 0)``.

        Examples
        --------
        >>> from fdatools import LDO
        >>> LDO.harmonic(period=1.0).weights[1] > 39.0
        True
        """
        omega = 2.0 * pi / float(period)
        return cls(weights=(0.0, omega * omega, 0.0))

    @property
    def order(self) -> int:
        """Highest derivative order ``m`` appearing in the operator.

        Returns
        -------
        int
            ``len(self.weights)``.

        Examples
        --------
        >>> from fdatools import LDO
        >>> LDO(weights=[1.0, 2.0]).order
        2
        """
        return len(self.weights)

    @property
    def is_derivative(self) -> bool:
        """Whether the operator is the plain derivative ``D^m``.

        Returns
        -------
        bool
            ``True`` when every weight is the scalar zero, so that closed-form
            derivative penalties apply.

        Examples
        --------
        >>> from fdatools import LDO
        >>> LDO(weights=[0.0, 1.0]).is_derivative
        False
        """
        return all(isinstance(w, (int, float)) and w == 0.0 for w in self.weights)

    @property
    def is_constant(self) -> bool:
        """Whether every weight is a plain number rather than a functional weight.

        Returns
        -------
        bool
            ``True`` when the operator is fully described by its numeric
            weights, and so compares and hashes by value.

        Examples
        --------
        >>> from fdatools import LDO
        >>> LDO(weights=[1.0, 2.0]).is_constant
        True
        """
        return all(isinstance(w, (int, float)) for w in self.weights)

    def apply(
        self,
        evaluate: Callable[[int], Any],
        t: Any,
        xp: ModuleType,
    ) -> Any:
        """Apply the operator to a family of functions evaluated at ``t``.

        Parameters
        ----------
        evaluate : callable
            ``evaluate(j)`` returns the ``j``-th derivative at ``t`` as an array
            whose first axis runs over ``t``.
        t : array
            Evaluation points, shape ``(n_t,)``.
        xp : module
            Array namespace to compute in.

        Returns
        -------
        array
            ``Σ_j w_j(t) · evaluate(j) + evaluate(m)``.

        Examples
        --------
        >>> import numpy as np
        >>> from fdatools import LDO
        >>> from fdatools._backend import default_namespace
        >>> t = np.linspace(0.0, 1.0, 3)
        >>> LDO(1).apply(lambda j: np.ones((3, 2)) * j, t, default_namespace()).tolist()
        [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]
        """
        out = evaluate(self.order)
        for j, weight in enumerate(self.weights):
            if isinstance(weight, (int, float)):
                if weight == 0.0:
                    continue
                out = out + float(weight) * evaluate(j)
                continue
            values = _weight_values(weight, t, xp)
            shape = (values.shape[0], *((1,) * (len(out.shape) - 1)))
            out = out + xp.reshape(values, shape) * evaluate(j)
        return out


def _weight_values(weight: Any, t: Any, xp: ModuleType) -> Any:
    """Evaluate a functional weight at ``t`` and flatten it to shape ``(n_t,)``."""
    values = asarray(weight(t), xp=xp)
    while len(values.shape) > 1:
        if values.shape[-1] != 1:
            raise ValueError("a functional LDO weight must hold exactly one curve")
        values = values[..., 0]
    return values
