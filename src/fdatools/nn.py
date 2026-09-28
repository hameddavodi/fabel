r"""PyTorch building blocks for functional data (the ``fdatools[torch]`` extra).

Three pieces connect fdatools' basis expansions to :mod:`torch`:

- :class:`BasisLayer` maps basis coefficients to curve values,
  ``x(t) = Φ(t) c`` -- a differentiable ``fd(t)`` for use inside a network,
  with gradients to the coefficients and, when the points are given at call
  time, to the points too.
- :class:`SmoothingLayer` maps raw observations on a grid to basis
  coefficients by penalised least squares, ``c = (ΦᵀWΦ + λR)⁻¹ΦᵀW y``, exactly
  as :func:`fdatools.smoothing.smooth` does; ``λ`` can be learned.
- :class:`FDataDataset` serves the curves of an :class:`~fdatools.core.FData`
  (coefficients, or values on a grid) and optional labels to a
  :class:`torch.utils.data.DataLoader`.

Basis and penalty matrices depend only on the basis and the grid, so they are
built once in float64 and stored as buffers: ``layer.to("cuda")`` or
``layer.to(torch.float32)`` moves them with the module.

This module imports PyTorch; ``import fdatools`` itself never does.

Examples
--------
>>> import numpy as np
>>> import torch
>>> import fdatools as fdt
>>> from fdatools.nn import BasisLayer, SmoothingLayer
>>> basis = fdt.BSpline(domain=(0.0, 1.0), n_basis=8)
>>> t = np.linspace(0.0, 1.0, 50)
>>> smoother, evaluate = SmoothingLayer(basis, t, lam=1e-6), BasisLayer(basis, t)
>>> y = torch.sin(2 * torch.pi * torch.tensor(t))[None, :]
>>> bool(torch.max(torch.abs(evaluate(smoother(y)) - y)) < 1e-2)
True
"""

from __future__ import annotations

from collections.abc import Sequence
from math import exp, log
from typing import Any

try:
    import torch
except ImportError as error:
    raise ImportError(
        "fdatools.nn requires PyTorch; install the torch extra: "
        'uv add "fdatools[torch] @ git+https://github.com/hameddavodi/fdatools"'
    ) from error

from fdatools import _linalg
from fdatools._backend import asarray, default_namespace, to_numpy
from fdatools._operator import LDO
from fdatools.basis import Basis
from fdatools.core import FData

__all__ = ["BasisLayer", "FDataDataset", "SmoothingLayer"]


def _buffer(values: Any, dtype: torch.dtype) -> torch.Tensor:
    """Return a NumPy-computed constant as a tensor of ``dtype`` on the CPU."""
    return torch.tensor(to_numpy(values), dtype=torch.float64).to(dtype)


def _grid(t: Any) -> Any:
    """Return evaluation points as a float64 NumPy vector."""
    points = asarray(to_numpy(t))
    if len(points.shape) != 1:
        raise ValueError(f"t must be one-dimensional, got shape {tuple(points.shape)}")
    return points


class BasisLayer(torch.nn.Module):
    r"""Evaluate basis expansions: map coefficients ``c`` to ``x(t) = Φ(t) c``.

    A differentiable ``fd(t)``: the input's last axis holds the ``n_basis``
    coefficients of one curve, the output's last axis the curve's values at
    the evaluation points.  The layer has no parameters of its own; put an
    :class:`torch.nn.Parameter` in front of it to learn a functional weight.

    Parameters
    ----------
    basis : Basis
        The expansion basis.
    t : array_like, optional
        Fixed evaluation points.  The basis matrix is precomputed and stored as
        the buffer ``design``.  When omitted, points must be passed to
        :meth:`forward`.
    deriv : int or LDO, optional
        Derivative order, or linear differential operator, to evaluate.
        Default ``0``.
    dtype : torch.dtype, optional
        Dtype of the stored basis matrix.  Default ``torch.float64``; use
        ``torch.float32`` on devices without double precision (Apple MPS).

    Examples
    --------
    >>> import numpy as np
    >>> import torch
    >>> import fdatools as fdt
    >>> from fdatools.nn import BasisLayer
    >>> layer = BasisLayer(fdt.BSpline(domain=(0.0, 1.0), n_basis=6), np.linspace(0.0, 1.0, 11))
    >>> layer(torch.ones(4, 6, dtype=torch.float64)).shape
    torch.Size([4, 11])
    """

    design: torch.Tensor | None

    def __init__(
        self,
        basis: Basis,
        t: Any = None,
        *,
        deriv: int | LDO = 0,
        dtype: torch.dtype = torch.float64,
    ) -> None:
        super().__init__()
        self.basis = basis
        self.deriv = deriv
        design = None if t is None else _buffer(basis(_grid(t), deriv), dtype)
        self.register_buffer("design", design)

    def extra_repr(self) -> str:
        """Describe the basis, the grid and the derivative."""
        points = "dynamic" if self.design is None else str(self.design.shape[0])
        return (
            f"basis={type(self.basis).__name__}(n_basis={self.basis.n_basis}), "
            f"n_points={points}, deriv={self.deriv!r}"
        )

    def _evaluate(self, t: Any, like: torch.Tensor) -> torch.Tensor:
        """Evaluate the basis at ``t`` (differentiably) and match ``like``'s device and dtype."""
        points = torch.as_tensor(t).to(device="cpu", dtype=torch.float64)
        design: torch.Tensor = self.basis(points, self.deriv)
        return design.to(device=like.device, dtype=like.dtype)

    def forward(self, coefs: torch.Tensor, t: Any = None) -> torch.Tensor:
        """Evaluate the curves whose coefficients are ``coefs``.

        Parameters
        ----------
        coefs : torch.Tensor
            Coefficients of shape ``(..., n_basis)``.
        t : array_like or torch.Tensor, optional
            Evaluation points of shape ``(n_t,)``, overriding the fixed grid.
            A tensor that requires grad receives gradients.

        Returns
        -------
        torch.Tensor
            Values of shape ``(..., n_t)``.

        Raises
        ------
        ValueError
            If there are no evaluation points, or the last axis of ``coefs``
            does not match the basis.

        Examples
        --------
        >>> import torch
        >>> import fdatools as fdt
        >>> from fdatools.nn import BasisLayer
        >>> layer = BasisLayer(fdt.Fourier(domain=(0.0, 1.0), n_basis=3))
        >>> t = torch.linspace(0.0, 1.0, 5, dtype=torch.float64)
        >>> layer(torch.ones(3, dtype=torch.float64), t).shape
        torch.Size([5])
        """
        if coefs.shape[-1] != self.basis.n_basis:
            raise ValueError(
                f"the last axis of coefs must hold {self.basis.n_basis} coefficients, "
                f"got {coefs.shape[-1]}"
            )
        if t is not None:
            design = self._evaluate(t, coefs)
        elif self.design is not None:
            design = self.design
        else:
            raise ValueError("no evaluation points: pass t to the constructor or to forward")
        return torch.matmul(coefs, design.transpose(0, 1))


class SmoothingLayer(torch.nn.Module):
    r"""Penalised least-squares smoothing: map observations to basis coefficients.

    Computes ``c = (ΦᵀWΦ + λR)⁻¹ ΦᵀW y`` for observations ``y`` on a fixed
    grid -- the fit of :func:`fdatools.smoothing.smooth` with a given ``λ`` --
    as a differentiable layer.  With a fixed ``λ`` this is one matrix product
    with a precomputed map; with ``trainable_lam=True`` the parameter
    ``log_lam = log λ`` is learned and the system is solved on every call.

    Parameters
    ----------
    basis : Basis
        The expansion basis.
    t : array_like
        Observation points, shape ``(n_t,)``.
    lam : float, optional
        Smoothing parameter ``λ ≥ 0`` (the initial value when trainable).
        Default ``0.0``.
    penalty : int or LDO, optional
        Roughness operator ``L`` of ``R = ∫ (Lφ)(Lφ)ᵀ``.  Default ``2``.
    weights : array_like, optional
        Observation weights, shape ``(n_t,)``.
    trainable_lam : bool, optional
        Learn ``log λ`` as a parameter.  Requires ``lam > 0``.
    dtype : torch.dtype, optional
        Dtype of the stored matrices and of ``log_lam``.  Default
        ``torch.float64``.

    Attributes
    ----------
    lam : float
        The current smoothing parameter.

    Raises
    ------
    ValueError
        If ``lam`` is negative (or not positive when trainable), or
        ``weights`` has the wrong shape.

    Examples
    --------
    >>> import numpy as np
    >>> import torch
    >>> import fdatools as fdt
    >>> from fdatools.nn import SmoothingLayer
    >>> t = np.linspace(0.0, 1.0, 30)
    >>> layer = SmoothingLayer(fdt.BSpline(domain=(0.0, 1.0), n_basis=8), t, lam=1e-4)
    >>> layer(torch.zeros(5, 30, dtype=torch.float64)).shape
    torch.Size([5, 8])
    """

    y2c: torch.Tensor | None
    gram: torch.Tensor | None
    roughness: torch.Tensor | None
    projection: torch.Tensor | None

    def __init__(
        self,
        basis: Basis,
        t: Any,
        *,
        lam: float = 0.0,
        penalty: int | LDO = 2,
        weights: Any = None,
        trainable_lam: bool = False,
        dtype: torch.dtype = torch.float64,
    ) -> None:
        super().__init__()
        if lam < 0.0:
            raise ValueError(f"lam must be non-negative, got {lam}")
        if trainable_lam and lam <= 0.0:
            raise ValueError("a trainable lam needs a positive initial value")
        xp = default_namespace()
        points = _grid(t)
        n_points = points.shape[0]
        w = xp.ones(n_points, dtype=xp.float64) if weights is None else asarray(to_numpy(weights))
        if tuple(w.shape) != (n_points,):
            raise ValueError(f"weights must have shape ({n_points},), got {tuple(w.shape)}")
        phi = asarray(basis(points))
        weighted = xp.matrix_transpose(w[:, None] * phi)
        gram = xp.matmul(weighted, phi)
        rough = asarray(basis.penalty(penalty))
        self.basis = basis
        self.n_points = n_points
        self.trainable_lam = trainable_lam
        if trainable_lam:
            self.log_lam = torch.nn.Parameter(torch.tensor(log(lam), dtype=dtype))
            self.register_buffer("y2c", None)
            self.register_buffer("gram", _buffer(gram, dtype))
            self.register_buffer("roughness", _buffer(rough, dtype))
            self.register_buffer("projection", _buffer(weighted, dtype))
        else:
            self._fixed_lam = float(lam)
            y2c = _linalg.solve_spd(gram + lam * rough, weighted)
            self.register_buffer("y2c", _buffer(y2c, dtype))
            self.register_buffer("gram", None)
            self.register_buffer("roughness", None)
            self.register_buffer("projection", None)

    @property
    def lam(self) -> float:
        """The current smoothing parameter ``λ``.

        Returns
        -------
        float
            ``exp(log_lam)`` for a trainable layer, else the fixed value.

        Examples
        --------
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.nn import SmoothingLayer
        >>> b = fdt.BSpline(n_basis=5)
        >>> SmoothingLayer(b, np.linspace(0, 1, 9), lam=0.5, trainable_lam=True).lam
        0.5
        """
        if self.trainable_lam:
            return float(exp(float(self.log_lam.detach().cpu())))
        return self._fixed_lam

    def extra_repr(self) -> str:
        """Describe the basis, the grid and the smoothing parameter."""
        return (
            f"basis={type(self.basis).__name__}(n_basis={self.basis.n_basis}), "
            f"n_points={self.n_points}, lam={self.lam:.6g}, trainable_lam={self.trainable_lam}"
        )

    def forward(self, y: torch.Tensor) -> torch.Tensor:
        """Smooth observations into basis coefficients.

        Parameters
        ----------
        y : torch.Tensor
            Observations of shape ``(..., n_t)``.

        Returns
        -------
        torch.Tensor
            Coefficients of shape ``(..., n_basis)``.

        Raises
        ------
        ValueError
            If the last axis of ``y`` does not match the observation grid.

        Examples
        --------
        >>> import numpy as np
        >>> import torch
        >>> import fdatools as fdt
        >>> from fdatools.nn import SmoothingLayer
        >>> layer = SmoothingLayer(fdt.Monomial(n_basis=2), np.linspace(0.0, 1.0, 5))
        >>> y = torch.linspace(1.0, 3.0, 5, dtype=torch.float64)  # y = 1 + 2t
        >>> [round(v, 10) for v in layer(y).tolist()]
        [1.0, 2.0]
        """
        if y.shape[-1] != self.n_points:
            raise ValueError(
                f"the last axis of y must hold the {self.n_points} observation points, "
                f"got {y.shape[-1]}"
            )
        if self.y2c is not None:
            y2c = self.y2c
        else:
            assert self.gram is not None
            assert self.roughness is not None
            assert self.projection is not None
            system = self.gram + torch.exp(self.log_lam) * self.roughness
            y2c = torch.linalg.solve(system, self.projection)
        return torch.matmul(y, y2c.transpose(0, 1))


class FDataDataset(torch.utils.data.Dataset[Any]):
    """Serve the curves of an :class:`~fdatools.core.FData` to a DataLoader.

    Item ``i`` is curve ``i`` -- its coefficients, or its values on a grid --
    optionally paired with its label.

    Parameters
    ----------
    fd : FData
        The curves.
    labels : array_like, optional
        One label per curve.  Numeric labels become a tensor; anything else
        (strings, say) is kept as a list.
    t : array_like, optional
        Evaluation grid.  When given, items are the curves' values at ``t``,
        shape ``(n_t,)`` (``(n_t, n_vars)`` for multivariate curves); otherwise
        they are coefficient vectors, shape ``(n_basis,)`` (``(n_basis,
        n_vars)``).
    dtype : torch.dtype, optional
        Dtype of the items.  Default ``torch.float64``.

    Raises
    ------
    TypeError
        If ``fd`` is not an :class:`~fdatools.core.FData`.
    ValueError
        If the number of labels differs from the number of curves.

    Examples
    --------
    >>> import numpy as np
    >>> import torch
    >>> import fdatools as fdt
    >>> from fdatools.nn import FDataDataset
    >>> fd = fdt.FData(np.eye(5), fdt.BSpline(n_basis=5))
    >>> ds = FDataDataset(fd, labels=[0, 1, 0, 1, 0])
    >>> x, y = next(iter(torch.utils.data.DataLoader(ds, batch_size=5)))
    >>> tuple(x.shape), y.tolist()
    ((5, 5), [0, 1, 0, 1, 0])
    """

    def __init__(
        self,
        fd: FData,
        labels: Sequence[Any] | Any = None,
        *,
        t: Any = None,
        dtype: torch.dtype = torch.float64,
    ) -> None:
        if not isinstance(fd, FData):
            raise TypeError(f"FDataDataset expects an FData, got {type(fd).__name__}")
        values = to_numpy(fd.coefs) if t is None else to_numpy(fd(_grid(t)))
        self.data = torch.as_tensor(values, dtype=torch.float64).movedim(1, 0).to(dtype)
        self.labels: torch.Tensor | list[Any] | None = None
        if labels is not None:
            if len(labels) != fd.n_curves:
                raise ValueError(f"got {len(labels)} labels for {fd.n_curves} curves")
            try:
                self.labels = torch.as_tensor(to_numpy(labels))
            except (TypeError, ValueError, RuntimeError):
                self.labels = list(labels)

    def __len__(self) -> int:
        """Return the number of curves."""
        return int(self.data.shape[0])

    def __getitem__(self, index: int) -> Any:
        """Return curve ``index``, with its label when the dataset has labels."""
        item = self.data[index]
        if self.labels is None:
            return item
        return item, self.labels[index]
