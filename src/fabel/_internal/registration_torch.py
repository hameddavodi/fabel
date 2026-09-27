r"""PyTorch path of :func:`fabel.registration.register`: autodiff Newton steps.

:func:`~fabel.registration.register` comes here when it receives curves (or a
continuous-registration target) with :class:`torch.Tensor` coefficients.  Two
pieces are provided:

- :class:`AutogradObjective` computes one curve's registration criterion in
  PyTorch (float64, CPU) and returns its gradient and Hessian by automatic
  differentiation.  :func:`fabel.registration._minimise` runs the same
  safeguarded Newton iteration with an Armijo line search on it as on the
  analytic NumPy derivatives, so both paths stop at the same optimum.
- :func:`to_torch_result` turns the finished registration into tensors in the
  input's dtype and device, with the registered curves differentiable with
  respect to the input coefficients.

The criterion is the one of :mod:`fabel.registration`: for the warp
parameters ``p`` (the free latent coefficients, plus the shift of a periodic
registration) it is the grid mean of ``(x_0 - x∘h)^2`` (least squares) or
twice the smaller eigenvalue of the grid-mean cross-product matrix of
``(x_0, x∘h)`` (eigen), plus ``λ cᵀRc``.  The smaller eigenvalue is computed
as the mean square of ``v_1 x_0 + v_2 x∘h`` with ``v`` the unit eigenvector
given by the rotation angle ``½ atan2(2b, a - d)``.  That expression is exact
for every ``p`` and free of the cancellation in ``(a + d)/2 - gap``, and
because ``v`` is built from ``p`` in the graph, autograd's Hessian includes
the eigenvector's own motion.  When the two eigenvalues coincide the angle is
held constant, which is the one-sided choice the analytic path makes too.

Gradients through the registered curves treat the optimal warps as fixed:
``registered`` is the least-squares projection of ``x_i(h_i(t))`` on the
registration grid, linear in the coefficients of ``x_i``, and that linear map
is what autograd sees.  The implicit dependence of the optimal ``h_i`` on
``x_i`` is not propagated.

This module imports PyTorch at the top; :mod:`fabel.registration` imports it
only when a tensor reaches :func:`~fabel.registration.register`, so
``import fabel`` never imports torch.

Examples
--------
>>> import numpy as np
>>> import torch
>>> from fabel import BSpline, FData
>>> from fabel.registration import register
>>> basis = BSpline(domain=(0.0, 1.0), n_basis=15)
>>> t = np.linspace(0.0, 1.0, 400)
>>> bumps = np.stack([np.exp(-(((t - s) / 0.1) ** 2)) for s in (0.45, 0.55)], axis=1)
>>> coefs = torch.tensor(np.linalg.lstsq(basis(t), bumps, rcond=None)[0], requires_grad=True)
>>> res = register(FData(coefs, basis), warp_basis=BSpline(n_basis=5), lam=1e-3)
>>> isinstance(res.registered.coefs, torch.Tensor)
True
>>> res.registered(t).sum().backward()
>>> coefs.grad.shape
torch.Size([15, 2])
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import torch
from fabel import _linalg
from fabel._backend import default_namespace, is_torch, to_numpy
from fabel.core import FData

if TYPE_CHECKING:
    from fabel.registration import RegistrationResult, _CurveProblem

__all__ = ["AutogradObjective", "as_tensor_like", "to_torch_result"]

Array = Any

#: The optimisation always runs in float64 on the CPU: registration stops at
#: gradients of order ``tol`` (1e-10 by default), and some accelerators (Apple
#: MPS) have no float64 at all.
_DTYPE = torch.float64


def _tensor(values: Any, dtype: torch.dtype = _DTYPE, device: Any = "cpu") -> torch.Tensor:
    """Copy a NumPy array (or array-like) into a new tensor of ``dtype`` on ``device``."""
    array = to_numpy(values).copy()
    return torch.from_numpy(array).to(dtype=dtype).to(device=device)


def as_tensor_like(values: Any, like: Any) -> torch.Tensor:
    """Return ``values`` as a constant tensor with the dtype and device of ``like``.

    Parameters
    ----------
    values : array_like
        Values to convert (NumPy or tensor; a tensor is detached and copied).
    like : torch.Tensor
        Reference tensor.

    Returns
    -------
    torch.Tensor
        A new tensor that does not track gradients.

    Examples
    --------
    >>> import numpy as np
    >>> import torch
    >>> from fabel._internal.registration_torch import as_tensor_like
    >>> as_tensor_like(np.arange(3.0), torch.zeros(1, dtype=torch.float64))
    tensor([0., 1., 2.], dtype=torch.float64)
    """
    return _tensor(values, like.dtype, like.device)


class AutogradObjective:
    """One curve's registration criterion with autograd derivatives.

    A drop-in replacement for :meth:`fabel.registration._CurveProblem.evaluate`:
    calling it with NumPy parameters returns the criterion as a float and its
    gradient and Hessian as NumPy arrays.  The criterion is evaluated in
    PyTorch; the gradient comes from one backward pass with
    ``create_graph=True`` and the Hessian from one more backward pass per
    parameter.

    Parameters
    ----------
    problem : _CurveProblem
        The curve, target values, warp quadrature, penalty and options of one
        curve, as built by :func:`fabel.registration.register`.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline, FData
    >>> from fabel.registration import _CurveProblem, _WarpQuadrature, _fine_grid
    >>> from fabel.registration import _penalty_matrix
    >>> from fabel._internal.registration_torch import AutogradObjective
    >>> curve = FData(np.linspace(0.0, 1.0, 6)[:, None] ** 2, BSpline(n_basis=6))
    >>> grid = _fine_grid((0.0, 1.0), 6)
    >>> wbasis = BSpline(n_basis=4)
    >>> problem = _CurveProblem(
    ...     curve=curve,
    ...     target=grid,
    ...     quadrature=_WarpQuadrature(wbasis, grid),
    ...     penalty=_penalty_matrix(wbasis, 2, 0.1),
    ...     lam=0.1,
    ...     criterion="least_squares",
    ...     periodic=False,
    ...     has_curvature=True,
    ... )
    >>> value, grad, hess = AutogradObjective(problem)(np.zeros(3))
    >>> auto = np.concatenate([[value], grad, hess.ravel()])
    >>> exact = problem.evaluate(np.zeros(3))
    >>> ref = np.concatenate([[exact[0]], exact[1], exact[2].ravel()])
    >>> bool(np.allclose(auto, ref, rtol=1e-12, atol=1e-14))
    True
    """

    def __init__(self, problem: _CurveProblem) -> None:
        quadrature = problem.quadrature
        self._phi = _tensor(quadrature.phi)
        self._weights = _tensor(quadrature.weights)
        self._index = torch.from_numpy(to_numpy(quadrature.index).astype("int64"))
        self._n_panels = quadrature.n_panels
        self._lower, self._upper = quadrature.span
        self._curve = FData(_tensor(problem.curve.coefs), problem.curve.basis)
        self._target = _tensor(problem.target)
        self._penalty = _tensor(to_numpy(problem.penalty)[1:, 1:])
        self._lam = float(problem.lam)
        self._eigen = problem.criterion == "eigen"
        self._periodic = problem.periodic
        self._free = problem.n_coefs - 1

    def criterion(self, params: torch.Tensor) -> torch.Tensor:
        """Evaluate the penalised registration criterion at ``params``.

        Parameters
        ----------
        params : torch.Tensor
            Free latent coefficients (the first one is pinned to zero),
            followed by the shift for a periodic registration; float64.

        Returns
        -------
        torch.Tensor
            The criterion as a 0-d tensor, differentiable in ``params``.
        """
        lower, upper = self._lower, self._upper
        free = self._free
        coefs = torch.cat([params.new_zeros(1), params[:free]])
        latent = torch.matmul(self._phi, coefs)
        # W only matters up to a constant: integrate exp(W - max W) so nothing
        # overflows.  The offset is detached; every ratio below is invariant.
        weighted = torch.exp(latent - torch.max(latent).detach()) * self._weights
        panels = torch.sum(torch.reshape(weighted, (self._n_panels, -1)), dim=1)
        cumulative = torch.cat([panels.new_zeros(1), torch.cumsum(panels, dim=0)])
        ratio = cumulative[self._index] / cumulative[-1]
        warp = torch.clamp(lower + (upper - lower) * ratio, lower, upper)
        if self._periodic:
            where = lower + torch.remainder(warp + params[free] - lower, upper - lower)
        else:
            where = warp
        values = self._curve(where)[:, 0]
        fit = self._eigen_fit(values) if self._eigen else self._least_squares_fit(values)
        rough = params[:free]
        return fit + self._lam * torch.dot(rough, torch.matmul(self._penalty, rough))

    def _least_squares_fit(self, values: torch.Tensor) -> torch.Tensor:
        """Grid mean of the squared distance to the target."""
        resid = self._target - values
        return torch.mean(resid * resid)

    def _eigen_fit(self, values: torch.Tensor) -> torch.Tensor:
        """Twice the smaller eigenvalue of the grid-mean cross-product matrix."""
        target = self._target
        aa = torch.mean(target * target)
        bb = torch.mean(target * values)
        dd = torch.mean(values * values)
        if float((0.25 * (aa - dd) ** 2 + bb * bb).detach()) > 0.0:
            angle = 0.5 * torch.atan2(2.0 * bb, aa - dd)
        else:
            # Equal eigenvalues: every direction is an eigenvector, and atan2
            # has no derivative at the origin.  Hold the angle constant.
            angle = torch.zeros((), dtype=_DTYPE)
        resid = -torch.sin(angle) * target + torch.cos(angle) * values
        return 2.0 * torch.mean(resid * resid)

    def __call__(self, params: Array) -> tuple[float, Array, Array]:
        """Return the criterion, its gradient and its Hessian at ``params``.

        Parameters
        ----------
        params : numpy.ndarray
            Parameter vector of shape ``(n_params,)``.

        Returns
        -------
        tuple
            ``(value, gradient, hessian)``: a float, a ``(n_params,)`` array
            and a ``(n_params, n_params)`` array, all NumPy float64.
        """
        xp = default_namespace()
        size = int(params.shape[0])
        point = _tensor(params).requires_grad_(True)
        value = self.criterion(point)
        if size == 0 or not value.requires_grad:
            return float(value.detach()), xp.zeros(size), xp.zeros((size, size))
        (grad,) = torch.autograd.grad(value, point, create_graph=True)
        rows = []
        for k in range(size):
            if grad.requires_grad:
                (row,) = torch.autograd.grad(
                    grad[k], point, retain_graph=True, allow_unused=True, materialize_grads=True
                )
            else:
                row = torch.zeros(size, dtype=_DTYPE)
            rows.append(row)
        hess = torch.stack(rows)
        return (
            float(value.detach()),
            xp.asarray(to_numpy(grad), dtype=xp.float64),
            xp.asarray(to_numpy(hess), dtype=xp.float64),
        )


def to_torch_result(
    result: RegistrationResult,
    fd: FData,
    target: FData | None,
    like: Any,
    *,
    periodic: bool,
) -> RegistrationResult:
    """Convert a finished registration into tensors, differentiable in ``fd``.

    ``registered`` is rebuilt as the least-squares projection (on the
    registration grid, as in :func:`fabel.registration.register`) of
    ``x_i(h_i(t))``, with the basis matrices at the warped times held
    constant and the coefficients of ``x_i`` taken from ``fd`` itself, so a
    loss on the registered curves back-propagates to ``fd.coefs``.  Every
    other field becomes a constant tensor.

    Parameters
    ----------
    result : RegistrationResult
        The registration computed on the NumPy copy of ``fd``.
    fd : FData
        The curves as passed to :func:`~fabel.registration.register`.
    target : FData or None
        The target as passed to :func:`~fabel.registration.register`.
    like : torch.Tensor
        Tensor whose dtype and device the outputs take.
    periodic : bool
        Whether the warps are periodic (warped times wrap around the domain).

    Returns
    -------
    RegistrationResult
        The same registration with tensor fields.

    Examples
    --------
    >>> import numpy as np
    >>> import torch
    >>> from fabel import BSpline, FData
    >>> from fabel.registration import register
    >>> from fabel._internal.registration_torch import to_torch_result
    >>> fd = FData(np.random.default_rng(0).standard_normal((8, 2)), BSpline(n_basis=8))
    >>> res = register(fd, landmarks=[0.45, 0.55])
    >>> out = to_torch_result(res, fd, None, torch.zeros(1, dtype=torch.float64), periodic=False)
    >>> type(out.warp_inverse.coefs).__name__, out.registered.coefs.dtype
    ('Tensor', torch.float64)
    >>> bool(torch.allclose(out.registered.coefs, torch.tensor(res.registered.coefs)))
    True
    """
    from fabel.registration import RegistrationResult, _fine_grid

    xp = default_namespace()
    dtype, device = like.dtype, like.device

    def constant(values: Any) -> torch.Tensor:
        return _tensor(values, dtype, device)

    def constant_fd(values: FData) -> FData:
        return FData(constant(values.coefs), values.basis)

    basis = fd.basis
    lower, upper = fd.domain
    grid = _fine_grid(fd.domain, basis.n_basis)
    where = to_numpy(result.warp_values(grid))
    if periodic:
        where = lower + xp.remainder(where - lower, upper - lower)
    n_curves = fd.n_curves
    design = xp.stack([basis(where[:, i]) for i in range(n_curves)])
    projection = _linalg.lstsq(basis(grid), xp.eye(int(grid.shape[0])))
    if is_torch(fd.coefs):
        curves = fd
        coefs = torch.reshape(fd.coefs, (basis.n_basis, n_curves)).to(dtype=dtype)
    else:
        curves = constant_fd(fd)
        coefs = constant(xp.reshape(fd.coefs, (basis.n_basis, n_curves)))
    values = torch.einsum("igk,ki->gi", constant(design), coefs)
    registered = FData(torch.matmul(constant(projection), values), basis)

    goal: FData | None = None
    if result.target is not None:
        if isinstance(target, FData) and is_torch(target.coefs):
            goal = target
        elif target is None and is_torch(fd.coefs):
            goal = fd.mean()
        else:
            goal = constant_fd(result.target)
    return RegistrationResult(
        registered=registered,
        warp=constant_fd(result.warp),
        unregistered=curves,
        latent=constant_fd(result.latent),
        shift=constant(result.shift),
        target=goal,
        warp_inverse=None if result.warp_inverse is None else constant_fd(result.warp_inverse),
        criterion=None if result.criterion is None else constant(result.criterion),
        n_iter=(
            None
            if result.n_iter is None
            else torch.from_numpy(to_numpy(result.n_iter).astype("int64")).to(device=device)
        ),
    )
