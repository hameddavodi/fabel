"""Array-API dispatch for Fabel.

This is the *only* module (together with :mod:`fabel._linalg`) that is allowed to
import NumPy directly.  Every other module obtains its array namespace through
:func:`array_namespace` and calls ``xp.*`` functions, so that a NumPy input
produces a NumPy output and a PyTorch input produces a PyTorch output with
gradients flowing through.

Examples
--------
>>> import numpy as np
>>> from fabel import _backend as be
>>> xp = be.array_namespace(np.zeros(3))
>>> be.asarray([1, 2, 3]).dtype
dtype('float64')
"""

from __future__ import annotations

from types import ModuleType
from typing import Any

import array_api_compat
import numpy as np

__all__ = [
    "array_namespace",
    "asarray",
    "default_namespace",
    "is_torch",
    "to_numpy",
]

_DEFAULT_NAMESPACE: ModuleType = array_api_compat.array_namespace(np.empty(0))


def default_namespace() -> ModuleType:
    """Return the fallback array namespace (NumPy).

    Returns
    -------
    module
        The array-API compatible NumPy namespace.

    Examples
    --------
    >>> from fabel._backend import default_namespace
    >>> default_namespace().__name__.endswith("numpy")
    True
    """
    return _DEFAULT_NAMESPACE


def array_namespace(*xs: Any) -> ModuleType:
    """Return the array-API namespace shared by ``xs``.

    Python scalars, ``None`` and plain sequences are ignored; if no argument is an
    array the NumPy namespace is returned.  If any argument is a
    :class:`torch.Tensor` the PyTorch namespace is returned, so tensors in give
    tensors out.

    Parameters
    ----------
    *xs : object
        Candidate arrays, scalars or ``None``.

    Returns
    -------
    module
        An array-API compatible namespace.

    Raises
    ------
    TypeError
        If ``xs`` mixes arrays from two different array libraries.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._backend import array_namespace
    >>> array_namespace(np.zeros(2), 1.0).__name__.endswith("numpy")
    True
    """
    arrays = [x for x in xs if array_api_compat.is_array_api_obj(x)]
    if not arrays:
        return _DEFAULT_NAMESPACE
    namespace: ModuleType = array_api_compat.array_namespace(*arrays)
    return namespace


def asarray(x: Any, xp: ModuleType | None = None, dtype: Any = float) -> Any:
    """Convert ``x`` to an array in namespace ``xp``.

    Parameters
    ----------
    x : array_like
        Value to convert.  Arrays are returned unchanged when the namespace and
        dtype already match, so autograd graphs are preserved.
    xp : module, optional
        Target namespace.  Defaults to the namespace of ``x``.
    dtype : dtype-like, optional
        Target dtype.  The default, the builtin :class:`float`, maps to the
        namespace's ``float64``.  Pass ``None`` to keep ``x``'s dtype.

    Returns
    -------
    array
        ``x`` as an array of namespace ``xp``.

    Examples
    --------
    >>> from fabel._backend import asarray
    >>> asarray([1, 2, 3]).dtype
    dtype('float64')
    """
    if xp is None:
        xp = array_namespace(x)
    resolved = xp.float64 if dtype is float else dtype
    if resolved is None:
        return xp.asarray(x)
    return xp.asarray(x, dtype=resolved)


def to_numpy(x: Any) -> np.ndarray[Any, np.dtype[Any]]:
    """Convert ``x`` to a NumPy array, detaching autograd tensors if needed.

    Parameters
    ----------
    x : array_like
        Array from any supported backend, or anything :func:`numpy.asarray`
        accepts.

    Returns
    -------
    numpy.ndarray
        A NumPy view or copy of ``x``.

    Examples
    --------
    >>> from fabel._backend import to_numpy
    >>> to_numpy([1.0, 2.0]).tolist()
    [1.0, 2.0]
    """
    if array_api_compat.is_torch_array(x):
        return np.asarray(x.detach().cpu().numpy())
    return np.asarray(x)


def is_torch(x: Any) -> bool:
    """Return whether ``x`` is a PyTorch tensor.

    Does not import PyTorch if it has not been imported already.

    Parameters
    ----------
    x : object
        Candidate value.

    Returns
    -------
    bool
        ``True`` if ``x`` is a :class:`torch.Tensor`.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._backend import is_torch
    >>> is_torch(np.zeros(2))
    False
    """
    return bool(array_api_compat.is_torch_array(x))
