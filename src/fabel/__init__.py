"""Fabel: functional data analysis for Python.

A clean-room, production-grade Python implementation of the functional data
analysis toolkit described in Ramsay & Silverman, *Functional Data Analysis*.

Examples
--------
>>> import numpy as np
>>> import fabel as fb
>>> basis = fb.BSpline(domain=(0.0, 1.0), n_basis=5)
>>> basis(np.array([0.0, 0.5, 1.0])).shape
(3, 5)
"""

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
from fabel.core import LDO, BiFData, FData, inprod

__all__ = [
    "LDO",
    "BSpline",
    "Basis",
    "BiFData",
    "Constant",
    "Exponential",
    "FData",
    "Fourier",
    "Monomial",
    "Polygonal",
    "Power",
    "__version__",
    "inprod",
]

__version__ = "1.0.0"
