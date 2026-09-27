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

from typing import TYPE_CHECKING, Any

from fabel import datasets, stats
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
from fabel.decomposition import FCCA, FPCA
from fabel.dynamics import PDA, phase_plane
from fabel.io import from_pandas, read_rds, to_pandas, to_xarray
from fabel.registration import Registrator, landmark_register, register
from fabel.regression import FRegress, fregress
from fabel.smoothing import Smoother, SmoothResult, smooth

if TYPE_CHECKING:
    from fabel import nn as nn

__all__ = [
    "FCCA",
    "FPCA",
    "LDO",
    "PDA",
    "BSpline",
    "Basis",
    "BiFData",
    "Constant",
    "Exponential",
    "FData",
    "FRegress",
    "Fourier",
    "Monomial",
    "Polygonal",
    "Power",
    "Registrator",
    "SmoothResult",
    "Smoother",
    "__version__",
    "datasets",
    "fregress",
    "from_pandas",
    "inprod",
    "landmark_register",
    "phase_plane",
    "read_rds",
    "register",
    "smooth",
    "stats",
    "to_pandas",
    "to_xarray",
]

__version__ = "1.0.0"


def __getattr__(name: str) -> Any:
    """Import the optional :mod:`fabel.nn` module on first access.

    ``import fabel`` never imports PyTorch; ``fabel.nn`` loads it lazily.
    """
    if name == "nn":
        import importlib

        module = importlib.import_module("fabel.nn")
        globals()["nn"] = module
        return module
    raise AttributeError(f"module 'fabel' has no attribute {name!r}")
