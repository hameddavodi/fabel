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

from fabel import datasets, density, profiling, sparse, stats
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
from fabel.density import DensityResult, IntensityResult, fit_density, fit_intensity
from fabel.dynamics import PDA, PDAStability, phase_plane
from fabel.io import from_pandas, read_rds, to_pandas, to_xarray
from fabel.profiling import ODEModel, ProfiledODE, profile_ode
from fabel.registration import Registrator, landmark_register, register
from fabel.regression import FRegress, LinmodResult, fregress, linmod
from fabel.smoothing import Smoother, SmoothResult, smooth
from fabel.sparse import PACE, SparseCov, sparse_cov, sparse_mean

if TYPE_CHECKING:
    from fabel import nn as nn

__all__ = [
    "FCCA",
    "FPCA",
    "LDO",
    "PACE",
    "PDA",
    "BSpline",
    "Basis",
    "BiFData",
    "Constant",
    "DensityResult",
    "Exponential",
    "FData",
    "FRegress",
    "Fourier",
    "IntensityResult",
    "LinmodResult",
    "Monomial",
    "ODEModel",
    "PDAStability",
    "Polygonal",
    "Power",
    "ProfiledODE",
    "Registrator",
    "SmoothResult",
    "Smoother",
    "SparseCov",
    "__version__",
    "datasets",
    "density",
    "fit_density",
    "fit_intensity",
    "fregress",
    "from_pandas",
    "inprod",
    "landmark_register",
    "linmod",
    "phase_plane",
    "profile_ode",
    "profiling",
    "read_rds",
    "register",
    "smooth",
    "sparse",
    "sparse_cov",
    "sparse_mean",
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
