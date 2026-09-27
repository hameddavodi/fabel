"""scikit-learn's own estimator checks on every Fabel estimator that takes ``X``.

SPEC 5.1: ``Smoother``, ``FPCA``, ``FCCA``, ``Registrator`` and ``FRegress``
implement the estimator API natively.  The four single-view estimators run the
full :func:`~sklearn.utils.estimator_checks.parametrize_with_checks` suite with
no expected failures.  ``FCCA.fit(X, y)`` takes two samples of *curves*
(``FData``), which the generic checks cannot build, so its API is covered in
``test_pipelines.py``.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import pytest
from sklearn.utils.estimator_checks import parametrize_with_checks

from fabel.decomposition import FPCA
from fabel.registration import Registrator
from fabel.regression import FRegress
from fabel.smoothing import Smoother

# sklearn's checks feed random blobs as coefficient rows; registering such
# curves to their mean often stops at ``max_iter`` and says so.  The checks
# test the API, not the registration, so that notice is expected here.
pytestmark = pytest.mark.filterwarnings("ignore:registration did not converge:RuntimeWarning")

ESTIMATORS = [Smoother(), FPCA(), FRegress(), Registrator()]


@parametrize_with_checks(ESTIMATORS)  # type: ignore[misc]
def test_sklearn_estimator_checks(estimator: Any, check: Callable[[Any], None]) -> None:
    """Every scikit-learn estimator check passes, with no exemptions."""
    check(estimator)
