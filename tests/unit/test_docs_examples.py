"""The code in README.md and the user guide must run as written.

Each ```python block of a page is executed in order in one shared namespace, so
later blocks may use names from earlier ones. A block whose first line is
``# requires: <module>`` is skipped when that optional module is missing.
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PAGES = ["README.md", "docs/index.md", "docs/quickstart.md"]
BLOCK = re.compile(r"^```python\n(.*?)^```", re.DOTALL | re.MULTILINE)
REQUIRES = re.compile(r"^# requires: (\w+)")


def python_blocks(page: str) -> list[str]:
    return BLOCK.findall((REPO_ROOT / page).read_text(encoding="utf-8"))


def test_block_extraction() -> None:
    text = "intro\n```python\nx = 1\n```\n\n```bash\nls\n```\n```python\ny = x\n```\n"
    assert BLOCK.findall(text) == ["x = 1\n", "y = x\n"]


@pytest.mark.parametrize("page", PAGES)
def test_page_code_runs(page: str) -> None:
    blocks = python_blocks(page)
    assert blocks, f"{page} has no python examples"
    namespace: dict[str, object] = {"__name__": "__docs__"}
    for block in blocks:
        needed = REQUIRES.match(block)
        if needed and importlib.util.find_spec(needed.group(1)) is None:
            continue
        exec(compile(block, page, "exec"), namespace)


def test_quickstart_numbers() -> None:
    """The values the README and quickstart describe are the values the code gives."""
    namespace: dict[str, object] = {}
    for block in python_blocks("docs/quickstart.md"):
        if REQUIRES.match(block):
            continue
        exec(compile(block, "quickstart", "exec"), namespace)
    assert namespace["girls"].n_curves == 54  # type: ignore[attr-defined]
    assert namespace["heights"].shape == (3, 54)  # type: ignore[attr-defined]
    assert round(namespace["fixed_df"].df, 6) == 10.0  # type: ignore[attr-defined]
    assert namespace["pca"].scores.shape == (54, 3)  # type: ignore[attr-defined]
    varprop = np.asarray(namespace["pca"].varprop)  # type: ignore[attr-defined]
    assert np.all(np.diff(varprop) <= 0.0)
    assert namespace["basis"].penalty(namespace["accel"]).shape == (12, 12)  # type: ignore[attr-defined]
