# Data release (`data-v1`) — human step

`src/fabel/datasets.py` ships three small datasets in-package
(`growth`, `gait`, `pinch`, under `src/fabel/_data/`) and lazy-downloads the
rest from a GitHub release on first use, caching them at `FABEL_DATA_DIR`
(default `~/.cache/fabel/`). No dataset or model weight is committed to git —
every non-in-package dataset is referenced by URL + SHA-256 checksum only
(`_CHECKSUMS` in `datasets.py`).

This document is the one manual step in that pipeline: publishing the actual
`.npz`/`.json` files as GitHub release assets. Everything upstream of it is
already scripted.

## Pipeline (automated parts)

1. `Rscript tools/export_datasets.R` — dumps every real `fda` 6.3.0 dataset to
   `data_export/<name>.json` (gitignored, full R precision, `digits=17`).
2. `.venv/bin/python tools/build_data_release.py` — reads `data_export/`,
   writes `data_release/<name>.npz` (float64 arrays, `savez_compressed`) +
   `data_release/<name>.json` (string/label metadata) for every dataset, and
   prints a `name -> {npz, json}` SHA-256 table.
3. Copy that printed table into `_CHECKSUMS` in `src/fabel/datasets.py`
   (already done for the current dataset set; re-run only if a dataset's
   content changes).

`data_release/` is gitignored — it is a regenerable staging directory, not a
place datasets live long-term.

## Manual step: publish the GitHub release

1. Confirm `data_release/*.npz` and `data_release/*.json` are up to date
   (re-run steps 1-2 above if `data_export/` changed).
2. Create (or update) the GitHub release tagged **`data-v1`** on
   `AISMAsrl/fabel`, either via the GitHub UI or:
   ```
   gh release create data-v1 data_release/*.npz data_release/*.json \
     --title "fabel dataset assets v1" \
     --notes "Dataset arrays for fabel.datasets lazy loaders. See docs/dev/data-release.md."
   ```
   If the release already exists and only some files changed, use
   `gh release upload data-v1 data_release/<name>.{npz,json} --clobber`
   for just the changed files.
3. Verify `_download_file`'s URL template matches the uploaded asset names:
   `https://github.com/AISMAsrl/fabel/releases/download/data-v1/<name>.<ext>`
   (`_RELEASE_URL` in `datasets.py`).
4. Verify `_CHECKSUMS` in `datasets.py` matches the just-uploaded bytes exactly
   (re-run `build_data_release.py` and diff its printed table against
   `_CHECKSUMS` if in doubt) — a mismatch makes every affected `load_*()` call
   raise `ValueError: checksum mismatch` for every user, by design (fail
   loudly on data corruption, no silent coercion).

## Never commit

- `data_export/`, `data_release/` — both gitignored. Datasets are large,
  regenerable from R, and not source code.
- Anything under `~/.cache/fabel/` (or `FABEL_DATA_DIR`) — that is the local
  download cache, machine-specific.

## Testing without the network

Tests never hit the real release URL: `tests/unit/test_datasets.py` and
`tests/parity/test_datasets.py` monkeypatch `datasets._download_file` to copy
from the local `data_release/` staging directory instead (skipped entirely if
`data_release/` doesn't exist — run steps 1-2 above first). The one test that
does hit the real URL,
`test_datasets.py::test_real_download_from_github_release`, is marked
`@pytest.mark.network` and additionally skipped unless
`FABEL_RUN_NETWORK_TESTS=1` is set — run it once after publishing a release to
confirm the public URL and checksums actually work:

```
FABEL_RUN_NETWORK_TESTS=1 .venv/bin/pytest -m network tests/unit/test_datasets.py -v
```
