# Flagstat search, filter and plot

Finds `multiqc_samtools_flagstat.txt` and `<run>_QC_Status.xlsx` in b38 CEN run projects (dxpy only), filters samples,
and plots supplementary alignments and mates mapped to a different chromosome as % of total reads per run.

## Setup and run

```bash
uv sync                  # Python 3.13 and the pinned packages in uv.lock
source .venv/bin/activate
dx login
python run.py
```

Dependencies and the Python version are pinned in [pyproject.toml](pyproject.toml) and [uv.lock](uv.lock). Format and
lint with `ruff format .` and `ruff check .` (config in `pyproject.toml`).

[run.py](run.py) runs the whole pipeline below and writes the tables, plots and `summary.txt` to `output/`. The other
modules are importable functions with no CLIs. To redraw the plots from the saved `output/flagstat_filtered.tsv` without
searching DNAnexus (e.g. after changing the plotting code), run `python plot.py` ([plot.py](plot.py)).

## 1. Search ([dx_search](dx_search))

[dx_search](dx_search) is generic: each config defines `project_search`, `file_search`, optional `project_filters` /
`file_filters` and `deduplicate`. The flagstat and QC status files are two separate searches.

Configs: [flagstat.py](dx_search/configs/flagstat.py) and [qc_status.py](dx_search/configs/qc_status.py). Both use the
same `project_search` (keep the two in step) and keep the newest file per project.

- `^002_26` keeps 2026 runs only and excludes projects renamed `003_`.
- The flagstat folder differs between runs (`eggd_MultiQC/multiqc_data`, `<run>-multiqc_data`, root), so it is not
  filtered on. The `test` filter drops ad-hoc copies such as `/test_themis_multiQC/...`.

## 2. Read

- [flagstat_reader.py](flagstat_reader.py): reads each flagstat file and concatenates, adding `run` (project name), `project_created`,
  `project_id`, `file_id`. Both readers raise an error, before reading anything, if any file is not `live` ([archival.py](archival.py)), listing every
  archived file, as all data for the projects found is needed. Unarchive them and rerun.
- [qc_status.py](qc_status.py): reads the xlsx in binary mode (`rb`), first sheet.
  - Headers are matched ignoring case, spaces and underscores, because the status header varies between runs
    (`QC_Status`, `QC Status`, `QC_status`). `Sample Name` and status are required, `Reason` is optional (blank header
    in one run, so its Reasons are not read).
  - Other columns and the run-level block to the right of the table are ignored.
  - `FAILED` is treated as `FAIL`.
  - Joined to the flagstat table on `project_id` + `Sample`. Samples not listed in a QC status file get no status.

## 3. Filter ([flagstat_filter.py](flagstat_filter.py))

Applied in this order:

1. **Duplicates:** key is the first two hyphen-separated fields of `Sample` (instrument ID-SPID). Keep the copy from
   the most recently created project.
2. **Controls:** remove samples with `Q` in the name (case sensitive, `control_pattern` argument).
3. **Failed QC:** remove `FAIL` and `CANCELLED`. Samples with no status are kept, as they cannot be judged.

Duplicates run first, so a sample whose newest copy failed would lose its older passing copy.

## 4. Plots ([flagstat_plot.py](flagstat_plot.py))

`% = count / flagstat_total * 100`, as MultiQC does for its flagstat "Percentage of total" table. `flagstat_total` is
passed + failed reads.

| Plot (`output/plots/`) | Count column |
|---|---|
| `supplementary_pct.html` | `supplementary_passed` |
| `mate_diff_chr_pct.html` | `with mate mapped to a different chr_passed` |
| `mate_diff_chr_mapq5_pct.html` | `with mate mapped to a different chr (mapQ >= 5)_passed` |

One box per run (by project creation date), every sample as a point (name and value on hover), overall mean and
±1 SD lines. No warning/fail thresholds. Metrics are set in `PLOT_METRICS` in [flagstat_plot.py](flagstat_plot.py).
