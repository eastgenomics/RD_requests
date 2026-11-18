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
modules are importable functions with no CLIs. To redraw the plots from the saved `output/flagstat_qc_merged_filtered.tsv` without
searching DNAnexus (e.g. after changing the plotting code), run `python plot.py` ([plot.py](plot.py)).

## Order of execution

`python run.py` ([run.py](run.py) `main`) runs these steps in order. Any step that raises stops the run.

| # | Step | Where | Output |
|---|---|---|---|
| 1 | Search for flagstat files: the newest `multiqc_samtools_flagstat.txt` per project | [`dx_search.search`](dx_search/search.py) + [flagstat.py](dx_search/configs/flagstat.py) | |
| 2 | Search for QC status files: the newest root-level `_QC_Status.xlsx` per project | [`dx_search.search`](dx_search/search.py) + [qc_status.py](dx_search/configs/qc_status.py) | |
| 3 | Read and concatenate the flagstat files | [`gather_flagstats`](flagstat_reader.py) | |
| 4 | Read and concatenate the QC status files | [`gather_qc_status`](qc_status.py) | |
| 5 | Write the raw tables, so they can be inspected if a later step fails | [`run.py`](run.py) | `flagstat_combined.tsv`, `qc_statuses_combined.tsv` |
| 6 | Add `project_sample_key`, remove controls, check keys are unique, compare the two tables | [`prepare_tables`](prepare_tables.py) | log only |
| 7 | Left join QC status onto flagstat | [`add_qc_status`](qc_status.py) | `flagstat_qc_merged_filtered.tsv`, `removed_rows.tsv` |
| 8 | Calculate the percentages, check them, split the samples into included in and excluded from the statistics, and plot | [`plot.make_plots`](plot.py) | `plots/flagstat_with_pct.tsv`, `plots/*.html` |
| 9 | Write the summary | [`run.py`](run.py) | `summary.txt` |

Everything is written to `output/`. Details of each step are in sections 1 to 4 below.

## Errors and warnings

What stops the run and what is only logged:

| Situation | Outcome |
|---|---|
| A search finds no files | Error |
| Any file is archived | Error listing every archived file. Unarchive them and rerun |
| A file cannot be read or parsed | Error naming the file and project |
| The same sample twice within a project (`project_sample_key` repeated in a table) | Error naming the keys. The same sample in different projects is fine |
| A metric or `flagstat_total` is missing, not finite, or `flagstat_total` is 0 or less | Error naming the run, sample and column |
| Flagstat sample with no QC row | Warning listing the `project_sample_key`s. Kept with no status, so treated as passing |
| QC row with no flagstat sample | Warning listing the `project_sample_key`s. Dropped |
| `QC_Status` not `PASS`, `WARNING`, `FAIL` or `CANCELLED` | Warning. Treated as no status |
| Project with neither file | Not noticed |

## 1. Search ([dx_search](dx_search))

[dx_search](dx_search) is generic: each config defines `project_search`, `file_search`, optional `project_filters` /
`file_filters` and `deduplicate`. It only searches for files; `search()` logs a summary of the files, excluded files and
duplicate files at the end. The flagstat and QC status files are two separate searches.

Configs: [flagstat.py](dx_search/configs/flagstat.py) and [qc_status.py](dx_search/configs/qc_status.py), which share
the `project_search` in [common.py](dx_search/configs/common.py). Both keep the newest file per project. See the
configs for the search parameters.

## 2. Read

Both readers share [gather.py](gather.py), which checks the files are not archived ([archival.py](archival.py)),
looks up each project's name and creation time, reads every file and concatenates the results. Failures stop the run
(see the table above).

- [flagstat_reader.py](flagstat_reader.py): reads each flagstat file and concatenates, adding `run` (project name),
  `project_created`, `project_id`, `file_id`.
- [qc_status.py](qc_status.py): reads the xlsx in binary mode (`rb`), first sheet, adding `run` (project name) and
  `project_id`.
  - Headers are matched ignoring case, spaces and underscores, because the status header varies between runs
    (`QC_Status`, `QC Status`, `QC_status`). `Sample Name` and status are required, `Reason` is optional (blank header
    in one run, so its Reasons are not read).
  - Other columns and the run-level block to the right of the table are ignored.
  - `FAILED` is treated as `FAIL`.

A project with no flagstat or no QC file simply has no rows in that table. A project with neither file is not
noticed. A run that should not be analysed is excluded by editing the `dx_search` configs, not in the app.

The raw tables are written to `output/` (`flagstat_combined.tsv`, `qc_statuses_combined.tsv`) before any checks, so they can be
inspected if a check stops the run.

## 3. Prepare tables and join ([prepare_tables.py](prepare_tables.py))

Flagstat is the left table: it is the source of the metrics, and QC status is only attached to it. In order, on the
flagstat and QC status tables:

1. **Sample key:** `project_sample_key` is `project_id:` plus the first two hyphen-separated fields of `Sample`
   (instrument ID-SPID).
2. **Controls:** remove samples with `Q` in the name from both tables (case sensitive, `control_pattern` argument).
3. **Duplicates:** within a project, a sample listed more than once is an **error** (in either table). Across
   projects it is fine: a sample sequenced again in another project is a separate sequencing event, so **all copies are
   kept** and there is no sample-level deduplication.
4. **Comparison:** log (no error) how the `project_sample_key` sets overlap: counts in flagstat, QC status, both and
   either, then a warning listing the sorted `project_sample_key`s of flagstat samples with no QC row and of QC rows
   with no flagstat sample (one per line, so grouped by project). Unrecognised `QC_Status` values (not `PASS`, `WARNING`, `FAIL` or `CANCELLED`) also log a warning.
5. **Join:** left join of QC status onto flagstat on `project_sample_key`
   ([qc_status.py](qc_status.py) `add_qc_status`). QC rows with no flagstat row are dropped. Flagstat samples with no
   QC row are kept with no status and treated as passing.

Failed (`FAIL`, `CANCELLED`) and `WARNING` samples are **not removed**. They stay in `flagstat_qc_merged_filtered.tsv` with
their `QC_Status` so they can be plotted, but are **excluded from the statistics** in section 4. "Removed" is only used
for controls, which are dropped from the data and listed in `removed_rows.tsv`.

## 4. Plots ([flagstat_plot.py](flagstat_plot.py))

`% = count / flagstat_total * 100`, as MultiQC does for its flagstat "Percentage of total" table. `flagstat_total` is
passed + failed reads. The step stops with an error naming the run, sample and column if a metric or `flagstat_total` is
missing or not finite, or `flagstat_total` is 0 or less.

| Plot (`output/plots/`) | Count column |
|---|---|
| `supplementary_pct.html` | `supplementary_passed` |
| `mate_diff_chr_pct.html` | `with mate mapped to a different chr_passed` |
| `mate_diff_chr_mapq5_pct.html` | `with mate mapped to a different chr (mapQ >= 5)_passed` |

One box per run (by project creation date), every sample as a point (name and value on hover), overall mean and
±1, ±2 and ±3 SD lines (`SD_LINES`). `FAIL`/`CANCELLED` samples are plotted as red dots and `WARNING` samples as
orange dots (`QC_COLOURS`), with the sample name, status and QC `Reason` (when there is one) on hover, but are excluded
from the statistics (the boxes, mean and SD, `split_by_qc`), as are the
figures in `summary.txt`. Samples with no QC status are treated as passing. No warning/fail thresholds. Metrics are set in `PLOT_METRICS` in [flagstat_plot.py](flagstat_plot.py).
