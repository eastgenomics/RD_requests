import logging
from pathlib import Path

from qc_status import add_qc_status, gather_qc_status

from dx_search import load_config, search
from flagstat_filter import filter_flagstats
from flagstat_plot import PLOT_METRICS
from flagstat_reader import gather_flagstats
from plot import make_plots

CONFIGS = Path("dx_search/configs")
OUTPUT = Path("output")


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUTPUT.mkdir(exist_ok=True)

    flagstat_files = search(load_config(CONFIGS / "flagstat.py"))["files"]
    qc_files = search(load_config(CONFIGS / "qc_status.py"))["files"]

    df, _ = gather_flagstats(flagstat_files)
    qc_df, _ = gather_qc_status(qc_files)
    kept, removed = filter_flagstats(add_qc_status(df, qc_df))

    for name, table in {
        "flagstat_combined": df,
        "qc_status": qc_df,
        "flagstat_filtered": kept,
        "flagstat_removed": removed,
    }.items():
        table.to_csv(OUTPUT / f"{name}.tsv", sep="\t", index=False)

    plot_df = make_plots(kept)

    summary = [
        f"Projects: {df['project_id'].nunique()}",
        f"Rows read: {len(df):,}",
        *(f"Removed, {reason}: {n}" for reason, n in removed["filter_reason"].value_counts().items()),
        f"Samples after filtering: {len(kept):,} ({kept['run'].nunique()} runs)",
        "",
    ]
    for column in PLOT_METRICS:
        values = plot_df[column]
        stats = f"mean {values.mean():.2f}%, SD {values.std():.2f}, range {values.min():.2f} - {values.max():.2f}%"
        summary.append(f"{column}: {stats}")
    (OUTPUT / "summary.txt").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
