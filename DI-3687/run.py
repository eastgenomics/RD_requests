import logging
from pathlib import Path

from dx_search import load_config, search
from flagstat_plot import PLOT_METRICS, split_by_qc
from flagstat_reader import gather_flagstats
from plot import make_plots
from prepare_tables import prepare_tables
from qc_status import add_qc_status, gather_qc_status

CONFIGS = Path("dx_search/configs")
OUTPUT = Path("output")


def main():
    """Run the whole pipeline: search, read, prepare, join, plot and summarise (see the README for the order)"""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUTPUT.mkdir(exist_ok=True)

    flagstat_files = search(load_config(CONFIGS / "flagstat.py"))["files"]
    qc_files = search(load_config(CONFIGS / "qc_status.py"))["files"]

    # Any unreadable file raises. The raw tables are written first so they can be inspected if a later check fails
    df = gather_flagstats(flagstat_files)
    qc_df = gather_qc_status(qc_files)
    df.to_csv(OUTPUT / "flagstat_combined.tsv", sep="\t", index=False)
    qc_df.to_csv(OUTPUT / "qc_statuses_combined.tsv", sep="\t", index=False)

    flagstat_df, qc_df, removed = prepare_tables(df, qc_df)
    kept = add_qc_status(flagstat_df, qc_df)
    kept.to_csv(OUTPUT / "flagstat_qc_merged_filtered.tsv", sep="\t", index=False)
    removed.to_csv(OUTPUT / "removed_rows.tsv", sep="\t", index=False)

    plot_df = make_plots(kept)
    included, excluded = split_by_qc(plot_df)

    summary = [
        f"Projects: {df['project_id'].nunique()}",
        f"Rows read: {len(df):,}",
        *(
            f"Removed (see removed_rows.tsv), {reason}: {n}"
            for reason, n in removed["filter_reason"].value_counts().items()
        ),
        f"Samples after removal: {len(kept):,} ({kept['run'].nunique()} runs)",
        f"Included in statistics: {len(included):,}",
        *(
            f"Excluded from statistics (plotted as dots), {status}: {n}"
            for status, n in excluded["QC_Status"].value_counts().items()
        ),
        "",
    ]
    for column in PLOT_METRICS:
        values = included[column]
        stats = f"mean {values.mean():.2f}%, SD {values.std():.2f}, range {values.min():.2f} - {values.max():.2f}%"
        summary.append(f"{column}: {stats}")
    (OUTPUT / "summary.txt").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
