from pathlib import Path

import pandas as pd

from flagstat_plot import PLOT_METRICS, add_pct_columns, make_plot

OUTPUT = Path("output")


def make_plots(kept):
    """Write the percentage table and one plot per metric from the filtered flagstat table, returning the table"""
    (OUTPUT / "plots").mkdir(parents=True, exist_ok=True)
    plot_df = add_pct_columns(kept)
    plot_df.to_csv(OUTPUT / "plots" / "flagstat_with_pct.tsv", sep="\t", index=False)
    for column, metric in PLOT_METRICS.items():
        make_plot(plot_df, column, metric["title"]).write_html(OUTPUT / "plots" / f"{column}.html")
    return plot_df


if __name__ == "__main__":
    # Replot from the saved output of run.py without searching DNAnexus
    make_plots(pd.read_csv(OUTPUT / "flagstat_filtered.tsv", sep="\t", parse_dates=["project_created"]))
