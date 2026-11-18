import textwrap

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# QC statuses plotted as coloured dots and excluded from the statistics (boxes, mean and SD)
QC_COLOURS = {"FAIL": "red", "CANCELLED": "red", "WARNING": "orange"}

# Standard deviations from the mean to draw lines at
SD_LINES = (1, 2, 3)

# Percentage columns to add and plot: name -> count column. Percentages are count / flagstat_total * 100,
# the same denominator MultiQC uses for its "Percentage of total" flagstat table.
PLOT_METRICS = {
    "supplementary_pct": {
        "count_column": "supplementary_passed",
        "title": "Supplementary alignments (% of total reads)",
    },
    "mate_diff_chr_pct": {
        "count_column": "with mate mapped to a different chr_passed",
        "title": "Mate mapped to a different chromosome (% of total reads)",
    },
    "mate_diff_chr_mapq5_pct": {
        "count_column": "with mate mapped to a different chr (mapQ >= 5)_passed",
        "title": "Mate mapped to a different chromosome, mapQ >= 5 (% of total reads)",
    },
}


def add_pct_columns(df, metrics=PLOT_METRICS):
    """
    Add each metric's percentage of flagstat_total as a new column.

    Parameters
    ----------
    df : pd.DataFrame
        Flagstat table with run, Sample, flagstat_total and each metric's count column
    metrics : dict
        Maps each new percentage column name to a dict with its count_column and title (default PLOT_METRICS)

    Returns
    -------
    pd.DataFrame
        Copy of df with one percentage column per metric

    Raises
    ------
    ValueError
        If a required column is missing, any value is missing or not finite (naming the run, sample and column),
        or flagstat_total is 0 or less
    """
    required = ["flagstat_total", *(metric["count_column"] for metric in metrics.values())]
    missing_columns = [column for column in required if column not in df.columns]
    if missing_columns:
        raise ValueError(f"Flagstat columns not found: {missing_columns}")

    # Blank cells, or a column absent from some runs' files, would otherwise become NaN percentages that the
    # mean/SD silently skip while the samples are still counted
    numeric = df[required].apply(pd.to_numeric, errors="coerce")
    invalid = numeric.isna() | numeric.isin([float("inf"), float("-inf")])
    if invalid.any().any():
        problems = [
            f"{df.at[index, 'run']} {df.at[index, 'Sample']}: {column}"
            for column in required
            for index in df.index[invalid[column]]
        ]
        raise ValueError(f"Missing or non-finite flagstat values ({len(problems)}): {problems[:20]}")

    if (numeric["flagstat_total"] <= 0).any():
        raise ValueError("flagstat_total must be greater than 0 for every sample")

    df = df.copy()
    for pct_column, metric in metrics.items():
        df[pct_column] = numeric[metric["count_column"]] / numeric["flagstat_total"] * 100
    return df


def split_by_qc(df):
    """
    Split a flagstat table into samples included in the statistics and samples excluded from them.

    Samples with a QC_Status in QC_COLOURS (failed, cancelled, warning) are excluded from the statistics (boxes,
    mean and SD).
    Samples with no status are kept in the statistics, as they can't be judged.

    Parameters
    ----------
    df : pd.DataFrame
        Flagstat table with a QC_Status column

    Returns
    -------
    included : pd.DataFrame
        Samples included in the statistics (boxes, mean and SD)
    excluded : pd.DataFrame
        Samples excluded from the statistics (failed, cancelled and warning), plotted as coloured dots only
    """
    is_excluded = df["QC_Status"].isin(QC_COLOURS)
    return df[~is_excluded], df[is_excluded]


def _hover_text(df, status, width=60):
    """
    Hover text for samples excluded from the statistics: the sample name, the QC status and, where there is one, the
    QC reason (wrapped so a long reason stays on screen).

    Parameters
    ----------
    df : pd.DataFrame
        Excluded samples with the same QC_Status, with Sample and (optionally) Reason columns
    status : str
        The QC_Status of these samples
    width : int
        Characters per line before the reason is wrapped

    Returns
    -------
    pd.Series
        Hover text for each row, using <br> for line breaks
    """
    text = df["Sample"] + " (" + status + ")"
    if "Reason" not in df.columns:
        return text
    reason = df["Reason"].fillna("").astype(str).str.strip()
    wrapped = reason.map(lambda r: "<br>Reason: " + textwrap.fill(r, width).replace("\n", "<br>") if r else "")
    return text + wrapped


def make_plot(df, column, title):
    """
    Box plot of a metric per run (ordered by project creation), with overall mean and +/-1, 2 and 3 SD lines.
    The boxes, mean and SD use only the samples from split_by_qc's included set. Failed and cancelled samples are
    plotted as red dots and warning samples as orange dots on top.

    Parameters
    ----------
    df : pd.DataFrame
        Flagstat table with the percentage column, run, project_created, Sample, QC_Status and (optionally) Reason
    column : str
        Percentage column to plot
    title : str
        Plot title, followed by the sample and run counts

    Returns
    -------
    plotly.graph_objects.Figure
    """
    run_order = df.sort_values(["project_created", "run"])["run"].unique().tolist()
    included, excluded = split_by_qc(df)

    fig = px.box(
        included,
        x="run",
        y=column,
        points="all",
        category_orders={"run": run_order},
        hover_data={"Sample": True, "run": True, column: ":.3f"},
    )

    for status, status_df in excluded.groupby("QC_Status"):
        fig.add_trace(
            go.Box(
                x=status_df["run"],
                y=status_df[column],
                name=status,
                hovertext=_hover_text(status_df, status),
                hoverinfo="text+y",
                boxpoints="all",
                pointpos=0,
                jitter=0.5,
                line_width=0,
                fillcolor="rgba(0,0,0,0)",
                marker={"color": QC_COLOURS[status], "size": 7},
                legendgroup=status,
                showlegend=False,
            )
        )
        # A box trace's legend symbol is its invisible box, so add a marker-only trace for a coloured legend dot
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                name=status,
                legendgroup=status,
                marker={"color": QC_COLOURS[status], "size": 9},
                showlegend=True,
            )
        )

    mean, std = included[column].mean(), included[column].std()
    fig.add_hline(
        y=mean,
        line_color="green",
        annotation_text=f"<b>Mean: {mean:.3f}<br>SD: {std:.3f}</b>",
        annotation_position="right",
    )
    for n_sd in SD_LINES:
        for sign in (1, -1):
            fig.add_hline(
                y=mean + sign * n_sd * std,
                line_dash="dot",
                annotation_text=f"<b>{'+' if sign > 0 else '-'}{n_sd}SD: {mean + sign * n_sd * std:.3f}</b>",
                annotation_position="right",
            )

    excluded_text = f", {len(excluded)} warning/failed excluded from statistics" if len(excluded) else ""
    fig.update_xaxes(tickangle=45, title_text="Run")
    fig.update_yaxes(title_text="% of total reads")
    fig.update_layout(
        title_text=f"{title} - {len(included)} samples{excluded_text}, {len(run_order)} runs",
        title_x=0.5,
        boxmode="overlay",
    )
    return fig
