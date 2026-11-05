import plotly.express as px

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
    """Add each metric's percentage of flagstat_total as a new column"""
    if (df["flagstat_total"] <= 0).any():
        raise ValueError("flagstat_total must be greater than 0 for every sample")

    df = df.copy()
    for pct_column, metric in metrics.items():
        df[pct_column] = df[metric["count_column"]] / df["flagstat_total"] * 100
    return df


def make_plot(df, column, title):
    """Box plot of a metric per run (ordered by project creation), with overall mean and +/-1 SD lines"""
    run_order = df.sort_values(["project_created", "run"])["run"].unique().tolist()

    fig = px.box(
        df,
        x="run",
        y=column,
        points="all",
        category_orders={"run": run_order},
        hover_data={"Sample": True, "run": True, column: ":.3f"},
    )

    mean, std = df[column].mean(), df[column].std()
    fig.add_hline(
        y=mean,
        line_color="green",
        annotation_text=f"<b>Mean: {mean:.3f}<br>SD: {std:.3f}</b>",
        annotation_position="right",
    )
    for sign, label in ((1, "+SD"), (-1, "-SD")):
        fig.add_hline(
            y=mean + sign * std,
            line_dash="dot",
            annotation_text=f"<b>{label}: {mean + sign * std:.3f}</b>",
            annotation_position="right",
        )

    fig.update_xaxes(tickangle=45, title_text="Run")
    fig.update_yaxes(title_text="% of total reads")
    fig.update_layout(title_text=f"{title} - {len(df)} samples, {len(run_order)} runs", title_x=0.5)
    return fig
