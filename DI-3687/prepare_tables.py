import logging

logger = logging.getLogger(__name__)


def add_project_sample_key(df):
    """
    Add project_sample_key: project_id:<instrument ID>-<SPID>, where the instrument ID and SPID are the first two
    hyphen-separated fields of Sample. The same sample in another project is a separate sequencing event.

    Parameters
    ----------
    df : pd.DataFrame
        Table with Sample and project_id columns

    Returns
    -------
    pd.DataFrame
        Copy of df with project_sample_key added

    Raises
    ------
    ValueError
        If any Sample name has fewer than two hyphen-separated fields
    """
    parts = df["Sample"].str.split("-")
    if parts.str.len().min() < 2:
        raise ValueError("Sample names must contain at least two hyphen-separated fields")
    df = df.copy()
    df["project_sample_key"] = df["project_id"] + ":" + parts.str[0] + "-" + parts.str[1]
    return df


def remove_controls(df, pattern="Q"):
    """
    Remove control samples, identified by pattern appearing in the Sample name (case sensitive, literal).

    Parameters
    ----------
    df : pd.DataFrame
        Table with a Sample column
    pattern : str
        Text identifying control samples in the Sample name

    Returns
    -------
    kept : pd.DataFrame
        Rows that are not controls
    removed : pd.DataFrame
        Control rows
    """
    is_control = df["Sample"].str.contains(pattern, regex=False)
    return df[~is_control], df[is_control]


def check_unique_sample_keys(df, description):
    """
    Raise if the same project_sample_key appears more than once, as each should be listed once per project.

    Parameters
    ----------
    df : pd.DataFrame
        Table with a project_sample_key column
    description : str
        Name of the table, used in the error message

    Raises
    ------
    ValueError
        If any project_sample_key is repeated, naming up to the first 20
    """
    keys = df["project_sample_key"]
    duplicated = keys[keys.duplicated(keep=False)].unique().tolist()
    if duplicated:
        raise ValueError(
            f"{len(duplicated)} project_sample_key listed more than once in the {description} table: {duplicated[:20]}"
        )


def compare_sample_keys(df, qc_df):
    """
    Log how the project_sample_key sets of the flagstat and QC status tables overlap. Does not raise: flagstat
    samples with no QC row are kept with no status, and QC rows with no flagstat row are dropped by the join.

    Parameters
    ----------
    df : pd.DataFrame
        Flagstat table with project_sample_key
    qc_df : pd.DataFrame
        QC status table with project_sample_key
    """
    flagstat_keys, qc_keys = set(df["project_sample_key"]), set(qc_df["project_sample_key"])
    logger.info(
        "project_sample_key comparison: %d in flagstat, %d in QC status, %d in both, %d in either",
        len(flagstat_keys),
        len(qc_keys),
        len(flagstat_keys & qc_keys),
        len(flagstat_keys | qc_keys),
    )

    unmatched = (
        (flagstat_keys - qc_keys, "flagstat samples have no QC status row (kept, with no status)"),
        (qc_keys - flagstat_keys, "QC status rows have no flagstat sample (dropped)"),
    )
    for only_keys, description in unmatched:
        if only_keys:
            logger.warning("%d %s:\n  %s", len(only_keys), description, "\n  ".join(sorted(only_keys)))


def check_qc_statuses(df, known_statuses=("PASS", "WARNING", "FAIL", "CANCELLED")):
    """
    Warn about QC_Status values that are not recognised (they are treated as no status).

    Parameters
    ----------
    df : pd.DataFrame
        Table with a QC_Status column
    known_statuses : tuple of str
        QC_Status values that are expected
    """
    unexpected = set(df["QC_Status"].dropna().unique()) - set(known_statuses)
    if unexpected:
        logger.warning("Unexpected QC_Status values: %s", sorted(unexpected))


def prepare_tables(df, qc_df, control_pattern="Q"):
    """
    Prepare the flagstat and QC status tables for joining: add project_sample_key, remove controls from both, check
    each has one row per project_sample_key, and log how the two tables compare. Samples sequenced in more than one
    project are all kept, as each is a separate sequencing event. Failed and warning samples are kept too, so they
    can be plotted but excluded from the statistics (see flagstat_plot.split_by_qc).

    Parameters
    ----------
    df : pd.DataFrame
        Flagstat table
    qc_df : pd.DataFrame
        QC status table
    control_pattern : str
        Text identifying control samples in the Sample name

    Returns
    -------
    df, qc_df : pd.DataFrame
        The tables without controls, with project_sample_key added, ready for qc_status.add_qc_status
    removed : pd.DataFrame
        Flagstat rows removed, with a filter_reason column

    Raises
    ------
    ValueError
        If a project_sample_key appears more than once in either table, or a Sample name has fewer than two
        hyphen-separated fields
    """
    df, controls = remove_controls(add_project_sample_key(df), control_pattern)
    qc_df, _ = remove_controls(add_project_sample_key(qc_df), control_pattern)
    removed = controls.assign(filter_reason=f"control ('{control_pattern}' in name)")

    for table, description in ((df, "flagstat"), (qc_df, "QC status")):
        check_unique_sample_keys(table, description)
    compare_sample_keys(df, qc_df)
    check_qc_statuses(qc_df)

    for reason, count in removed["filter_reason"].value_counts().items():
        logger.info("Removed %d rows: %s", count, reason)
    logger.info("%d flagstat rows remaining", len(df))

    return df, qc_df, removed
