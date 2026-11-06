import logging

import pandas as pd

logger = logging.getLogger(__name__)


def add_sample_key(df):
    """Add sample_key: the first two hyphen-separated fields of Sample (instrument ID-SPID)"""
    parts = df["Sample"].str.split("-")
    if parts.str.len().min() < 2:
        raise ValueError("Sample names must contain at least two hyphen-separated fields")
    df = df.copy()
    df["sample_key"] = parts.str[0] + "-" + parts.str[1]
    return df


def drop_duplicate_samples(df):
    """
    Keep one row per sample_key: the one from the project created most recently.
    Ties on creation time are broken by run name.

    Returns
    -------
    kept, removed : pd.DataFrame
    """
    ordered = df.sort_values(["project_created", "run"])
    is_duplicate = ordered.duplicated("sample_key", keep="last")
    return ordered[~is_duplicate].sort_index(), ordered[is_duplicate].sort_index()


def remove_controls(df, pattern="Q"):
    """
    Remove control samples, identified by pattern appearing in the Sample name (case sensitive, literal).

    Returns
    -------
    kept, removed : pd.DataFrame
    """
    is_control = df["Sample"].str.contains(pattern, regex=False)
    return df[~is_control], df[is_control]


def drop_failed_samples(df, failed_statuses=("FAIL", "CANCELLED")):
    """
    Remove samples whose QC_Status is a failure. Samples with no status are kept, as they can't be judged.

    Returns
    -------
    kept, removed : pd.DataFrame
    """
    if "QC_Status" not in df.columns:
        raise ValueError("QC_Status column not found - join QC status with qc_status.add_qc_status first")

    known = {"PASS", "WARNING", *failed_statuses}
    unexpected = set(df["QC_Status"].dropna().unique()) - known
    if unexpected:
        logger.warning("Unexpected QC_Status values (kept): %s", sorted(unexpected))

    is_failed = df["QC_Status"].isin(failed_statuses)
    return df[~is_failed], df[is_failed]


def filter_flagstats(df, control_pattern="Q"):
    """
    Apply the filtering rules in order, recording why rows were removed.

    Parameters
    ----------
    df : pd.DataFrame
        Flagstat table with QC_Status joined on
    control_pattern : str
        Text identifying control samples in the Sample name

    Returns
    -------
    kept : pd.DataFrame
    removed : pd.DataFrame
        Removed rows with a filter_reason column
    """
    df = add_sample_key(df)
    removed = []

    df, dropped = drop_duplicate_samples(df)
    removed.append(dropped.assign(filter_reason="duplicate sample (older project)"))

    df, dropped = remove_controls(df, control_pattern)
    removed.append(dropped.assign(filter_reason=f"control ('{control_pattern}' in name)"))

    df, dropped = drop_failed_samples(df)
    removed.append(dropped.assign(filter_reason="failed QC (" + dropped["QC_Status"].astype(str) + ")"))

    removed = pd.concat(removed).sort_index()
    for reason, count in removed["filter_reason"].value_counts().items():
        logger.info("Removed %d rows: %s", count, reason)
    logger.info("%d rows remaining", len(df))

    return df, removed
