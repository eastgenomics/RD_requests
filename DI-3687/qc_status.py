import io
import logging

import dxpy
import pandas as pd

from gather import gather_files

logger = logging.getLogger(__name__)

# Header names vary between runs (QC_Status, QC Status, QC_status), so headers are matched ignoring case,
# spaces and underscores. Reason is optional as some files leave its header blank.
REQUIRED_HEADERS = {"samplename": "Sample", "qcstatus": "QC_Status"}
OPTIONAL_HEADERS = {"reason": "Reason"}

# Statuses spelt differently between runs
STATUS_ALIASES = {"FAILED": "FAIL"}


def _normalise_header(header):
    """
    Lower-case a header and remove spaces and underscores, so header spellings can be compared.

    Parameters
    ----------
    header : object
        Column header, converted to a string

    Returns
    -------
    str
    """
    return str(header).lower().replace(" ", "").replace("_", "")


def read_qc_status(file, project_info):
    """
    Read one <run>_QC_Status.xlsx from DNAnexus into a cleaned df of Sample, QC_Status, Reason, project_id and run
    (project name).

    Called by gather_files, which passes project_info: the name and creation time of the file's project, from
    gather.get_project_info. Only the name is used here.

    Parameters
    ----------
    file : dict
        File dict from dx_search.search()["files"] (id, project, describe)
    project_info : dict
        Name and creation time (ms since the epoch) of the file's project

    Returns
    -------
    pd.DataFrame
        Sample, QC_Status, Reason, project_id and run for every sample listed

    Raises
    ------
    ValueError
        If the Sample or QC_Status column cannot be found
    """
    contents = dxpy.open_dxfile(file["id"], project=file["project"], mode="rb").read()
    df = pd.read_excel(io.BytesIO(contents), engine="openpyxl")

    wanted = {**REQUIRED_HEADERS, **OPTIONAL_HEADERS}
    rename = {}
    for column in df.columns:
        target = wanted.get(_normalise_header(column))
        if target and target not in rename.values():
            rename[column] = target

    missing = set(REQUIRED_HEADERS.values()) - set(rename.values())
    if missing:
        raise ValueError(f"Columns not found: {sorted(missing)}. Headers present: {list(df.columns)}")

    df = df[list(rename)].rename(columns=rename).dropna(subset=["Sample"])
    if "Reason" not in df.columns:
        df["Reason"] = pd.NA

    status = df["QC_Status"].astype("string").str.strip().str.upper()
    df["QC_Status"] = status.replace(STATUS_ALIASES)
    df["project_id"] = file["project"]
    df["run"] = project_info["name"]
    return df


def gather_qc_status(qc_files):
    """
    Read and concatenate QC status files.

    Parameters
    ----------
    qc_files : list
        File dicts from dx_search.search()["files"] using configs/qc_status.py

    Returns
    -------
    qc_df : pd.DataFrame
        Sample, QC_Status, Reason, project_id and run (project name) for every sample

    Raises
    ------
    RuntimeError
        If there are no files, or any file cannot be read or parsed
    """
    return gather_files(qc_files, read_qc_status, "QC status")


def add_qc_status(df, qc_df):
    """
    Left join QC_Status and Reason onto a flagstat table on project_sample_key (see
    prepare_tables.add_project_sample_key). QC rows with no flagstat row are dropped.

    Parameters
    ----------
    df : pd.DataFrame
        Flagstat table with project_sample_key
    qc_df : pd.DataFrame
        QC status table with project_sample_key, QC_Status and Reason

    Returns
    -------
    pd.DataFrame
        df with QC_Status and Reason added. Samples with no QC row have no status

    Raises
    ------
    pandas.errors.MergeError
        If project_sample_key is not unique in either table
    """
    merged = df.merge(
        qc_df[["project_sample_key", "QC_Status", "Reason"]],
        on="project_sample_key",
        how="left",
        validate="1:1",
    )
    logger.info(
        "QC status: %s; %d of %d samples have no status",
        merged["QC_Status"].value_counts().to_dict(),
        merged["QC_Status"].isna().sum(),
        len(merged),
    )
    return merged
