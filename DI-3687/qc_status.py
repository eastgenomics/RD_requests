import io
import logging

import dxpy
import pandas as pd

from archival import raise_if_archived

logger = logging.getLogger(__name__)

# Header names vary between runs (QC_Status, QC Status, QC_status), so headers are matched ignoring case,
# spaces and underscores. Reason is optional as some files leave its header blank.
REQUIRED_HEADERS = {"samplename": "Sample", "qcstatus": "QC_Status"}
OPTIONAL_HEADERS = {"reason": "Reason"}

# Statuses spelt differently between runs
STATUS_ALIASES = {"FAILED": "FAIL"}


def _normalise_header(header):
    return str(header).lower().replace(" ", "").replace("_", "")


def read_qc_status(file_id, project_id):
    """Read one <run>_QC_Status.xlsx from DNAnexus into a df of Sample, QC_Status, Reason and project_id"""
    contents = dxpy.open_dxfile(file_id, project=project_id, mode="rb").read()
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
    df["project_id"] = project_id
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
        Sample, QC_Status, Reason and project_id for every sample
    skipped : list
        (file_id, project_id, reason) for files that could not be read
    """
    raise_if_archived(qc_files, "QC status")

    dfs, skipped = [], []
    for file in qc_files:
        file_id, project_id = file["id"], file["project"]
        try:
            dfs.append(read_qc_status(file_id, project_id))
        except (dxpy.exceptions.DXError, ValueError) as error:
            skipped.append((file_id, project_id, str(error)))

    if not dfs:
        raise RuntimeError(f"No QC status files could be read. Skipped: {skipped}")

    qc_df = pd.concat(dfs, ignore_index=True)
    duplicated = qc_df.duplicated(["project_id", "Sample"], keep=False)
    if duplicated.any():
        raise ValueError(f"Samples listed more than once in a QC status file: {qc_df[duplicated]['Sample'].tolist()}")

    for file_id, project_id, reason in skipped:
        logger.warning("Skipped %s (%s): %s", file_id, project_id, reason)

    return qc_df, skipped


def add_qc_status(df, qc_df):
    """Left join QC_Status and Reason onto a flagstat table on project_id and Sample"""
    projects_without_qc = sorted(set(df["project_id"]) - set(qc_df["project_id"]))
    if projects_without_qc:
        logger.warning("%d projects have no QC status: %s", len(projects_without_qc), projects_without_qc)

    merged = df.merge(qc_df, on=["project_id", "Sample"], how="left", validate="m:1")
    logger.info(
        "QC status: %s; %d of %d samples have no status",
        merged["QC_Status"].value_counts().to_dict(),
        merged["QC_Status"].isna().sum(),
        len(merged),
    )
    return merged
