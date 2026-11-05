import io
import logging

import dxpy
import pandas as pd

from archival import raise_if_archived

logger = logging.getLogger(__name__)


def get_project_info(project_ids):
    """Map project IDs to their name and creation time"""
    return {
        project_id: dxpy.DXProject(project_id).describe(input_params={"fields": {"name": True, "created": True}})
        for project_id in project_ids
    }


def read_flagstat(file_id, project_id, sep="\t"):
    """Read one multiqc_samtools_flagstat.txt from DNAnexus into a df"""
    contents = dxpy.DXFile(file_id, project=project_id).read()
    return pd.read_csv(io.StringIO(contents), sep=sep)


def gather_flagstats(files):
    """
    Read and concatenate flagstat files into a single df

    Parameters
    ----------
    files : list
        File dicts from dx_search.search()["files"] (id, project, describe)

    Returns
    -------
    df : pd.DataFrame
        One row per sample with added run, project_created, project_id and file_id columns
    skipped : list
        (file_id, project_id, reason) for files that could not be read
    """
    raise_if_archived(files, "flagstat")

    project_info = get_project_info({f["project"] for f in files})
    dfs, skipped = [], []

    for file in files:
        file_id, project_id = file["id"], file["project"]

        try:
            df = read_flagstat(file_id, project_id)
        except (dxpy.exceptions.DXError, pd.errors.ParserError, pd.errors.EmptyDataError) as error:
            skipped.append((file_id, project_id, str(error)))
            continue

        df["run"] = project_info[project_id]["name"]
        df["project_created"] = pd.to_datetime(project_info[project_id]["created"], unit="ms")
        df["project_id"] = project_id
        df["file_id"] = file_id
        dfs.append(df)

    if not dfs:
        raise RuntimeError(f"No flagstat files could be read. Skipped: {skipped}")

    df = pd.concat(dfs, ignore_index=True)

    duplicated = df.duplicated(subset=["run", "Sample"], keep=False)
    if duplicated.any():
        logger.warning("%d rows share a run and Sample", duplicated.sum())

    for file_id, project_id, reason in skipped:
        logger.warning("Skipped %s (%s): %s", file_id, project_id, reason)

    return df, skipped
