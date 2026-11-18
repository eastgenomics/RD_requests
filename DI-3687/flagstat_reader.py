import io

import dxpy
import pandas as pd

from gather import gather_files


def read_flagstat(file, project_info):
    """
    Read one multiqc_samtools_flagstat.txt from DNAnexus into a df, adding run (project name), project_created,
    project_id and file_id.

    Called by gather_files, which passes project_info: the name and creation time of the file's project, from
    gather.get_project_info.

    Parameters
    ----------
    file : dict
        File dict from dx_search.search()["files"] (id, project, describe)
    project_info : dict
        Name and creation time (ms since the epoch) of the file's project

    Returns
    -------
    pd.DataFrame
        One row per sample, with run, project_created, project_id and file_id added
    """
    contents = dxpy.DXFile(file["id"], project=file["project"]).read()
    df = pd.read_csv(io.StringIO(contents), sep="\t")
    df["run"] = project_info["name"]
    df["project_created"] = pd.to_datetime(project_info["created"], unit="ms")
    df["project_id"] = file["project"]
    df["file_id"] = file["id"]
    return df


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

    Raises
    ------
    RuntimeError
        If there are no files, or any file cannot be read or parsed
    """
    return gather_files(files, read_flagstat, "flagstat")
