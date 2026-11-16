import dxpy
import pandas as pd

from archival import raise_if_archived


def get_project_info(project_ids):
    """
    Map project IDs to their name and creation time.

    Parameters
    ----------
    project_ids : iterable of str
        DNAnexus project IDs

    Returns
    -------
    dict
        Project ID -> {"name": ..., "created": ...} (creation time in ms since the epoch)
    """
    return {
        project_id: dxpy.DXProject(project_id).describe(input_params={"fields": {"name": True, "created": True}})
        for project_id in project_ids
    }


def gather_files(files, file_to_df_reader, description):
    """
    Read every file with file_to_df_reader and concatenate the results.

    Parameters
    ----------
    files : list
        File dicts from dx_search.search()["files"] (id, project, describe)
    file_to_df_reader : callable(file, project_info)
        Takes one file dict and the info for that file's project ({"name", "created"}, from get_project_info) and
        returns a df
    description : str
        Name of the kind of file, used in error messages

    Returns
    -------
    pd.DataFrame

    Raises
    ------
    RuntimeError
        If there are no files, any file is archived, or any file cannot be read or parsed
    """
    if not files:
        raise RuntimeError(f"No {description} files to read")
    raise_if_archived(files, description)

    project_infos = get_project_info({file["project"] for file in files})

    dfs = []
    for file in files:
        project_info = project_infos[file["project"]]
        try:
            dfs.append(file_to_df_reader(file, project_info))
        except (dxpy.exceptions.DXError, ValueError) as error:
            raise RuntimeError(
                f"Could not read {description} file {file['id']} in {project_info['name']} ({file['project']}): {error}"
            ) from error

    return pd.concat(dfs, ignore_index=True)
