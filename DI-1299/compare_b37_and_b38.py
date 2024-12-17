import argparse
import json
import re
import sys
from time import sleep

import dxpy
import pandas as pd
import plotly.express as px


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Information required to find Somalier files"
    )

    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=True,
        help="Path to the config file with variables to plot",
    )

    return parser.parse_args()


def read_in_json(file_name):
    """
    Read in JSON file to a dict

    Parameters
    ----------
    file_name : str
        name of JSON file to read in

    Returns
    -------
    json_dict : dict
        the JSON converted to a Python dictionary
    """
    with open(file_name, "r", encoding="utf8") as json_file:
        json_dict = json.load(json_file)

    return json_dict


def get_config_info(config_dict):
    """
    Get required keys from our config

    Parameters
    ----------
    config_dict : dict
        dict of JSON contents

    Returns
    -------
    assay : str
        The assay we're looking at
    search_term : str
        The search term to find projects in DX
    number_of_projects : int (optional)
        The number of projects to retrieve
    filename : str
        The name of the file to search for in the projects
    column_to_compare : str
        The column to compare between GRCh37 and GRCh38
    sample_column : str
        The column containing the sample IDs
    variables_to_plot : list
        List of lists of variables to plot
    """
    keys = [
        "assay",
        "search_term",
        "number_of_projects",
        "filename",
        "column_to_compare",
        "sample_column",
        "variables_to_plot",
    ]

    (
        assay,
        search_term,
        number_of_projects,
        filename,
        column_to_compare,
        sample_column,
        variables_to_plot,
    ) = list(map(config_dict.get, keys))

    return (
        assay,
        search_term,
        number_of_projects,
        filename,
        column_to_compare,
        sample_column,
        variables_to_plot,
    )


def find_projects(search_term, number_of_projects=None):
    """
    Find projects in DNAnexus using specified search term

    Parameters
    ----------
    search_term : str
        term to use to search DNAnexus projects
    number_of_projects : int, optional
        number of projects to retrieve - obtains the last n projects by name,
        default None

    Returns
    -------
    projects : list
        list of dicts, each representing a project
    """
    projects = list(
        dxpy.find_projects(
            name=search_term,
            name_mode="glob",
            describe={"fields": {"name": True}},
        )
    )
    if not projects:
        raise ValueError(
            f"No projects found with the search term {search_term}"
        )
    projects = sorted(projects, key=lambda x: x["describe"]["name"])

    if number_of_projects is not None:
        projects = projects[-int(number_of_projects) :]

    return projects


def find_files_in_project(search_term, project_id):
    """
    Find files in a DNAnexus project by name

    Parameters
    ----------
    search_term : str
        Term to use to find files
    project_id : str
        ID of the DX project to search in

    Returns
    -------
    files : list
        list of dicts, each representing a file in DNAnexus
    """
    files = list(
        dxpy.find_data_objects(
            project=project_id,
            name=search_term,
            name_mode="glob",
            classname="file",
            describe={"fields": {"name": True, "archivalState": True}},
        )
    )

    return files


def get_run_name_from_project_name(b38_project_name):
    """
    Extract the run name from a b38 project name

    Parameters
    ----------
    b38_project_name : str
        name of the GRCh38 project

    Returns
    -------
    match.group(0) : str
        name of the sequencing run
    """
    match = re.search(r"\d{6}_[A-Za-z0-9]+_\d{4}_[A-Z0-9]+", b38_project_name)
    if not match:
        raise ValueError(
            f"Error - no sequencing run name extracted from {b38_project_name}"
        )

    return match.group(0)


def read_dnanexus_file_to_df(file_id, project_id):
    """
    Reads the contents of TSV to pandas DataFrame and add project column

    Parameters
    ----------
    file_id : str
        DNAnexus file ID of report to read
    project_id : str
        DNAnexus project id

    Returns
    -------
    df : pd.DataFrame
        Pandas df containing content of report
    """
    with dxpy.open_dxfile(file_id, project=project_id) as dx_file:
        df = pd.read_csv(dx_file, sep="\t")
    df["project"] = project_id

    return df


def unarchive_non_live_files(file_list):
    """
    Unarchive all non-live files and exit

    Parameters
    ----------
    file_list : list
        list of dicts, each with info about a file
    """
    non_live_files = [
        file
        for file in file_list
        if file["describe"]["archivalState"] != "live"
    ]

    if non_live_files:
        for non_live_file in non_live_files:
            print(f"Requesting unarchiving for file {non_live_file['id']}")
            file_object = dxpy.DXFile(
                non_live_file["id"], project=non_live_file["project"]
            )
            file_object.unarchive()
            sleep(5)
        print("Exiting now. Please re-run once files are unarchived")
        sys.exit()
    else:
        print("All files are live")


def find_files_in_all_projects(b38_projects, filename, assay):
    """
    Find all files in GRCh37 and GRCh38 projects

    Parameters
    ----------
    b38_projects : list
        list of dicts, each with info about a DX GRCh38 project
    filename : str
        name of the file to search for in the projects
    assay : str
        the assay being investigated

    Returns
    -------
    b37_files : list
        list of dicts, each a file found in GRCh37
    b38_files : list
        list of dicts, each a file found in GRCh38
    """
    b37_files = []
    b38_files = []
    for proj in b38_projects:
        b38_file = find_files_in_project(filename, proj["id"])
        b38_files.extend(b38_file)
        run_name = get_run_name_from_project_name(proj["describe"]["name"])
        b37_project = find_projects(f"002_{run_name}_{assay}")
        b37_file = find_files_in_project(filename, b37_project[0]["id"])
        b37_files.extend(b37_file)

    return b37_files, b38_files


def read_in_b37_and_b38_files_and_merge(b38_files, b37_files):
    """
    Read in all GRCh37 files and all GRCh38 files to df then merge

    Parameters
    ----------
    b38_files : list
        list of dicts, each with info about a file in GRCh38
    b37_files : list
        list of dicts, each with info about a file in GRCh37

    Returns
    -------
    final_df : pd.DataFrame
        dataframe with results in GRCh38 and GRCh37 for each sample
    """
    all_b38_dfs = []
    for file in b38_files:
        all_b38_dfs.append(
            read_dnanexus_file_to_df(file["id"], file["project"])
        )
    b38_df = pd.concat(all_b38_dfs, ignore_index=True)

    all_b37_dfs = []
    for file in b37_files:
        all_b37_dfs.append(
            read_dnanexus_file_to_df(file["id"], file["project"])
        )
    b37_df = pd.concat(all_b37_dfs, ignore_index=True)

    all_results = pd.merge(
        b37_df,
        b38_df,
        on="sample_id",
        suffixes=("_GRCh37", "_GRCh38"),
    )

    # Reorder df so sample_id column is first + remove whitespace from headers
    columns_to_sort = [
        col for col in all_results.columns if col != "sample_id"
    ]
    reordered_df = all_results.reindex(
        ["sample_id"] + sorted(columns_to_sort), axis=1
    )
    final_df = reordered_df.rename(columns=lambda x: x.strip())

    return final_df


def write_to_file(df, filename):
    """
    Write a pandas DataFrame to a TSV file

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to write to file
    filename : str
        Name of the file to write to
    """
    df.to_csv(filename, sep="\t", index=False)


def convert_mismatches_to_wide_format(mismatches, sample_column):
    """
    Convert df to wide format and melt so we can plot easily with Plotly
    facets

    Parameters
    ----------
    mismatches : pd.DataFrame
        df with rows for any samples with mismatch in GRCh37 vs GRCh38
    sample_column : str
        name of the column containing sample IDs

    Returns
    -------
    df_melted : pd.DataFrame
        melted df with multiple rows per sample - all variables are in
        'variable' column and all values are in 'value' column
    """
    col_stubnames = list(
        set(
            col.removesuffix("_GRCh37").removesuffix("_GRCh38").strip()
            for col in mismatches.columns
            if col != sample_column
        )
    )

    wide_format = pd.wide_to_long(
        mismatches,
        i=sample_column,
        j="Genome",
        stubnames=col_stubnames,
        suffix="_(GRCh37|GRCh38)",
    ).reset_index()

    wide_format["Genome"] = wide_format["Genome"].apply(
        lambda s: s.removeprefix("_")
    )

    df_melted = wide_format.melt(
        id_vars=[sample_column, "Genome"],
        value_vars=col_stubnames,
        var_name="variable",
        value_name="value",
    )

    df_melted[sample_column] = (
        df_melted[sample_column].str.split("-").str[:2].str.join("-")
    )

    return df_melted


def create_and_save_scatter_plot(df, variable_set, file_name, sample_column):
    """
    Create and save a scatter plot for a given subset of variables.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to plot
    variable_set : list
        List of variables to include in the plot.
    file_name : str
        The name of the file to save the plot to
    sample_column : str
        the name of the column with sample ID
    """
    plot_subset = df.loc[df["variable"].isin(variable_set)]
    fig = px.scatter(
        plot_subset,
        x=sample_column,
        y="value",
        color="Genome",
        facet_row="variable",
    )
    fig.update_yaxes(matches=None)
    fig.update_xaxes(tickangle=45)
    fig.write_html(file_name)


def main():
    args = parse_args()
    config_dict = read_in_json(args.config)
    (
        assay,
        search_term,
        number_of_projects,
        filename,
        column_to_compare,
        sample_column,
        variables_to_plot,
    ) = get_config_info(config_dict)
    b38_projs = find_projects(search_term, number_of_projects)
    print(f"Found {len(b38_projs)} {assay} GRCh38 projects:")
    print("\n".join([proj["describe"]["name"] for proj in b38_projs]))

    b37_files, b38_files = find_files_in_all_projects(
        b38_projs, filename, assay
    )
    unarchive_non_live_files(b37_files + b38_files)

    merged_df = read_in_b37_and_b38_files_and_merge(b38_files, b37_files)
    write_to_file(merged_df, f"{assay}_all_results.tsv")
    print(f"There are {len(merged_df)} samples found for {assay}")

    mismatches = merged_df.query(
        f"{column_to_compare}_GRCh37 != {column_to_compare}_GRCh38"
    )

    print(
        f"There are {len(mismatches)} samples with mismatches in"
        f" {column_to_compare} between GRCh37 and GRCh38"
    )
    write_to_file(mismatches, f"{assay}_all_mismatches.tsv")

    if not mismatches.empty:
        mismatches_wide = convert_mismatches_to_wide_format(
            mismatches, sample_column
        )
        if variables_to_plot:
            for index, variable_list in enumerate(variables_to_plot):
                create_and_save_scatter_plot(
                    mismatches_wide,
                    variable_list,
                    f"{assay}_discrepancies_{index}.html",
                    sample_column,
                )
        else:
            print("No variables to plot in config file")


if __name__ == "__main__":
    main()
