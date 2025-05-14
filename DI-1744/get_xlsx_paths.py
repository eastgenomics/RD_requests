"""
This script retrieves and processes data from DNAnexus and Clarity.
"""
import pandas as pd
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

import argparse
import dxpy
import swifter


def parse_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--clarity_extract",
        type=str,
        help=("Path to Clarity extract file"),
    )

    return parser.parse_args()


def open_files(clarity):
    """
    Open files and read in file contents to DataFrames
    """
    with open(clarity) as f:
        clarity_df = pd.read_csv(f, delimiter=',')

    return clarity_df


def get_matching_projects():
    """
    Retrieve project IDs matching specific name patterns.
    This function searches for projects with names starting with '002' and ending with either 'CEN' or 'TWE'.

    Returns
    -------
    matching_projects_tuple_list: list
        A list of tuples containing project IDs and their names.
        Each tuple is in the format (project_id, project_name).
    """
    pattern = '^002.*_(CEN|TWE)$'
    matching_projects = list(dxpy.find_projects(
        name={'regexp': pattern}, describe=True))
    matching_projects_tuple_list = [(proj['id'], proj['describe']['name']) for proj in matching_projects]

    return matching_projects_tuple_list


def query_reports_for_project(project_id, sample_ids):
    """
    Query reports for a given project ID and sample IDs.
    This function retrieves file names matching specific patterns and returns a list of records.
    It handles errors for files not found or multiple matches.
    It also extracts the sample ID from the file name.

    Parameters
    ----------
    project_id : str
        The ID of the project to query.
    sample_ids : list
        List of sample IDs to search for in the project.

    Returns
    -------
    records: list
        A list of dictionaries containing sample ID, project ID, and file name.
    """
    records = []
    # pattern = rf'^\d+-({"|".join(sample_ids)})-[\w-]+_R\d+\.\d_(?:SNV|CNV)_1\.xlsx$'
    pattern = rf'.*({"|".join(sample_ids)}).*xlsx'
    try:
            matching_files = dxpy.find_data_objects(
                project=project_id,
                name=pattern,
                name_mode='regexp',
                describe={'fields': {'name': True}}
            )
            for file in matching_files:
                file_name = file['describe']['name']
                sample_id = file_name.split('-')[1]
                records.append({
                    'sample_id': sample_id,
                    'project_id': project_id,
                    'file_name': file_name
                })
    except Exception as e:
        print(
            f"Error fetching files for sample in project {project_id}: {e}")
    return records

def fetch_all_reports(df, chunk_size=100, max_workers=64):
    """
    Fetch all reports for the given DataFrame of sample IDs.
    This function retrieves project IDs matching specific name patterns,
    queries reports for each project, and merges the results with the original DataFrame.
    It also handles errors for samples found in multiple projects.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing sample IDs to search for.
    chunk_size : int, optional
        set number of ids to process at a time, by default 100
    max_workers : int, optional
        number of maximum workers to destribute across, by default 16

    Returns
    -------
    final_df : pd.DataFrame
        DataFrame containing the merged results with additional columns.
    """
    sample_ids = df['sample_id'].tolist()
    project_info = get_matching_projects()
    project_dict = dict(project_info)

    # Chunk sample IDs to manage search load
    chunks = [sample_ids[i:i + chunk_size]
              for i in range(0, len(sample_ids), chunk_size)]

    all_records = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for project_id in project_dict.keys():
            for chunk in chunks:
                futures.append(executor.submit(
                    query_reports_for_project, project_id, chunk))

        for future in as_completed(futures):
            try:
                records = future.result()
                all_records.extend(records)
            except Exception as e:
                print(f"Error in future: {e}")

    # Add project names to records
    for record in all_records:
        record['project_name'] = project_dict.get(
            record['project_id'], 'Unknown')

    # Create a DataFrame from the records
    records_df = pd.DataFrame(all_records)
    print(records_df.columns)
    # Merge with the original df to retain additional columns
    merged_df = pd.merge(df, records_df, on='sample_id', how='left')

    # Identify samples with multiple projects
    project_counts = merged_df.groupby('sample_id')['project_id'].nunique()
    multiple_projects = project_counts[project_counts > 1].index.tolist()

    # Log errors for samples with multiple projects
    for sample in multiple_projects:
        print(f"Error: Sample {sample} found in multiple projects.")

    # Filter out samples with multiple projects
    final_df = merged_df[~merged_df['sample_id'].isin(multiple_projects)]

    return final_df


def create_path(filename, assay, run):
    """
    Create a path to a specific file on clingen
    Parameters
    ----------
        assay (str): either CEN or WES
        filename (str): filename for the xlsx report
        run (str): sequencing run name
    Returns
    -------
        path (str): path to the given filename on clingen
    """
    base_path = "/appdata/clingen/cg/Regional Genetics Laboratories/Molecular Genetics/Data archive/Sequencing HT/"
    if assay == "CEN":
        path = base_path + f"{assay}/Run folders/{run}_{assay}/{filename}"
    elif assay == "WES":
        path = base_path + f"{assay}/{run}_TWE/{filename}"
    return path


def find_file_name(search_query):
    """
    Find the file name for a given search query on DNAnexus
    Parameters
    ----------
        search_query (str): search query for file name when searching DNAnexus
    Returns
    -------
        filename (str): a file name, or None if no or multiple matches found
    """
    files = list(dxpy.find_data_objects(
        name=search_query,
        name_mode="glob",
        describe=True
    ))
    for file in files:
        filenames = []
        filenames.append(file.get('describe').get('name'))
        if len(filenames) != 1:
            print(f"{search_query} returned multiple/no files")
            return None
        else:
            return filenames[0]


def main():
    """
    Script entry point
    """
    # Parse args
    args = parse_args()

    # Read data into dataframes
    clarity_df= open_files(args.clarity_extract)

    # Process data to construct a path for each specimen
    clarity_df = clarity_df.astype(str)
    # create df with column by splitting the Beaker Procedure Name to create a new column for assay
    clarity_df['Assay'] = clarity_df['Beaker Procedure Name'].str.split(' ').str[0]
    clarity_df['sample_id'] = clarity_df['Specimen Identifier'].str.split('-').str[1]

    report_df = fetch_all_reports(clarity_df)
    # Add the path to the report
    report_df['path'] = report_df.apply(
        lambda x: create_path(x['file_name'], x['Assay'], x['project_name']),
        axis=1
    )
    # Add specimen ID to the report_df
    report_df['specimen_id'] = report_df['file_name'].str.split('-').str[0]
    # Full sample ID
    report_df['full_sample_id'] = report_df['sample_id'].str.split('-').str[0] + "-" + report_df['sample_id'].str.split('-').str[1]
    # Create output files
    report_df.to_csv("all_data.csv")

if __name__ == "__main__":
    main()
