"""
This script retrieves and processes data from DNAnexus and Clarity.
"""
import pandas as pd
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

import argparse
import dxpy
import swifter

# Function to retrieve project IDs matching specific name patterns
def get_matching_projects():
    pattern = '^002.*_(CEN|TWE)$'
    matching_projects = list(dxpy.find_projects(name={'regexp': pattern}, describe=True))
    return [(proj['id'], proj['describe']['name']) for proj in matching_projects]

# Function to search for files matching the sample ID pattern in a given project
def query_reports_for_project(project_id, sample_ids):
    records = []
    for sample_id in sample_ids:
        pattern = rf'.*{re.escape(sample_id)}.*_(SNV|CNV)\.xlsx$'
        try:
            matching_files = dxpy.find_data_objects(
                project=project_id,
                name=pattern,
                name_mode='regexp',
                describe={'fields': {'name': True}}
            )
            for file in matching_files:
                file_name = file['describe']['name']
                records.append({
                    'sample_id': sample_id,
                    'project_id': project_id,
                    'file_name': file_name
                })
        except Exception as e:
            print(f"Error fetching files for sample {sample_id} in project {project_id}: {e}")
    return records

# Main function to fetch all reports
def fetch_all_reports(df, chunk_size=100, max_workers=16):
    sample_ids = df['sample_id'].tolist()
    project_info = get_matching_projects()
    project_dict = dict(project_info)

    # Chunk sample IDs to manage search load
    chunks = [sample_ids[i:i + chunk_size] for i in range(0, len(sample_ids), chunk_size)]

    all_records = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for project_id in project_dict.keys():
            for chunk in chunks:
                futures.append(executor.submit(query_reports_for_project, project_id, chunk))

        for future in as_completed(futures):
            try:
                records = future.result()
                all_records.extend(records)
            except Exception as e:
                print(f"Error in future: {e}")

    # Add project names to records
    for record in all_records:
        record['project_name'] = project_dict.get(record['project_id'], 'Unknown')

    # Create a DataFrame from the records
    records_df = pd.DataFrame(all_records)

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

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--clarity_extract",
        type=str,
        help=("Path to Clarity extract file"),
    )

    return parser.parse_args()

def open_files(clarity):
    '''
    Open files and read in file contents to DataFrames
    '''
    with open(clarity) as f:
        clarity_df = pd.read_csv(f, delimiter=',')

    return clarity_df

def create_path(filename, assay, run):
    '''
    Create a path to a specific file on clingen
    Inputs:
        assay (str): either CEN or WES
        filename (str): filename for the xlsx report
        run (str): sequencing run name
    Outputs:
        path (str): path to the given filename on clingen
    '''
    base_path = "/appdata/clingen/cg/Regional Genetics Laboratories/Molecular Genetics/Data archive/Sequencing HT/"
    if assay == "CEN":
        path = base_path + f"{assay}/Run folders/{run}_{assay}/{filename}"
    elif assay == "WES":
        path = base_path + f"{assay}/{run}_TWE/{filename}"
    return path

def find_file_name(search_query):
    '''
    Find the file name for a given search query on DNAnexus
    Inputs:
        search_query (str): search query for file name when searching DNAnexus
    Outputs:
        filename (str): a file name, or None if no or multiple matches found
    '''
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
    '''
    Script entry point
    '''
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
    report_df['full_sample_id'] = report_df['sample_id'].str.split('-').str[0]
    # Create output files
    report_df.to_csv("all_data.csv")

if __name__ == "__main__":
    main()
