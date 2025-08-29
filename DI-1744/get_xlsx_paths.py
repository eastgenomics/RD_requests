"""
This script retrieves and processes data from DNAnexus and Clarity.
"""

import pandas as pd
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

import argparse
import dxpy


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
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="all_data.csv",
        help=("Output file name for the processed data"),
    )
    args = parser.parse_args()

    if not args.clarity_extract:
        parser.error("The --clarity_extract argument is required.")

    return args


def open_files(clarity):
    """
    Open files and read in file contents to DataFrames
    """
    with open(clarity) as f:
        clarity_df = pd.read_csv(f, delimiter=",")

    return clarity_df


def get_matching_projects(assay):
    """
    Retrieve project IDs matching specific name patterns.
    This function searches for projects with names starting with '002' and ending with either 'CEN' or 'TWE'.

    Returns
    -------
    matching_projects_tuple_list: list
        A list of tuples containing project IDs and their names.
        Each tuple is in the format (project_id, project_name).
    """
    pattern = rf"^002.*_{assay}$"
    matching_projects = list(
        dxpy.find_projects(name={"regexp": pattern}, describe=True)
    )
    matching_projects_tuple_list = [
        (proj["id"], proj["describe"]["name"]) for proj in matching_projects
    ]
    print(f"Found {len(matching_projects_tuple_list)} matching projects.")
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
    pattern = rf".*({'|'.join(sample_ids)}).*xlsx"
    try:
        matching_files = dxpy.find_data_objects(
            project=project_id,
            name=pattern,
            name_mode="regexp",
            describe={"fields": {"name": True}},
        )
        # If no files found, return empty records
        if not matching_files:
            print(f"No files found for project {project_id} with pattern {pattern}")
            return records
        for file in matching_files:
            file_name = file["describe"]["name"]
            sample_id = file_name.split("-")[1]
            records.append(
                {
                    "sample_id": sample_id,
                    "project_id": project_id,
                    "file_name": file_name,
                }
            )
    except Exception as e:
        print(f"Error fetching files for sample in project {project_id}: {e}")
    return records


def fetch_all_reports_for_assay(df, assay, chunk_size=100, max_workers=64):
    """
    Fetch all reports for the given DataFrame of sample IDs.
    This function retrieves project IDs matching specific name patterns,
    queries reports for each project, and merges the results with the original DataFrame.
    It also handles errors for samples found in multiple projects.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing sample IDs to search for.
    assay : str
        The assay type to filter projects (e.g., 'CEN' or 'TWE').
    chunk_size : int, optional
        set number of ids to process at a time, by default 100
    max_workers : int, optional
        number of maximum workers to destribute across, by default 16

    Returns
    -------
    final_df : pd.DataFrame
        DataFrame containing the merged results with additional columns.
    """
    sample_ids = df["sample_id"].tolist()
    project_info = get_matching_projects(assay)
    project_dict = dict(project_info)

    # Chunk sample IDs to manage search load
    chunks = [
        sample_ids[i : i + chunk_size] for i in range(0, len(sample_ids), chunk_size)
    ]

    all_records = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for project_id in project_dict.keys():
            for chunk in chunks:
                futures.append(
                    executor.submit(query_reports_for_project, project_id, chunk)
                )

        for future in as_completed(futures):
            try:
                records = future.result()
                all_records.extend(records)
            except Exception as e:
                print(f"Error in future: {e}")

    # Add project names to records
    for record in all_records:
        record["project_name"] = project_dict.get(record["project_id"], "Unknown")

    # Create a DataFrame from the records
    records_df = pd.DataFrame(all_records)
    print(records_df.columns)
    if records_df.empty:
        print("No records found for the given sample IDs.")
        return df
    # Merge with the original df to retain additional columns
    merged_df = pd.merge(df, records_df, on="sample_id", how="left")

    # Identify samples with multiple projects
    project_counts = merged_df.groupby("sample_id")["project_id"].nunique()
    multiple_projects = project_counts[project_counts > 1].index.tolist()
    print(f"Samples with multiple projects: {multiple_projects}")
    # Log errors for samples with multiple projects
    for sample in multiple_projects:
        print(f"Error: Sample {sample} found in multiple projects.")

    # Filter out samples with multiple projects
    final_df = merged_df[~merged_df["sample_id"].isin(multiple_projects)]

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
    # Handle NaN values
    if pd.isna(run) or run == "nan":
        # print(f"Warning: run is NaN or 'nan' for filename {filename}. Returning None.")
        # print(assay)
        return None

    # Convert to string if it's not already
    run = str(run)

    base_path = r"/appdata/clingen/cg/Regional Genetics Laboratories/Molecular Genetics/Data archive/Sequencing HT/"

    run_without_prefix = re.sub(r"^002_", "", run)
    if assay == "CEN":
        path = base_path + rf"{assay}/Run\ folders/{run_without_prefix}/{filename}"
    elif assay == "WES":
        path = base_path + rf"{assay}/{run_without_prefix}/{filename}"
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
    files = list(
        dxpy.find_data_objects(name=search_query, name_mode="glob", describe=True)
    )
    for file in files:
        filenames = []
        filenames.append(file.get("describe").get("name"))
        if len(filenames) != 1:
            print(f"{search_query} returned multiple/no files")
            return None
        else:
            return filenames[0]


def filter_duplicate_files(df):
    """
    Filter duplicate files based on specific rules:
    - For samples with exactly 2 reports: keep files with _CNV_ or _SNV_,
      but filter out files that don't contain _CNV_ or _SNV_ and end in _2
    """
    # Group by sample_id
    grouped = df.groupby("sample_id")
    filtered_rows = []

    for sample_id, group in grouped:
        if len(group) == 2:  # Exactly 2 reports
            # Check each file in the group
            files_to_keep = []
            for _, row in group.iterrows():
                filename = row["file_name"]
                # Keep files that contain _CNV_ or _SNV_
                if "_CNV_" in filename or "_SNV_" in filename:
                    files_to_keep.append(row)
                # Filter out files that don't contain _CNV_ or _SNV_ and end in _2
                elif not (
                    "_CNV_" in filename or "_SNV_" in filename
                ) and filename.endswith("_2.xlsx"):
                    continue  # Skip this file
                else:
                    files_to_keep.append(row)  # Keep other files

            filtered_rows.extend(files_to_keep)
        else:
            # Keep all rows for samples with != 2 reports (will be handled separately)
            filtered_rows.extend([row for _, row in group.iterrows()])

    return pd.DataFrame(filtered_rows).reset_index(drop=True)


def main():
    """
    Script entry point
    """
    # Parse args
    args = parse_args()

    # Read data into dataframes
    clarity_df = open_files(args.clarity_extract)

    # Process data to construct a path for each specimen
    clarity_df = clarity_df.astype(str)
    # create df with column by splitting the Beaker Procedure Name to create a new column for assay
    clarity_df["Assay"] = clarity_df["Beaker Procedure Name"].str.split(" ").str[0]
    clarity_df["sample_id"] = clarity_df["Specimen Identifier"].str.split("-").str[1]
    assays = ["TWE", "CEN"]
    report_df = pd.DataFrame()
    for assay in assays:
        assay_samples = clarity_df[clarity_df["Assay"] == assay]
        assay_df = fetch_all_reports_for_assay(assay_samples, assay)
        report_df = pd.concat([report_df, assay_df], ignore_index=True)

    print(report_df.head())
    print(f"Total reports fetched: {report_df.shape[0]}")
    print(report_df.iloc[0:2, :])

    # Split R codes into a list
    report_df["R_codes"] = report_df["Test Directory Test Code"].str.split("|")
    # remove decimal points from R codes
    report_df["R_codes"] = report_df["R_codes"].apply(
        lambda x: [re.sub(r"\.\d+", "", code) for code in x if code.startswith("R")]
    )

    # Add the path to the processed reports
    report_df["path"] = report_df.apply(
        lambda x: create_path(x["file_name"], x["Assay"], x["project_name"]), axis=1
    )

    # Add specimen ID to the processed report_df
    report_df["specimen_id"] = report_df["file_name"].str.split("-").str[0]
    report_df["full_sample_id"] = (
        report_df["specimen_id"] + "-" + report_df["sample_id"]
    )
    # create report_r_code column from file_name
    report_df["report_r_code"] = report_df["file_name"].str.extract(r"_(R\d+\.\d+)_")[0]

    # Stratify into multiple files for different outcomes
    # Filter out all rows where filename contains _CNV_ or _mosaic_
    report_df = report_df[~report_df["file_name"].str.contains("_CNV_", na=False)]
    report_df = report_df[~report_df["file_name"].str.contains("_mosaic_", na=False)]
    print(
        f"After filtering out CNV and mosaic files: {report_df.shape[0]} rows remaining"
    )

    missing_data = report_df["file_name"].isna() | report_df["R_codes"].isna()
    mising_data_df = report_df[missing_data]
    if missing_data.any():
        print(
            f"Warning: {missing_data.sum()} rows have missing data in 'file_name' or 'R_codes'."
        )
    mising_data_df.to_csv(f"{args.output}_missing_data.csv", index=False)

    # Mask to filter out rows with NaN R codes and empty file names which aren't '' just blank
    report_df = report_df.dropna(subset=["file_name", "R_codes"])
    print(
        f"After filtering out NaN R codes and empty file names: {report_df.shape[0]} rows remaining"
    )

    # Drop duplicates on all columns except certain ones
    cols_to_check = [col for col in report_df.columns if col not in ["R_codes"]]
    report_df = report_df.drop_duplicates(subset=cols_to_check)
    print(f"After dropping duplicates: {report_df.shape[0]} rows remaining")

    # Save rows with multiple reports per assay to a separate file
    multiple_reports = (
        report_df.groupby(["sample_id", "Assay", "report_r_code"])
        .size()
        .reset_index(name="report_count")
    )
    multiple_reports = multiple_reports[multiple_reports["report_count"] > 1]
    if not multiple_reports.empty:
        print(f"Samples with multiple reports per assay: {multiple_reports.shape[0]}")
        # Filter the original report_df to keep only samples with single reports
        multiple_reports_df = pd.merge(
            report_df,
            multiple_reports[["sample_id", "Assay", "report_r_code"]],
            on=["sample_id", "Assay", "report_r_code"],
            how="inner",
        )
        multiple_reports_df.to_csv(f"{args.output}_multiple_reports.csv", index=False)
    else:
        print("No samples with multiple reports per assay found.")

    # Remove rows with multiple reports per assay
    report_df_grouped = (
        report_df.groupby(["sample_id", "Assay", "report_r_code"])
        .size()
        .reset_index(name="report_count")
    )
    single_reports_df = report_df_grouped[report_df_grouped["report_count"] == 1]
    # Filter the original report_df to keep only samples with single reports
    report_df_filtered = pd.merge(
        report_df,
        single_reports_df[["sample_id", "Assay", "report_r_code"]],
        on=["sample_id", "Assay", "report_r_code"],
        how="inner",
    )
    print(f"After multiple reports: {report_df_filtered.shape[0]} rows remaining")

    # Create output files with no index
    report_df_filtered.to_csv(f"{args.output}_to_process.csv", index=False)


if __name__ == "__main__":
    main()
