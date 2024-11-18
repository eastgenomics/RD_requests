import argparse
import dxpy
import pandas as pd
import re
import sys

from collections import Counter
from time import sleep


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description=(
            "Information required to find VCFs in DNAnexus for creation of a "
            "pop AF VCF for a specific assay"
        )
    )

    parser.add_argument(
        "-a",
        "--assay",
        type=str,
        required=True,
        help="DNAnexus project suffix to search for",
    )

    parser.add_argument(
        "-s",
        "--start",
        type=str,
        default=None,
        help="Start date to search for corresponding b37 projects",
    )

    parser.add_argument(
        "-e",
        "--end",
        type=str,
        default=None,
        help="End date to search for corresponding b37 projects",
    )

    parser.add_argument(
        "-o",
        "--outfile_prefix",
        type=str,
        required=True,
        help="Prefix to name the output file",
    )

    return parser.parse_args()


def find_projects(project_name, start=None, end=None):
    """
    Find DNANexus projects by name

    Parameters
    ----------
    project_name : str
        project name
    start : str (optional)
        start date to look for projects from
    end: str (optional)
        end date to look for projects until

    Returns
    -------
    projects : list
        list of DNAnexus projects
    """
    projects = list(
        dxpy.find_projects(
            name=project_name,
            created_before=end,
            created_after=start,
            name_mode="glob",
            describe=True
        )
    )

    return projects


def find_data(file_name, project_id):
    """
    Find files in DNAnexus project by name

    Parameters
    ----------
    file_name : str
        file name to search for
    project_id : str
        DX project ID

    Returns
    -------
    files : list
        list of files found in project
    """
    files = list(
        dxpy.find_data_objects(
            name=file_name,
            name_mode="glob",
            project=project_id,
            describe=True
        )
    )
    return files


def read_in_qc_file_to_df(qc_file, b37_proj_id):
    """
    Read in QC status file to a pandas dataframe

    Parameters
    ----------
    qc_file : dict
        dict with info about a qc status file in DNAnexus
    b37_proj : dict
        dict with info about the b37 project it's present in

    Returns
    -------
    qc_df : pd.DataFrame
        the QC status file read in as a dataframe
    """
    file = dxpy.open_dxfile(
        qc_file["id"], project=b37_proj_id, mode='rb'
    )

    file_contents = file.read()
    params = {
        "engine": "openpyxl",
        "usecols": range(8),
        "names": [
            "Sample",
            "M Reads Mapped",
            "Contamination (S)",
            "% Target Bases 20X",
            "% Aligned",
            "Insert Size",
            "QC_status",
            "Reason",
        ],
    }
    try:
        qc_df = pd.read_excel(file_contents, **params)
    # One QC status file weirdly has two sheets so read in from the second
    except ValueError:
        qc_df = pd.read_excel(
            file_contents, sheet_name="Sheet2", **params
        )

    return qc_df


def get_qc_files(b38_projects, start=None, end=None):
    """
    Find QC status files and read in to create list of data frames

    Parameters
    ----------
    b38_projects : list
        list of dicts, each representing a DNAnexus project
    start : str (optional)
        start date to look for projects from
    end: str (optional)
        end date to look for projects until


    Returns
    -------
    all_qc_files : list
        list of dicts, each representing a QC status file in DX
    missing_projects : list
        list of projects that are missing QC files
    b38_project_subset : list
        list of b38 projects that have corresponding b37 projects
    """
    b37_projects = []
    b37_project_dict = {}
    all_qc_files = []
    missing_projects = []
    b38_project_subset = []
    for b38_proj in b38_projects:
        folder_002 = (
            b38_proj["describe"]["name"]
            .rsplit("_", maxsplit=1)[0]
            .split("_", maxsplit=1)[1]
        )
        run_name = f"002_{folder_002}*"
        b37_project = find_projects(run_name, start, end)

        for i, b37_proj in enumerate(b37_project):
            if i > 0:
                print("More than one b37 project found")
                sys.exit()  # checks if multiple b37 projects are found.
            qc_files = find_data("*QC*.xlsx", b37_proj["describe"]["id"])
            print(
                f"Found {len(qc_files)} QC files in {b37_proj['id']} - "
                f"{b37_proj['describe']['name']}"
            )
            if len(qc_files) > 1:
                print(
                    f"\n{len(qc_files)} QC files found in {b37_proj['id']}. "
                    "Taking latest QC status file"
                )
                qc_file = max(
                    qc_files, key=lambda x: x['describe']['created']
                )
                # Append the b38 project to the subset list of projects
                # that have corresponding b37 projects and QC file.
                b38_project_subset.append(b38_proj)
                # Append the b37 project to the list of b37 projects
                # that have QC files.
                b37_projects.append(b37_proj)
                b37_project_dict[b37_proj["describe"]["name"]] = b37_proj["id"]
            elif len(qc_files) == 1:
                qc_file = qc_files[0]
                # Append the b38 project to the subset list of projects
                # that have corresponding b37 projects and QC file.
                b38_project_subset.append(b38_proj)
                # Append the b37 project to the list of b37 projects
                # that have QC files.
                b37_projects.append(b37_proj)
                b37_project_dict[b37_proj["describe"]["name"]] = b37_proj["id"]
            else:
                print(
                    f"No QC files found for this project: "
                    f"{b37_proj['id']} - {b37_proj['describe']['name']}"
                )
                missing_project_info = {
                    "b37_project_name": b37_proj["describe"]["name"],
                    "b37_project_id": b37_proj["id"],
                    "b38_project_name": b38_proj["describe"]["name"],
                    "b38_project_id": b38_proj["id"]
                }
                missing_projects.append(missing_project_info)
                continue
            all_qc_files.append(qc_file)
    print(len(all_qc_files), "QC files found in total")
    print(len(b37_projects), "b37 projects found in total")

    return all_qc_files, missing_projects, b38_project_subset


def unarchive_qc_status_files(all_qc_files):
    """
    Unarchive any QC status files that are not live

    Parameters
    ----------
    all_qc_files : list
        list of dicts, each representing a QC status file in DX
    """
    non_live_files = [
        qc_file for qc_file in all_qc_files
        if qc_file['describe']['archivalState'] != 'live'
    ]

    if non_live_files:
        for non_live_file in non_live_files:
            print(
                "Requesting unarchiving for QC status file "
                f"{non_live_file['id']}"
            )
            file_object = dxpy.DXFile(
                non_live_file["id"], project=non_live_file["project"]
            )
            file_object.unarchive()
            sleep(5)
        print(
            "Exiting now. Please re-run once QC status files are unarchived"
        )
        sys.exit()


def read_in_qc_files_to_df(all_qc_files):
    """
    Read in all QC status files to a single dataframe

    Parameters
    ----------
    all_qc_files : list
        list of dicts, each representing a QC status file in DX

    Returns
    -------
    merged_qc_df : pd.DataFrame
        a single pandas df with all QC status files merged
    """
    qc_file_dfs = []
    for qc_file in all_qc_files:
        qc_df = read_in_qc_file_to_df(qc_file, qc_file["project"])
        qc_file_dfs.append(qc_df)

    print(f"Read in {len(qc_file_dfs)} QC status files")
    merged_qc_df = pd.concat(qc_file_dfs)

    return merged_qc_df


def get_failed_samples(qc_status_df):
    """
    Get any failed samples

    Parameters
    ----------
    qc_status_df : df
        pandas df which is a merge of all QC status files

    Returns
    -------
    fail_sample_names : list
        list of names of samples which have failed
    """
    df_fail = qc_status_df.loc[
        qc_status_df['QC_status'].str.upper() == 'FAIL'
    ]
    fail_samples = list(df_fail['Sample'])
    fail_sample_names = list(set([
        sample.split("-")[0] + "-" + sample.split("-")[1]
        for sample in fail_samples
    ]))

    return fail_sample_names


def get_sample_types(projects):
    """
    Get validation and non-validation samples

    Parameters
    ----------
    projects : list
        list of dicts, each representing a DNAnexus project

    Returns
    -------
    all_non_validation_samples : list
        list of dicts, each with info about a non-validation sample
    all_validation_samples : list
        list of dicts, each with info about a validation sample
    """
    all_validation_samples = []
    all_non_validation_samples = []

    for project in projects:
        vcf_files = find_data(
            "*_markdup_recalibrated_Haplotyper.vcf.gz",
            project["describe"]["id"]
        )
        non_validation_samples_in_run = []

        for vcf in vcf_files:
            instrument_id = vcf["describe"]["name"].split("-")[0]
            sample_id = vcf["describe"]["name"].split("-")[1]
            file_id = vcf["describe"]["id"]

            if (
                re.match(r"^\d{9}$", instrument_id, re.IGNORECASE)
                or re.match(r"^[X]\d{6}$", instrument_id, re.IGNORECASE)
            ) and (
                re.match(r"^[G][M]\d{6,7}$", sample_id, re.IGNORECASE)
                or re.match(r"^\d{5}[R]\d{4}$", sample_id, re.IGNORECASE)
            ):
                all_non_validation_samples.append(
                    {
                        "sample": instrument_id + "-" + sample_id,
                        "project": project["describe"]["id"],
                        "file_id": file_id
                    }
                )
                non_validation_samples_in_run.append(
                    instrument_id + "-" + sample_id
                )

            else:
                all_validation_samples.append(
                    {
                        "sample": instrument_id + "-" + sample_id,
                        "project": project["describe"]["id"],
                        "file_id": file_id
                    }
                )

        if (
            len(non_validation_samples_in_run)
            != len(list(set(non_validation_samples_in_run)))
        ):
            print("Sample duplication in the same run", project['id'])

    return all_non_validation_samples, all_validation_samples


def main():
    args = parse_args()

    b38_projects = find_projects(args.assay)
    projects_to_print = '\n\t'.join([
        f"{x['describe']['name']} - {x['id']}" for x in b38_projects
    ])
    print(f"\n{len(b38_projects)} projects found:\n\t{projects_to_print}")

    # Get QC status files from b37 projects and read them in
    all_qc_files, missing_projects, b38_project_subset = get_qc_files(
        b38_projects, args.start, args.end
    )

    unarchive_qc_status_files(all_qc_files)
    merged_qc_file_df = read_in_qc_files_to_df(all_qc_files)

    # Get any failed samples from QC status reports
    fail_sample_names = get_failed_samples(merged_qc_file_df)
    print("\nFailed samples:")
    print("\n".join(sample for sample in fail_sample_names))

    # Get validation and duplicated samples
    non_validation_samples, validation_samples = get_sample_types(
        b38_project_subset
        )

    # Check duplicated samples from all b38 folders
    sample_names = [item['sample'] for item in non_validation_samples]
    duplicated_samples = [
        item for item, count in Counter(sample_names).items()
        if count > 1
    ]
    print("\nDuplicated_samples:")
    print("\n".join(sample for sample in duplicated_samples))

    # Create a list of all missing projects
    if missing_projects:
        missing_projects_filename = f"{args.outfile_prefix}_projects_missing_QC.csv"
        print(
            f"Outputting projects missing QC files: "
            f"{missing_projects_filename}"
        )
        # convert list to pd.DataFrame and then CSV output.
        df_missing_projects = pd.DataFrame(missing_projects)
        df_missing_projects.to_csv(missing_projects_filename, index=False)

    # Create dfs
    df_validation_samples = pd.DataFrame(validation_samples)
    df_all_non_validation_samples = pd.DataFrame(non_validation_samples)
    df_validation_samples.to_csv(
        f"{args.outfile_prefix}_validation_samples.csv", index=False
    )

    # Drop the duplicated samples and keep once
    df_non_duplicated = df_all_non_validation_samples.drop_duplicates(
        subset=["sample"], keep="last"
    )
    df_non_duplicated.reset_index(drop=True, inplace=True)
    print(
        f"\n{len(df_all_non_validation_samples)} non-validation samples found"
    )
    print(
        f"{len(df_non_duplicated)} non-duplicated non-validation samples "
        "found"
    )

    # Get list of non-failed non-validation samples to merge
    print("Removing failed samples")
    df_file_to_merge = df_non_duplicated[
        ~df_non_duplicated['sample'].isin(fail_sample_names)
    ]
    df_file_to_merge.to_csv(
        f"{args.outfile_prefix}_files_to_merge.txt", sep="\t", header=False
    )
    print("Number of final VCF files to merge:", len(df_file_to_merge))

    # Simple check we don't have any failed or duplicated samples left
    for fail_sample in fail_sample_names:
        if fail_sample in list(df_file_to_merge["sample"]):
            print(f"Failed file found: {fail_sample}")

    for dup_sample in duplicated_samples:
        if dup_sample in list(df_file_to_merge["sample"]):
            print(f"Duplicated sample found: {dup_sample}")


if __name__ == '__main__':
    main()
