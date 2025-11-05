import dxpy
import pandas as pd
import argparse
import re


def parse_arguments():
    """
    Parse command-line arguments.

    Returns
    -------
    args: argparse.Namespace
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Search for VCF files in DNAnexus projects."
    )
    parser.add_argument(
        "--before_date", type=str, required=False,
        help="Date to search before (YYYY-MM-DD)."
    )
    parser.add_argument(
        "--after_date", type=str, required=True,
        help="Date to search after (YYYY-MM-DD)."
    )
    parser.add_argument(
        "--project_search_term", type=str, required=True,
        help="Search term for project names."
    )
    parser.add_argument(
        "--file_search_term", type=str, required=True,
        help="Regex pattern for file names."
    )
    parser.add_argument(
        "--unarchive", action="store_true", default=False,
        help="Whether to unarchive files."
    )
    parser.add_argument(
        "--output_prefix", type=str, required=False, default="",
        help="Prefix for output file names."
    )
    parser.add_argument(
        "--exclude_projects", nargs='+', type=str, required=False, default=[],
        help="List of projects IDs to exclude if found in search, e.g. validation runs." \
        "example format \"--exclude_projects project-xxxx" \
        "project-xxxx\""
    )
    parser.add_argument(
        "--exclude_samples", nargs='+', type=str, required=False, default=[],
        help="List of samples to exclude if found in search, ie top-up samples. " \
        "This will remove samples from all runs" \
        "example format \" --exclude_samples 12345K0067 12345K0089\""
    )
    parser.add_argument(
        "--exclude_sample_on_run", nargs='+', type=str, required=False, default=[],
        help="List of samples in specific runs to exclude if found in search, ie failed samples. " \
        "This will remove samples from specific runs" \
        "example format \" --exclude_sample_on_run 002_240221_A01303_0346_AHYJC3DRX3_38_CEN:23341R0046\""
    )
    return parser.parse_args()


def find_projects(name: str,
                  name_mode: str = "regexp",
                  created_after: str = None,
                  created_before: str = None
) -> list:
    """
    Find projects by name and created date.

    Parameters
    ----------
    name : str
        The name or glob pattern to search for.
    name_mode : str
        The mode for name matching.
    created_after : str, optional
        The date to search for projects created after (YYYY-MM-DD).
    created_before : str, optional
        The date to search for projects created before (YYYY-MM-DD).

    Returns
    -------
    list
        A list of project objects (dicts) matching the name and date criteria.
    """
    projects = list(
        dxpy.find_projects(
            name=name,
            name_mode=name_mode,
            created_after=created_after,
            created_before=created_before,
            describe={"fields": {"name": True}},
        )
    )

    if not projects:
        raise RuntimeError("No projects found matching criteria")

    return projects


def find_files_in_project(project: dict, name: str,
                          name_mode: str = "regexp") -> list:
    """
    Find files in a project by name.

    Parameters
    ----------
    project : dict
        The project object.
    name : str
        The name to search for.
    name_mode : str
        The mode for name matching.

    Returns
    -------
    list
        A list of VCF file objects matching the name.
    """
    vcfs = list(
        dxpy.find_data_objects(
            project=project["id"],
            name=name,
            name_mode=name_mode,
            classname="file",
            describe={
                "fields": {
                    "archivalState": True,
                    "name": True,
                    "created": True,
                }
            },
        )
    )

    return vcfs


def bulk_unarchive_per_project(df: pd.DataFrame):
    """
    Bulk unarchive files per project.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing file information.
    """
    for project_id, group in df.groupby("project_id"):
        file_ids = group["file_id"].tolist()
        print(f"Unarchiving {len(file_ids)} files in project: {project_id}")

        try:
            dxpy.api.project_unarchive(
                project_id,
                input_params={"files": file_ids},
            )
        except Exception as error:
            print(f"Error unarchiving files for {project_id}: {error}")
            raise RuntimeError("Error unarchiving files") from error


def convert_to_df(
    vcf_list: list,
    exclude_projects: list,
    exclude_samples: list,
    exclude_sample_on_run: list
) -> pd.DataFrame:
    """
    Convert a list of VCF file metadata to a pandas DataFrame, exclude rows
    from DF using exclude_projects and exclude_samples lists.

    Parameters
    ----------
    vcf_list : list
        list of VCF file metadata dicts.
    exclude_projects : list
        List of project IDs to exclude from the DataFrame.
    exclude_samples : list
        List of sample names or patterns to exclude from the DataFrame.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the VCF file metadata.
    """
    flat_list = [
        {
            "file_id": item["id"],
            "project_id": item["project"],
            "project_name": item["project_name"],
            "name": item["describe"]["name"],
            "created": item["describe"]["created"],
            "archive_state": item["describe"]["archivalState"],
        }
        for item in vcf_list
    ]
    df = pd.DataFrame(flat_list)

    # Add useful fields for later
    df["sample"] = df["name"].str.split("-").str[0:2].str.join("-")
    df["project_file"] = df["project_id"] + ":" + df["file_id"]
    df["assay"] = df["project_name"].str.split("_").str[-1]

    if exclude_projects:
        print(f"Excluding the following projects: {exclude_projects}")
        # remove unwanted samples and projects
        df = df[~df["project_id"].isin(exclude_projects)]

    if exclude_samples:
        print(f"Excluding the following samples from all runs: {exclude_samples}")
        # Only know part of top up samples so need to set up a regex
        escaped_samples = [re.escape(sample) for sample in exclude_samples]
        pattern = "|".join(escaped_samples)
        df = df[~df["sample"].str.contains(pattern, regex=True)]

    if exclude_sample_on_run:
        print(f"Excluding the following samples from specific runs: {exclude_sample_on_run}")
        for pair in exclude_sample_on_run:
            sample_to_remove = pair.split(":")[1]
            run_to_remove_from = pair.split(":")[0]
            df = df[~(df["sample"].str.contains(sample_to_remove) & df["project_name"].str.contains(run_to_remove_from))]
        # remove unwanted samples and projects
    return df


def remove_controls_and_dups(df: pd.DataFrame) -> tuple:
    """
    Remove control samples (containing "Q") and duplicates (by sample)
    from the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame containing VCF file information.

    Returns
    -------
    tuple
        A tuple containing two DataFrames:
        - df_no_control (pd.DataFrame): A DataFrame with control samples removed,
          including additional columns for project_list (list of all projects
          the sample is found in) and project_count (number of projects the
          sample is found in).
        - df_no_dups (pd.DataFrame): A DataFrame with duplicates removed,
          keeping the latest VCF for each sample based on the "created" column.
    """
    # Remove controls and X numbers
    valid_pattern = r"^\d{9}-\d{5}[RSK]\d{4}"
    mask = df["sample"].str.contains(valid_pattern, regex=True) & ~df[
        "sample"
    ].str.contains(r"0--|NA|Oncospan|ctrl|Q", case=False, regex=True)
    df_no_control = df[mask].copy()

    # Add column containing a list of all projects the sample is found in
    project_lists = (
        df_no_control.groupby("sample")["project_name"]
        .apply(lambda x: list(x.unique()))
        .reset_index(name="project_list")
    )

    # Merge the project_list back into df_no_control
    df_no_control = df_no_control.merge(project_lists, on="sample", how="left")

    # Add a column with the count of projects each sample is found in
    df_no_control["project_count"] = df_no_control["project_list"].apply(len)

    # Remove dups, keeping the latest VCF for each sample
    df_no_control = df_no_control.sort_values(by="created", ascending=False)
    df_no_dups_no_control = df_no_control.drop_duplicates(subset="sample", keep="first")

    return df_no_control, df_no_dups_no_control


def main() -> None:
    args = parse_arguments()

    vcf_projects = find_projects(
        name=args.project_search_term,
        created_after=args.after_date,
        created_before=args.before_date
    )
    vcfs = []
    for project in vcf_projects:
        vcfs_in_project = find_files_in_project(
            project, args.file_search_term
        )
        for vcf in vcfs_in_project:
            vcf["project_name"] = project["describe"]["name"]
        vcfs.extend(vcfs_in_project)

    all_vcfs = convert_to_df(vcfs, args.exclude_projects, args.exclude_samples, args.exclude_sample_on_run)
    df_no_control, df_no_dups_no_control = remove_controls_and_dups(all_vcfs)

    # Write out summary TSV of files found
    all_vcfs.to_csv(
        f"{args.output_prefix}VCFs_all.tsv",
        sep="\t",
        index=False,
    )

    df_no_control.to_csv(
        f"{args.output_prefix}VCFs_no_controls.tsv",
        sep="\t",
        index=False,
    )

    df_no_dups_no_control.to_csv(
        f"{args.output_prefix}VCFs_no_dups_no_control.tsv",
        sep="\t",
        index=False,
    )

    non_live = df_no_dups_no_control[df_no_dups_no_control["archive_state"] != "live"]

    # Write out VCFs to unarchive and call unarchiving
    if not non_live.empty:
        non_live["project_file"].to_csv(
            f"{args.output_prefix}VCFs_to_unarchive.txt",
            index=False,
            header=False,
        )
        if args.unarchive:
            print("Unarchiving files...")
            bulk_unarchive_per_project(non_live)

    # Write out file IDs for each assay
    for assay in df_no_dups_no_control["assay"].unique():
        assay_vcfs = df_no_dups_no_control[df_no_dups_no_control["assay"] == assay]

        assay_vcfs["project_file"].to_csv(
            f"{args.output_prefix}{assay}_VCFs_ids.txt",
            index=False,
            header=False,
        )


if __name__ == "__main__":
    main()
