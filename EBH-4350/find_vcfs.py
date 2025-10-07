import dxpy
import pandas as pd
import argparse


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


def convert_to_df(vcf_list: list):
    """
    Convert a list of VCF file metadata to a pandas DataFrame.

    Parameters
    ----------
    vcf_list : list
        list of VCF file metadata dicts.

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

    return df


def remove_controls_and_dups(df: pd.DataFrame):
    """
    Remove control samples (containing "Q") and duplicates (by sample)
    from the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame containing VCF file information.

    Returns
    -------
    pd.DataFrame
        A DataFrame with control samples and duplicates removed.
    """
    # Remove controls and X numbers
    valid_pattern = r"^\d{9}-\d{5}R\d{4}"
    mask = df["sample"].str.contains(valid_pattern, regex=True) & ~df[
        "sample"
    ].str.contains(r"0--|NA|Oncospan|ctrl|Q", case=False, regex=True)
    df_no_control = df[mask]

    # Remove dups, keeping the latest VCF for each sample
    df_no_control = df_no_control.sort_values(by="created", ascending=False)
    df_no_dups = df_no_control.drop_duplicates(subset="sample", keep="first")

    return df_no_dups


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

    all_vcfs = convert_to_df(vcfs)
    all_vcf_no_dups = remove_controls_and_dups(all_vcfs)
    all_vcf_no_dups.to_csv(
        "GRCh38_CEN_TWE_Sentieon_VCFs_no_controls_or_dups.tsv",
        sep="\t",
        index=False,
    )

    non_live = all_vcf_no_dups[all_vcf_no_dups["archive_state"] != "live"]

    # Write out summary CSV of files found
    all_vcf_no_dups.to_csv(
        f"{args.output_prefix}all_VCFs_summary.csv",
        index=False,
    )

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
    for assay in all_vcf_no_dups["assay"].unique():
        assay_vcfs = all_vcf_no_dups[all_vcf_no_dups["assay"] == assay]

        assay_vcfs["project_file"].to_csv(
            f"{args.output_prefix}{assay}_VCFs_ids.txt",
            index=False,
            header=False,
        )


if __name__ == "__main__":
    main()
