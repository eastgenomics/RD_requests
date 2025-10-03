import dxpy
import pandas as pd


def find_projects(name: str, name_mode: str) -> list:
    """
    Find projects by name.

    Parameters
    ----------
    name : str
        The name or glob pattern to search for.
    name_mode : str
        The mode for name matching.

    Returns
    -------
    list
        A list of project objects (dicts) matching the name.
    """
    projects = list(
        dxpy.find_projects(
            name=name,
            name_mode=name_mode,
            describe={"fields": {"name": True}},
        )
    )
    return projects


def find_files_in_project(project: dict, name: str, name_mode: str) -> list:
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
    df_no_control = df[~df["sample"].str.contains("Q")]
    # Remove dups, keeping the latest VCF for each sample
    df_no_control = df_no_control.sort_values(by="created", ascending=False)
    df_no_dups = df_no_control.drop_duplicates(subset="sample", keep="first")

    return df_no_dups


def main():
    vcf_projects = find_projects("(^002_.*_38_(CEN|TWE)$)", name_mode="regexp")

    vcfs = []
    for project in vcf_projects:
        vcfs_in_project = find_files_in_project(
            project, "*recalibrated_Haplotyper.vcf.gz", "glob"
        )
        for vcf in vcfs_in_project:
            vcf["project_name"] = project["describe"]["name"]
        vcfs.extend(vcfs_in_project)

    all_vcfs = convert_to_df(vcfs)
    all_vcf_no_dups = remove_controls_and_dups(all_vcfs)

    non_live = all_vcf_no_dups[all_vcf_no_dups["archive_state"] != "live"]

    # Write out VCFs to unarchive and call unarchiving
    if not non_live.empty:
        non_live["project_file"].to_csv(
            "GRCh38_CEN_WES_Sentieon_VCFs_to_unarchive.txt",
            index=False,
            header=False,
        )
        bulk_unarchive_per_project(non_live)

    # Write out file IDs for each assay
    for assay in ["CEN", "TWE"]:
        assay_vcfs = all_vcf_no_dups[all_vcf_no_dups["assay"] == assay]
        assay_vcfs["project_file"].to_csv(
            f"GRCh38_{assay}_Sentieon_VCFs_no_controls_or_dups.txt",
            index=False,
            header=False,
        )


if __name__ == "__main__":
    main()
