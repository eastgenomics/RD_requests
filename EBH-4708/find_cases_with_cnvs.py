import argparse
import dxpy
import pandas as pd
import io


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Information required to find workbooks"
    )

    parser.add_argument(
        "-n",
        "--project_name",
        required=True,
        type=str,
        help="Name regex to filter 002 projects with, e.g. '^002.*_CEN$'",
    )

    parser.add_argument(
        "-r",
        "--r_code",
        required=True,
        type=str,
        help="R code to use to check for CNV reports, e.g. 'R228.1'",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=str,
        help="Name of the output Excel file",
    )

    return parser.parse_args()


def find_projects(name: str, date=None) -> list:
    """
    Find DNAnexus projects matching the name and date provided.

    Parameters
    ----------
    name: str
        name used to search projects
    date : str (optional)
        date in the format 'YYYY-MM-DD'

    Returns
    -------
    list[dict]
        projects found
    """
    projects = list(
        dxpy.find_projects(
            name=name,
            name_mode="regexp",
            created_after=date,
            describe={"fields": {"name": True}},
        )
    )

    return projects


def find_cnv_reports(project: list, test_code: str) -> list:
    """
    Find specific files in the given projects.

    Parameters
    ----------
    projects : list
        list of projects to search in
    test_code : str
        the test code to look for

    Returns
    -------
    list[dict]
        list of files found in the projects.
    """
    reports = list(
        dxpy.find_data_objects(
            project=project["id"],
            name=f".*_{test_code}\_CNV_\d+\.xlsx$",
            name_mode="regexp",
            describe={
                "fields": {
                    "name": True,
                    "archivalState": True,
                    "created": True,
                }
            },
        )
    )

    # Add in project name so we can tell if GRCh37 or GRCh38 later
    for report in reports:
        report["project_name"] = project["describe"]["name"]

    return reports


def get_variant_details(list_of_files: list) -> list:
    """
    Add in number of CNVs detected per workbook by querying the file
    details metadata.

    Parameters
    ----------
    list_of_files : list[dict]
        list of files to query

    Returns
    -------
    list[dict]
        list of files, with new key 'variants'
    """
    for report in list_of_files:
        details = dxpy.DXFile(
            dxid=report["id"], project=report["project"]
        ).get_details()
        report["variants"] = details.get("variants")

    return list_of_files


def unarchive_files(list_of_files) -> None:
    """
    Unarchive any files not in the live state.

    Parameters
    ----------
    list_of_files : list
        list of files with the archivalState key
    """
    for report in list_of_files:
        if report["describe"]["archivalState"] != "live":
            print(
                f"Trying to access {report['id']} but it is not in the live"
                " state. Now requesting unarchiving"
            )
            file_object = dxpy.DXFile(report["id"], project=report["project"])
            file_object.unarchive()


def format_cases(list_of_files: list) -> pd.DataFrame:
    """
    Format the CNV reports found into a pandas dataframe.

    Parameters
    ----------
    list_of_files : list[dict]
        list of dicts, each representing a file

    Returns
    -------
    pd.DataFrame
        dataframe with one row per file found
    """
    df = pd.json_normalize(list_of_files, sep=".")

    # Add in genome assembly from project name
    df["assembly"] = (
        df["project_name"]
        .str.extract(r"_(\d{2})_", expand=False)
        .fillna("37")
        .astype(int)
    )

    # Rename columns from the describe nested dict
    df.rename(
        columns={
            "describe.id": "file_id",
            "describe.name": "file_name",
            "describe.archivalState": "archive_state",
            "describe.created": "created",
        },
        inplace=True,
    )

    # Get specimen ID from file name
    df["sample_id"] = df["file_name"].str.extract(
        r"-([^-\s]*R[^-\s]*)-", expand=False
    )
    # If nothing,  get the X number from the first part
    df["sample_id"] = df["sample_id"].fillna(
        df["file_name"].str.split("-", n=1).str[0]
    )

    return df


def get_cases_to_read_in(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter for cases which are live and we know have variants or
    can't get that info from metadata so must parse the CNV workbook

    Parameters
    ----------
    df : pd.DataFrame
        dataframe of cases

    Returns
    -------
    pd.DataFrame
        dataframe of cases we are going to read in
    """
    df = df[
        ((df["variants"] > 0) | (df["variants"].isna()))
        & (df["archive_state"] == "live")
    ]

    return df


def read_in_report(report: dict) -> pd.DataFrame:
    """
    Read in the CNV excel workbook variants sheet to a pandas dataframe.

    Parameters
    ----------
    report : dict
        dict with info about a CNV report

    Returns
    -------
    pd.DataFrame
        dataframe of the variants sheet info
    """
    try:
        dxfile = dxpy.open_dxfile(
            report["id"], project=report["project"], mode="rb"
        )
        df = pd.read_excel(
            io.BytesIO(dxfile.read()),
            sheet_name="variants",
            engine="openpyxl",
        )

        df["file_id"] = report["file_id"]
        df["project"] = report["project"]
    except ValueError:
        print(f"Cant read in {report['file_id']}")
        return None

    return df


def get_variant_counts_from_excel(
    cases_to_read_in: pd.DataFrame,
) -> pd.DataFrame:
    """
    Get the counts from the Excel workbook for any cases with no
    details metadata

    Parameters
    ----------
    cases_to_read_in : pd.DataFrame
        rows of cases which we need to read in

    Returns
    -------
    pd.DataFrame
        info about the counts of variants for each case
    """
    variant_counts_list = []
    for _, row in cases_to_read_in.iterrows():
        df = read_in_report(row)
        if df is not None:
            count = len(df)
            variant_counts_list.append(
                {
                    "file_id": row["file_id"],
                    "project": row["project"],
                    "variant_count_from_excel": count,
                }
            )
    variant_counts_df = pd.DataFrame(variant_counts_list)

    return variant_counts_df


def merge_variant_counts(
    file_df: pd.DataFrame, variant_counts_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge variant counts into main file dataframe and fill missing variants.

    Parameters
    ----------
    file_df : pd.DataFrame
        info for each case
    variant_counts_df : pd.DataFrame
        info of variant counts from reading in Excel files

    Returns
    -------
    pd.DataFrame
        dataframe with variant info from the Excel workbook
    """
    merged = pd.merge(
        file_df,
        variant_counts_df,
        on=["file_id", "project"],
        how="left",
    )
    merged["variants"] = (
        merged["variants"]
        .fillna(merged["variant_count_from_excel"])
        .astype(int)
    )

    return merged


def remove_duplicates(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove report duplicates for each sample for the same genome assembly
    with the same number of variants detected.

    Parameters
    ----------
    merged_df : pd.DataFrame
        dataframe with duplicate reports per sample

    Returns
    -------
    pd.DataFrame
        dataframe with duplicates removed
    """
    merged_clean = (
        merged_df.sort_values(by="created", ascending=False)
        .drop_duplicates(
            subset=["sample_id", "assembly", "variants"], keep="first"
        )
        .reset_index(drop=True)
    )

    return merged_clean


def main():
    args = parse_args()
    projects = find_projects(args.project_name)
    files_found = []
    for project in projects:
        files_found.extend(find_cnv_reports(project, args.r_code))
    files_found = get_variant_details(files_found)
    unarchive_files(files_found)
    file_df = format_cases(files_found)
    cases_to_read_in = get_cases_to_read_in(file_df)
    variant_counts_df = get_variant_counts_from_excel(cases_to_read_in)
    merged = merge_variant_counts(file_df, variant_counts_df)

    merged_clean = remove_duplicates(merged)
    merged_clean = merged_clean.sort_values(by=["assembly", "sample_id"])
    merged_clean.drop(
        columns=["variant_count_from_excel", "created", "archive_state"],
        inplace=True,
    )
    merged_clean.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
