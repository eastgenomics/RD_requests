import dxpy
import pandas as pd
import subprocess
import argparse
import re
import io as io

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
        A list of excel file objects matching the name.
    """
    excels = list(
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

    return excels

# convert file information into dataframe
def convert_to_df(
    excel_list: list
) -> pd.DataFrame:
    """
    Convert a list of excel file metadata to a pandas DataFrame, exclude rows
    from DF using exclude_projects and exclude_samples lists.

    Parameters
    ----------
    excel_list : list
        list of excel file metadata dicts.
    exclude_projects : list
        List of project IDs to exclude from the DataFrame.
    exclude_samples : list
        List of sample names or patterns to exclude from the DataFrame.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the excel file metadata.
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
        for item in excel_list
    ]
    df = pd.DataFrame(flat_list)

    # Add useful fields for later
    df["sample"] = df["name"].str.split("-").str[0:2].str.join("-")
    df["project_file"] = df["project_id"] + ":" + df["file_id"]
    df["assay"] = df["project_name"].str.split("_").str[-1]

    return df



def main() -> None:
    args = parse_arguments()

    found_projects = find_projects(
        name=args.project_search_term,
        created_after=args.after_date,
        created_before=args.before_date
    )
    print("number of projects in time period: ", len(found_projects))
    text_files = []
    for project in found_projects:

        files_in_project = find_files_in_project(
            project, args.file_search_term
        )
        for file in files_in_project:
            file["project_name"] =  project["describe"]["name"]
        text_files.extend(files_in_project)

    all_files = convert_to_df(text_files)
    
    all_files.to_csv(
        f"file_status.tsv",
        sep="\t",
        index=False,
    )

    all_read_counts = pd.DataFrame() 
    runs_used = 0
    for file_id in all_files["file_id"]:
        print(file_id)
        try :
            cmd = (
            f"dx cat {file_id} "
            )
            output = subprocess.run(cmd, shell=True,
                             capture_output=True, check=False)
            df = pd.read_csv(io.BytesIO(output.stdout), sep="\t")
            # print(df[["mapped_passed", "total_passed"]])
            to_add = df[["mapped_passed", "total_passed"]]
            all_read_counts = pd.concat([all_read_counts, to_add], axis=0, ignore_index= True)
            runs_used+=1
        except pd.errors.EmptyDataError as error:
            print(f"Error reading file {file_id}: {error}")
            print(f"Check archival status of {file_id}")
            pass
    # print(all_read_counts)
    mean_mapped_passed = round(all_read_counts["mapped_passed"].mean()/1000000, 2)
    mean_total_passed = round(all_read_counts["total_passed"].mean()/1000000,2)
    print(f"runs used in calculation: {runs_used}")
    print(f"Mean N of Mapped passed reads per sample: {mean_mapped_passed} Mb")
    print(f"Mean N of total passed reads per sample: {mean_total_passed} Mb")

if __name__ == "__main__":
    main()


