"""
Created on Mon Dec 16 11:50:04 2024

@author: arun

Code logic:
    1. get a list of bam files from the project of Staging area 055?
    2. get the list of jsons from the same project.
    3. For each JSON, determine the sample name using regex.
    4. For each BAM, determine the sample name using regex
    5. Merge JSON and BAM files into a dataframe where the common samplename

"""
import re
import pandas as pd
import dxpy


# Set Default project ID
PROJECT_ID = "project-GXx8Gg8468bZ072Py1Bb1p0F"

def extract_sample_name(file_name: str) -> str:
    """
    Extracts a sample name from a file name based on a regex.
    """
    pattern = r'(SP[a-zA-Z0-9]+|[gG][mM][0-9]+_?[0-9]+)'
    match = re.search(pattern, file_name)
    return match.group() if match else None


def find_dx_files(project_id: str, extension: str, outfile: str) -> pd.DataFrame:
    """ Finds dx files matching extension in project ID
    """
    print(f"Looking for {extension} files in {project_id}")
    result = list(dxpy.bindings.search.find_data_objects(
        classname="file",
        name=f"*.{extension}",
        name_mode="glob",
        project=project_id,
        describe={"fields": {"name": True,
                             "folder": True}}
    ))

    print(f"Found {len(result)} files with extension: {extension}")

    df = pd.DataFrame([
        {
            "sample_name": extract_sample_name(x["describe"]["name"]),
            f"{extension}_file_name": x["describe"]["name"],
            f"{extension}_file_id": x["id"],
            f"{extension}_folder_path": x["describe"]["folder"]
        } for x in result
    ])
    df.to_csv(outfile, sep="\t", index=False)
    print(f"Files written to {outfile}\n")

    return clean_df(df, extension, outfile)


def clean_df(df: pd.DataFrame, extension: str, outfile: str) -> pd.DataFrame:
    """_
    Log and handle duplicate and missing sample_name(s) in a df.
    """
    sample_col = "sample_name"
    files = f"{extension}_files"
    outfile = outfile.removesuffix(".tsv")
    outfile = f"{outfile}_cleaned.tsv"
    print(f"Cleaning {files} dataframe ...\n")

    rows_with_na = df[df[sample_col].isna()]
    if rows_with_na.empty:
        print(f"All {len(df)} {files} matched a samplename.")
    else:
        print(
            f"Removing the following {len(rows_with_na)} {files} "
            f"with no matched samples:"
        )
        print(rows_with_na)
        df = df.dropna(subset=[sample_col])

    duplicate_rows = df[df.duplicated(subset=sample_col, keep='first')]
    if duplicate_rows.empty:
        print("\nNo duplicate sample name in {files}.")
    else:
        print(
            f"\nRemoving the following {len(duplicate_rows)} rows "
            f"with duplicated sample names in {files}, keeping only "
            f"the first instance."
        )
        print(duplicate_rows)
        df = df.drop_duplicates(subset=sample_col, keep='first')

    df.to_csv(outfile, sep="\t", index=False)
    print(f"Saved cleaned data in {outfile} \n")

    return df


def main():
    """Entry point
    """
    outfile_without_json = "medicover_without_json_file.tsv"
    outfile_with_json = "medicover_with_json_file.tsv"

    bam_df = find_dx_files(PROJECT_ID, "bam", "medicover_bam_files.tsv")
    json_df = find_dx_files(PROJECT_ID, "json", "medicover_json_files.tsv")

    samples_without_json = bam_df[~bam_df.sample_name.isin(json_df.sample_name)]
    samples_without_json.to_csv(outfile_without_json, sep="\t", index=False)

    samples_with_json = pd.merge(bam_df, json_df,
                                 on="sample_name", how="left").dropna()
    samples_with_json.to_csv(outfile_with_json, sep="\t", index=False)

    samples_with_json_in_json_missing = samples_with_json[
        samples_with_json.json_folder_path.str.contains(
            "JSON_MISSING|json_missing", regex=True
        )
    ]

    print(
        f"Found {len(samples_with_json_in_json_missing)} Medicover samples "
         "with a json folder in a JSON_MISSING folder. "
         )

    print(samples_with_json_in_json_missing)

    print("\nSUMMARY...\n")
    print(f"Found {len(bam_df)} unique Medicover samples with sequence data")
    print(
        f"{len(samples_without_json)} "
         "Medicover samples have no JSON result "
        f"files.\n"
        f"{len(samples_with_json)} "
        "Medicover samples have a JSON result file. Where "
        f"{len(samples_with_json_in_json_missing)} "
         "Medicover samples have a JSON result file in a JSON_MISSING folder.\n"
        f"See {outfile_without_json} and {outfile_with_json} for full "
         "information.\n"
    )

if __name__ == "__main__":
    main()
