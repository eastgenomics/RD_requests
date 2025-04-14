import argparse
import concurrent.futures
import dxpy as dx
import numpy as np
import pandas as pd
import subprocess

from collections import defaultdict
from mergedeep import merge
from pathlib import Path


def parse_args():
    """
    Parse the command line arguments inputs given

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Information necessary to compare filtering"
    )

    parser.add_argument(
        "--dx_old",
        type=str,
        required=True,
        help=(
            "Path to folder with first set of filtered VCFs in DNAnexus. "
            "Format: project-XXX:/folder"
        ),
    )

    parser.add_argument(
        "--dx_new",
        type=str,
        required=True,
        help=(
            "Path to folder with second set of filtered VCFs in DNAnexus. "
            "Format: project-XXX:/folder"
        ),
    )

    parser.add_argument(
        "--old_path",
        type=str,
        required=True,
        help="Full path where the old filtering VCFs will be downloaded to",
    )

    parser.add_argument(
        "--new_path",
        type=str,
        required=True,
        help="Full path where the old filtering VCFs will be downloaded to",
    )

    parser.add_argument(
        "--old_name",
        type=str,
        required=True,
        help="Name of the column for the old filtering",
    )

    parser.add_argument(
        "--new_name",
        type=str,
        required=True,
        help="Name of the column for the new filtering",
    )

    parser.add_argument(
        "--fields",
        type=str,
        required=True,
        help=(
            "Comma-separated string of fields to get from the PASS variants,"
            " e.g. 'CSQ_HGVSc,MOI,CSQ_gnomADe_AF'"
        ),
    )

    parser.add_argument(
        "--outfile",
        type=str,
        required=True,
        help="Name of output file with filtering comparison per sample",
    )

    args = parser.parse_args()

    return args


def find_files_in_project(search_folder):
    """
    Find all of the soft-filtered VCF files in the relevant DNAnexus folder

    Parameters
    ----------
    search_folder : str
        the folder to search for VCFs in, format project-XX:/folder
    Returns
    -------
    files : list
        list of dicts, each dict containing info about one file
    """
    project_id, folder = search_folder.split(":")
    files = list(
        dx.find_data_objects(
            project=project_id,
            folder=folder,
            name="*optimised_filtered.vcf.gz",
            name_mode="glob",
            classname="file",
            describe={
                "fields": {
                    "name": True,
                    "createdBy": True,
                }
            },
        )
    )

    print(f"Found {len(files)} files in folder {search_folder}")

    return files


def make_vcf_dict(filtered_vcfs, vcf_filter_type):
    """
    Make a dictionary of VCFs with sample as key and files for that sample
    for each clinical indication as values

    Parameters
    ----------
    filtered_vcfs : list
        list of dicts, each dict containing info about one VCF file
    vcf_filter_type : str
        'old' or 'new' depending on whether the VCFs have old or new filtering
        applied

    Returns
    -------
    vcf_file_dict : dict
        dict with sample as key and file info as nested dict
    Example (fake sample IDs):
    {
        '25010R0123': {
            'R58.4_Adult onset neurodegenerative disorder_P': {
                'old': [
                    {
                        'file_id': 'file-XYZ',
                        'name': '123456789-25010R0123-25NGWES5-9526-F-103698_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                        'name_with_panel': '123456789-25010R0123-25NGWES5-9526-F-103698_markdup_recalibrated_Haplotyper_annotated_R158.1.optimised_filtered.vcf.gz'
                    }
                ]
            }
        },
        '25020R0123'...
    }
    """
    vcf_file_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    # For each file, get panel it was tested for and add file
    # info to dict under the relevant sample and panel
    for vcf_file in filtered_vcfs:
        filename = vcf_file["describe"]["name"]
        specimen_id = filename.split("-")[1]

        # Get the panel the sample was tested for via the job input
        job_id = vcf_file["describe"]["createdBy"]["job"]
        panel = dx.bindings.dxjob.DXJob(dxid=job_id).describe()["runInput"][
            "panel_string"
        ]

        # Get R code of panel to update filename for the VCF
        r_code = panel.split("_")[0]
        file_base = filename.removesuffix(".optimised_filtered.vcf.gz")
        updated_filename = f"{file_base}_{r_code}.optimised_filtered.vcf.gz"

        vcf_file_dict[specimen_id][panel][vcf_filter_type].append(
            {
                "file_id": vcf_file["id"],
                "name": filename,
                "name_with_panel": updated_filename,
            }
        )

    print(
        f"Found {len(vcf_file_dict.keys())} unique samples for "
        f"{vcf_filter_type}"
    )

    return vcf_file_dict


def create_sample_file_dict(old_filter_path, new_filter_path):
    """
    Create a single dict containing info for each sample, listing VCF
    files for the sample and panel for both old and new filtering
    Parameters
    ----------
    old_filter_path : str
        name of folder in DNAnexus for the old filter VCFs
    new_filter_path : str
        name of folder in DNAnexus for the new filter VCFs

    Returns
    -------
    old_vcf_dict : dict
        dict with sample as key and file info for old + new filtering
        as nested dict
    Example:
    {
        '25010R0123': {
            'R58.4_Adult onset neurodegenerative disorder_P': {
                'old': [
                    {
                        'file_id': 'file-XYZ',
                        'name': '123456789-25010R0123-25NGWES5-9526-F-103698_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                        'name_with_panel': '123456789-25010R0123-25NGWES5-9526-F-103698_markdup_recalibrated_Haplotyper_annotated_R158.1.optimised_filtered.vcf.gz'
                    }
                ],
                'new': [
                    {
                        'file_id': 'file-ABC',
                        'name': '123456789-25010R0123-25NGWES5-9526-F-103698_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                        'name_with_panel': '123456789-25010R0123-25NGWES5-9526-F-103698_markdup_recalibrated_Haplotyper_annotated_R158.1.optimised_filtered.vcf.gz'
                    }
                ]
            }
        },
        '25020R0123'...
    }
    """
    # Find all the files from the folder containing the VCFs with old
    # filtering and then new filtering
    old_vcfs = find_files_in_project(old_filter_path)
    new_vcfs = find_files_in_project(new_filter_path)

    # Make a dict from each and then merge the dicts
    old_vcf_dict = make_vcf_dict(old_vcfs, "old")
    new_vcf_dict = make_vcf_dict(new_vcfs, "new")
    merge(old_vcf_dict, new_vcf_dict)

    return old_vcf_dict


def check_file_duplicates(vcf_dict):
    """
    Check for samples with duplicates files for the same clinical
    indication, only keep one set of old and new files per clinical indication
    per sample

    Parameters
    ----------
    vcf_dict : dict
        dict with sample as key, with nested dicts for each clinical indication
        per sample containing lists of files with old and new filtering

    Returns
    -------
    updated_dict : dict
        dict with sample as key and nested dicts for each clinical indication
        for that sample (a single dict per old and new filtering)
    """
    updated_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for sample, sample_info in vcf_dict.items():
        for clin_ind in sample_info:
            # Check how many 'old' filtering files there are for that clin ind
            no_of_old_files = len(sample_info[clin_ind]["old"])
            if no_of_old_files > 1:
                print(
                    f"Warning - {sample} has {no_of_old_files} VCFs for old"
                    f" filtering found for clinical indication {clin_ind}."
                    " Keeping first instance"
                    f" {sample_info[clin_ind]['old'][0]['file_id']}"
                )

            updated_dict[sample][clin_ind]["old"] = sample_info[clin_ind][
                "old"
            ][0]
            # Check how many 'new' filtering files there are for that clin ind
            no_of_new_files = len(sample_info[clin_ind]["new"])
            if no_of_new_files > 1:
                print(
                    f"Warning - {sample} has {no_of_new_files} VCFs for new "
                    f"filtering found for clinical_indication {clin_ind}."
                    " Keeping first instance"
                    f" {sample_info[clin_ind]['new'][0]['file_id']}"
                )

            updated_dict[sample][clin_ind]["new"] = sample_info[clin_ind][
                "new"
            ][0]

    return updated_dict


def download_files(old_path, new_path, sample_info):
    """
    If each file isn't already downloaded in the specified local directory,
    download there

    Parameters
    ----------
    old_path : str
        full path to where the VCFs with old filtering will be downloaded to
    new_path : str
        full path to where the VCFs with new filtering will be downloaded to
    sample_info : dict
        the dictionary holding file info for a single sample
    """
    Path(old_path).mkdir(parents=True, exist_ok=True)
    Path(new_path).mkdir(parents=True, exist_ok=True)

    for clin_ind in sample_info:
        old_file_id = sample_info[clin_ind]["old"]["file_id"]
        old_file_name = sample_info[clin_ind]["old"]["name_with_panel"]
        file_path = Path(f"{old_path}/{old_file_name}")

        # Rename any files to include panel when compared to DNAnexus
        if not file_path.exists():
            dx.download_dxfile(old_file_id, f"{old_path}/{old_file_name}")

        # Do same for the new filtering files in the separate folder
        new_file_id = sample_info[clin_ind]["new"]["file_id"]
        new_file_name = sample_info[clin_ind]["new"]["name_with_panel"]
        file_path = Path(f"{new_path}/{new_file_name}")

        if not file_path.exists():
            dx.download_dxfile(new_file_id, f"{new_path}/{new_file_name}")


def concurrent_download(vcf_dict, workers, old_path, new_path):
    """
    Concurrently download VCF files to respective folders locally

    Parameters
    ----------
    vcf_dict : dict
        final dict with GM number as key and dict with single VCF and sample
        info as value
    workers : int
        number of workers
    old_path : str
        full path to where the VCFs with old filtering will be downloaded to
    new_path : str
        full path to where the VCFs with new filtering will be downloaded to
    """
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=workers
    ) as executor:
        concurrent_jobs = {
            executor.submit(
                download_files, old_path, new_path, sample_info
            ): sample_info
            for sample_info in vcf_dict.values()
        }
        for future in concurrent.futures.as_completed(concurrent_jobs):
            try:
                data = future.result()
            except Exception as exc:
                print(
                    "Error downloading data for "
                    f"{concurrent_jobs[future]}: {exc}"
                )


def get_PASS_variants(folder, name, fields):
    """
    Get INFO fields for the PASS variants so we can understand why they
    are filtered in

    Parameters
    ----------
    folder : str
        name of the folder the file is in
    name : str
        name of the file
    fields: str
        comma-separated string of fields to get from the PASS variants

    Returns
    -------
    variants : str
        a string containing all of the PASS variants, each separated by a
        newline character so they display nicely later in an Excel spreadsheet

    Raises
    ------
    RuntimeError
        If there is an error with the bcftools query command
    """
    filepath = f"{folder}/{name}"
    fields_list = ["%" + item for item in fields.split(",")]
    formatted_fields = "\t".join(fields_list)

    variant_output = subprocess.run(
        "bcftools query -i'FILTER=\"PASS\"' -f"
        f"'%CHROM\t%POS\t%REF\t%ALT\t{formatted_fields}\n' {filepath}",
        shell=True,
        capture_output=True,
    )

    if variant_output.returncode != 0:
        raise RuntimeError(
            f"Bcftools query error for {filepath}. Error:"
            f" {variant_output.stderr.decode()}"
        )

    variants = variant_output.stdout.decode()

    return variants


def get_total_variants(folder, vcf_name):
    """
    Get total number of variants in a VCF file

    Parameters
    ----------
    folder : str
        The full path of the folder containing the VCFs (old or new filtering)
    vcf_name : str
        The name of the VCF file

    Returns
    -------
    variants : int
        The number of variants in the VCF

    Raises
    ------
    RuntimeError
        If there is an error with the grep command
    """
    filepath = f"{folder}/{vcf_name}"

    variant_output = subprocess.run(
        f"zgrep -v ^# {filepath} | wc -l", shell=True, capture_output=True
    )

    if variant_output.returncode != 0:
        raise RuntimeError(
            f"Grep error for {filepath}. Error:"
            f" {variant_output.stderr.decode()}"
        )

    variants = variant_output.stdout.decode()

    return variants


def get_variant_info(vcf_dict, old_path, new_path, fields):
    """
    Adds info about total number of variants in the panel and the
    variants and counts filtered in with old and new filtering

    Parameters
    ----------
    vcf_dict : dict
        dict with sample as key and nested dicts for each clinical indication
        for that sample with a single dict per old and new filtering
    old_path : str
        path to folder holding all the VCFs with old filtering
    new_path : str
        path to folder holding all the VCFs with new filtering
    fields: str
        comma-separated string of fields to get from the PASS variants

    Returns
    -------
    vcf_dict : dict
        dict with sample as key and nested dicts for each clinical indication
        for that sample with added variants and counts
    """
    for sample, sample_info in vcf_dict.items():
        for clin_ind in sample_info:
            # Get PASS variants with old filtering and add to dict
            old_filename = sample_info[clin_ind]["old"]["name_with_panel"]
            old_variants = get_PASS_variants(old_path, old_filename, fields)
            sample_info[clin_ind]["old"]["variants"] = old_variants

            # Add count from this to dict
            old_variant_list = old_variants.split("\n")
            old_var_list = list(filter(None, old_variant_list))
            sample_info[clin_ind]["old"]["variant_count"] = len(old_var_list)
            total_old = get_total_variants(old_path, old_filename)
            sample_info[clin_ind]["old"]["total_variants"] = int(
                total_old.rstrip("\n")
            )

            # Get PASS variants with new filtering and add to dict
            new_filename = sample_info[clin_ind]["new"]["name_with_panel"]
            new_variants = get_PASS_variants(new_path, new_filename, fields)
            sample_info[clin_ind]["new"]["variants"] = new_variants
            # Add count from this to dict
            new_variant_list = new_variants.split("\n")
            new_var_list = list(filter(None, new_variant_list))
            sample_info[clin_ind]["new"]["variant_count"] = len(new_var_list)

            total_new = get_total_variants(new_path, new_filename)
            sample_info[clin_ind]["new"]["total_variants"] = int(
                total_new.rstrip("\n")
            )

    return vcf_dict


def create_df_one_row_per_sample(vcf_dict, old_name, new_name):
    """
    Create a df of each sample, panel and variants

    Parameters
    ----------
    vcf_dict : dict
        dict of each sample and the routine and optimised variants and counts
    old_name : str
        name of the column for the old filtering
    new_name : str
        name of the column for the new filtering

    Returns
    -------
    variant_df : pd.DataFrame
        dataframe from the above dict
    """
    df_rows = []
    for sample, file_info in vcf_dict.items():
        for clin_ind in file_info:
            old_included = file_info[clin_ind]["old"]["variants"].rstrip("\n")
            old_count = file_info[clin_ind]["old"]["variant_count"]
            old_total = file_info[clin_ind]["old"]["total_variants"]

            new_included = file_info[clin_ind]["new"]["variants"].rstrip("\n")
            new_count = file_info[clin_ind]["new"]["variant_count"]
            new_total = file_info[clin_ind]["new"]["total_variants"]

            df_rows.append(
                {
                    "Sample": sample,
                    "CI": clin_ind,
                    f"{old_name}_panel_total": old_total,
                    f"{old_name}_included_vars": old_included,
                    f"{old_name}_included_count": old_count,
                    f"{new_name}_panel_total": new_total,
                    f"{new_name}_included_vars": new_included,
                    f"{new_name}_included_count": new_count,
                    "Difference": new_count - old_count,
                }
            )

    variant_df = pd.DataFrame(df_rows)

    return variant_df


def create_tuple_for_multiindex(col_name, old_name, new_name):
    """
    Create tuples of index and column names for multiindex

    Parameters
    ----------
    col_name : str
        name of a column in the dataframe
    old_name : str
        name for old filtering, e.g. 'GRCh37'
    new_name : str
        name for new filtering, e.g. 'GRCh38'

    Returns
    -------
    tuple
        tuple of two strings, the first is the name of the filtering
        (old or new) and the second is the name of the column
    """
    if col_name.startswith(f"{old_name}_"):
        return (old_name, col_name.removeprefix(f"{old_name}_"))
    elif col_name.startswith(f"{new_name}_"):
        return (new_name, col_name.removeprefix(f"{new_name}_"))


def create_df_multiple_rows_per_sample(variant_df, old_name, new_name, fields):
    """
    Create a dataframe with one row per sample, with counts and all
    variants filtered in

    Parameters
    ----------
    variant_df : pd.DataFrame
        dataframe of each sample, panel and variant info
    old_name : str
        name for old filtering, e.g. 'GRCh37'
    new_name : str
        name for new filtering, e.g. 'GRCh38'
    fields: str
        comma-separated string of fields to get from the PASS variants

    Returns
    -------
    df_merged : pd.DataFrame
        dataframe where the sample is the index and the same variants per
        sample are merged by HGVSc notation
    """
    old_column_name = f"{old_name}_included_vars"
    new_column_name = f"{new_name}_included_vars"

    # Get dfs of old variants and new variants per sample separately
    old_subset = variant_df[["Sample", old_column_name]]
    new_subset = variant_df[["Sample", new_column_name]]

    # Explode out the variants so we have one row per variant per sample
    old_explode = old_subset.assign(
        **{
            old_column_name: (
                old_subset[f"{old_name}_included_vars"].str.split("\n")
            )
        }
    ).explode(old_column_name)
    new_explode = new_subset.assign(
        **{
            new_column_name: (
                new_subset[f"{new_name}_included_vars"].str.split("\n")
            )
        }
    ).explode(new_column_name)

    fields_list = [field.replace("CSQ_", "") for field in fields.split(",")]
    fields_list = ["CHROM", "POS", "REF", "ALT"] + fields_list
    old_fields = [f"{old_name}_{field}" for field in fields_list]
    new_fields = [f"{new_name}_{field}" for field in fields_list]
    old_explode[old_fields] = old_explode[old_column_name].str.split(
        "\t", expand=True
    )
    new_explode[new_fields] = new_explode[new_column_name].str.split(
        "\t", expand=True
    )
    old_explode.replace(["", None], np.nan, inplace=True)
    new_explode.replace(["", None], np.nan, inplace=True)

    sample_merged_variants = pd.merge(
        old_explode,
        new_explode,
        left_on=["Sample", f"{old_name}_HGVSc"],
        right_on=["Sample", f"{new_name}_HGVSc"],
        how="outer",
    )
    sample_merged_variants.sort_values(
        by="Sample", inplace=True, ignore_index=True
    )

    # Remove cases where an empty row is duplicated due to a mismatch
    non_sample_cols = [
        col for col in sample_merged_variants.columns if col != "Sample"
    ]

    # Check cases where all rows empty
    sample_merged_variants["is_empty"] = (
        sample_merged_variants[non_sample_cols].isna().all(axis=1)
    )

    # Remove rows per sample where there are duplicated rows per sample
    # and they are flagged as empty
    df_merged = sample_merged_variants.loc[
        ~(
            sample_merged_variants.duplicated(subset=["Sample"], keep=False)
            & sample_merged_variants["is_empty"]
        )
    ]
    # Remove
    df_merged.drop(
        columns=[
            "is_empty",
            f"{old_name}_included_vars",
            f"{new_name}_included_vars",
        ],
        inplace=True,
    )

    # Index on sample so we can see variants grouped by sample
    df_merged["idx"] = df_merged.groupby("Sample").cumcount()
    df_merged.set_index(["Sample", "idx"], inplace=True)

    # Create multiindex using tuples and set that as the columns
    # this means we end up with a top level index of genome build and
    # the second level is the field name (e.g. GRCh37 top level and HGVSc
    # as lower level)
    new_columns = pd.MultiIndex.from_tuples(
        [
            create_tuple_for_multiindex(col, old_name, new_name)
            for col in df_merged.columns
        ]
    )
    df_merged.columns = new_columns
    # Create new empty Comment column to fill in manually
    df_merged["Comment"] = ""

    return df_merged


def write_to_excel_workbook(df1, df2, sheetname1, sheetname2, outfile_name):
    """
    Write out each dataframe to a separate sheet in an Excel workbook

    Parameters
    ----------
    df1 : pd.DataFrame
        dataframe to write to first sheet
    df2 : pd.DataFrame
        dataframe to write to second sheet
    sheetname1 : str
        name of first sheet
    sheetname2 : str
        name of second sheet
    outfile_name : str
        name of the output Excel file
    """
    with pd.ExcelWriter(outfile_name, mode="w", engine="openpyxl") as writer:
        df1.to_excel(writer, sheet_name=sheetname1, index=False)
        df2.to_excel(writer, sheet_name=sheetname2)


def main():
    args = parse_args()
    vcf_dict = create_sample_file_dict(args.dx_old, args.dx_new)
    vcf_dict_no_dups = check_file_duplicates(vcf_dict)
    concurrent_download(vcf_dict_no_dups, 8, args.old_path, args.new_path)
    vcf_dict_with_variants = get_variant_info(
        vcf_dict_no_dups, args.old_path, args.new_path, args.fields
    )
    one_row_per_sample_df = create_df_one_row_per_sample(
        vcf_dict_with_variants, args.old_name, args.new_name
    )
    multiple_rows_per_sample_df = create_df_multiple_rows_per_sample(
        one_row_per_sample_df, args.old_name, args.new_name, args.fields
    )
    write_to_excel_workbook(
        one_row_per_sample_df,
        multiple_rows_per_sample_df,
        "one_row_per_sample",
        "per_sample_variants_grouped",
        args.outfile,
    )


if __name__ == "__main__":
    main()
