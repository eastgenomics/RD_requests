import argparse
import pandas as pd
import uuid

def parse_args():
    """
    Parse command line arguments

    Returns:
        args (Namespace): Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--path",
        type=str,
        help=("Path to workbook"),
    )

    return parser.parse_args()

def make_variants_df_from_wgs_excel(path):
    """
    Args:
        path (str): path to workbook
    """
    df = pd.read_excel(path, usecols=[
        "LastModifiedDate",
        "Clinical_indication",
        "Build",
        "Chromosome",
        "Start",
        "Reference",
        "Alternate",
        "Classification",
        "Comment"
        ])
    df = df.rename(columns={
        "Chromosome": "chromosome",
        "Start": "start",
        "Reference": "reference_allele",
        "Alternate": "alternate_allele",
        "Classification": "germline_classification",
        "LastModifiedDate": "date_last_evaluated",
        "Clinical_indication": "preferred_condition_name",
    })
    df = convert_classification_to_valid_classification_for_clinvar(df)
    print(df["germline_classification"])
    df = add_consistent_columns(df)
    print(df["collection_method"])
    df = add_uuid(df)
    print(df['local_id'])
    print(df['linking_id'])

def add_consistent_columns(df):
    """
    Add columns to df, that have the same value for all variants

    Args:
        df (pd.DataFrame): variant df
    Returns:
        df (pd.DataFrame): variant df, with added columns.
    """
    df["collection_method"] = 'clinical testing'
    df["allele_origin"] = 'germline'
    df["affected_status"] = "yes"
    df["organisation_id"] = "288359"
    return df


def convert_classification_to_valid_classification_for_clinvar(df):
    """
    Replace values in the classification column of the df with valid
    classifications for ClinVar submission
    Args:
        df (pd.DataFrame): variant df
    Returns:
        df (pd.DataFrame): variant df, with updated Classification column.
    """
    valid = {
        "pathogenic_variant": "Pathogenic",
        "likely_pathogenic_variant": "Likely pathogenic",
        "variant_of_unknown_clinical_significance": "Uncertain significance",
        "likely_benign_variant": "Likely benign",
        "benign_variant": "Benign"
    }
    df = df.replace({"germline_classification": valid})
    return df


def add_uuid(df):
    df['local_id'] = [f"uid_{uuid.uuid1().time}" for _ in range(len(df.index))]
    df['linking_id'] = df.loc[:, 'local_id']
    return df

def main():
    args = parse_args()
    make_variants_df_from_wgs_excel(args.path)



if __name__ == "__main__":
    main()
