import pandas as pd
import argparse
import plotly.express as px

def read_file(mismatches_tsv, variant_validator_batch):
    '''
    Read mismatches tsv from compare_vcf_hgvs script and VariantValidator batch output.
    '''
    mismatches_df = pd.read_csv(
        mismatches_tsv,
        sep='\t'
    )
    
    vv_batch_df = pd.read_csv(
        variant_validator_batch,
        header=2,
        sep='\t',
        usecols=['Input', 'HGVS_transcript', 'HGVS_Predicted_Protein']
    )
    vv_batch_df.columns=["Variant", "HGVSc_validator", "HGVSp_validator"]

    return mismatches_df, vv_batch_df


def parse_validator_batch(vv_batch_df):
    '''
    Extract the relevant information from the VariantValidator batch output and reformat
    to similar format as any_mismatches.tsv
    '''
    # get only Feature column
    vv_batch_df['Feature'] = vv_batch_df['HGVSc_validator'].str.split(':').str[0]

    # split variant output back to CHROM, POS, REF, ALT
    variant_split = vv_batch_df["Variant"].str.split("-", expand=True)
    variant_split.columns=["CHROM", "POS", "REF", "ALT"]

    # concat them
    batch_reformatted = pd.concat([variant_split, vv_batch_df.drop(columns=["Variant"])], axis=1)
    batch_reformatted['POS'] = batch_reformatted['POS'].astype(str).astype(int)

    return batch_reformatted


def merge_mismatches_batch(mismatches_df, vv_batch):
    '''
    Merge any_mismatches df with the parsed and reformatted VariantValidator batch df
    '''
    validator_merged = pd.merge(mismatches_df, vv_batch, on=["CHROM", "POS", "REF", "ALT", "Feature"], how='inner')

    return validator_merged


def convention_changes(validator_merged):
    '''
    Changing the convention of the merged df because the output of VariantValidator differs to the output of VEP
    '''
    # if ends with p.? set to .
    validator_merged.loc[validator_merged["HGVSp_validator"].str.endswith("p.?"), "HGVSp_validator"] = "."

    # remove brackets from HGVSp validator
    validator_merged["HGVSp_validator"] = validator_merged["HGVSp_validator"].str.replace(r"[()]", "", regex=True)

    # covert = in HGVSp to %3D
    validator_merged["HGVSp_validator"] = validator_merged["HGVSp_validator"].str.replace("=", "%3D")

    return validator_merged


def find_vep_validator_mismatches(merged_df, old_version, new_version):
    '''
    Compare HGVSc and HGVSp for each VEP version with VariantValidator batch output
    '''
    merged_df["HGVSc_" + old_version + "_mismatch"] = merged_df["HGVSc_" + old_version] != merged_df["HGVSc_validator"]
    merged_df["HGVSc_" + new_version + "_mismatch"] = merged_df["HGVSc_" + new_version] != merged_df["HGVSc_validator"]
    merged_df["HGVSp_" + old_version + "_mismatch"] = merged_df["HGVSp_" + old_version] != merged_df["HGVSp_validator"]
    merged_df["HGVSp_" + new_version + "_mismatch"] = merged_df["HGVSp_" + new_version] != merged_df["HGVSp_validator"]
    merged_df["HGVSc_both_mismatch"] = merged_df["HGVSc_" + old_version + "_mismatch"] & merged_df["HGVSc_" + new_version + "_mismatch"]
    merged_df["HGVSp_both_mismatch"] = merged_df["HGVSp_" + old_version + "_mismatch"] & merged_df["HGVSp_" + new_version + "_mismatch"]

    mismatch_df = merged_df[[
    "HGVSc_" + old_version + "_mismatch",
    "HGVSc_" + new_version + "_mismatch",
    "HGVSc_both_mismatch",
    "HGVSp_" + old_version + "_mismatch",
    "HGVSp_" + new_version + "_mismatch",
    "HGVSp_both_mismatch"
    ]].sum().reset_index()
    mismatch_df.columns = ["Mismatch_Type", "Count"]

    mismatches_fig = px.bar(
        mismatch_df,
        x="Mismatch_Type",
        y="Count",
        title="Mismatch Counts",
    )

    mismatches_fig.write_image("mismatches_fig.png")


def find_new_changes(merged_df, old_version, new_version):
    '''
    Find newly wrong and newly corrected changes
    '''
    # Newly corrected in v113
    merged_df["HGVSc_newly_corrected"] = (merged_df["HGVSc_" + old_version + "_mismatch"] == True) & (merged_df["HGVSc_" + new_version + "_mismatch"] == False)
    merged_df["HGVSp_newly_corrected"] = (merged_df["HGVSp_" + old_version + "_mismatch"] == True) & (merged_df["HGVSp_" + new_version + "_mismatch"] == False)

    # Newly wrong in v113
    merged_df["HGVSc_newly_wrong"] = (merged_df["HGVSc_" + old_version + "_mismatch"] == False) & (merged_df["HGVSc_" + new_version + "_mismatch"] == True)
    merged_df['HGVSp_newly_wrong'] = (merged_df["HGVSp_" + old_version + "_mismatch"] == False) & (merged_df["HGVSp_" + new_version + "_mismatch"] == True)

    # Get counts
    between_version_counts = merged_df[[
        "HGVSc_newly_wrong",
        "HGVSp_newly_wrong",
        "HGVSc_newly_corrected",
        "HGVSp_newly_corrected"
    ]].sum().reset_index()
    between_version_counts.columns = ["Mismatch_Type", "Count"]

    # Write to files
    merged_df.loc[merged_df['HGVSc_newly_wrong'] == True].to_csv('hgvsc_newly_wrong.tsv', sep='\t', index=False)
    merged_df.loc[merged_df['HGVSp_newly_wrong'] == True].to_csv('hgvsp_newly_wrong.tsv', sep='\t', index=False)
    merged_df[merged_df['HGVSc_newly_corrected'] == True].to_csv('hgvsc_newly_corrected.tsv', sep='\t', index=False)
    merged_df[merged_df['HGVSp_newly_corrected'] == True].to_csv('hgvsp_newly_corrected.tsv', sep='\t', index=False)

    # Plot
    changes_fig = px.bar(
        between_version_counts,
        x="Mismatch_Type",
        y="Count",
        title="Mismatch Counts",
    )

    changes_fig.write_image("new_changes_fig.png")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    ----------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Python script comparing any mismatches in HGVSc, HGVSp, or Consequence between 2 VEP versions produced by compare_vcf_hgvs.py script"
    )

    parser.add_argument(
        '-tsv', '--mismatches_tsv', type=str, help="full path of tsv output of compare_vcf_hgvs script (any_mismatches.tsv)"
    )
    parser.add_argument(
        '-vv', '--variant_validator', type=str, help="full path of VariantValidator batch job txt file"
    )
    parser.add_argument(
        '-old', '--old_version', type=str, help="Old VEP version"
    )
    parser.add_argument(
        '-new', '--new_version', type=str, help="New VEP version"
    )
    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    mismatches_df, vv_batch_df = read_file(args.mismatches_tsv, args.variant_validator)
    parsed_vv_batch = parse_validator_batch(vv_batch_df)
    validator_merged = merge_mismatches_batch(parsed_vv_batch, mismatches_df)
    merged_df = convention_changes(validator_merged)

    merged_df.to_csv('merged_validator_mismatches.tsv', sep='\t', header=True, index=False)

    find_vep_validator_mismatches(merged_df, args.old_version, args.new_version)
    find_new_changes(merged_df, args.old_version, args.new_version)


if __name__ == "__main__":
    main()
