import pandas as pd
import argparse
import sys
import os
import re
import numpy as np

def get_sample_pairs(directory):
    '''
    Read sample pairs in directory.

    Parameters
    ----------
    directory : str
        Path to directory with VCFs of both VEP versions

    Returns
    -------
    sample_pairs : list
        List of filenames of sample pairs (new and old VEP)
    '''
    files = sorted(f for f in os.listdir(directory) if f.endswith("_filtered.vcf"))
    sample_ids = {}

    # find files with regex pattern and group
    for file in files:
        match = re.match(r"^(\d{9}-\d{5}[A-Z]\d{4})-.*_vep(\d+)\_filtered.vcf$", file) 
        if match:
            sample_name = match.group(1)
            sample_ids.setdefault(sample_name, []).append(file)

    # make sure 2 versions (old and new vep) per sample
    sample_pairs = [(files[0], files[1]) for sample, files in sample_ids.items() if len(files) == 2]
    
    return sample_pairs


def read_file(old_vcf, new_vcf, old_version, new_version, directory):
    '''
    Read old and new VEP version VCFs into pandas dataframes.

    Parameters
    ----------
    old_vcf : str
        Filename of the sample annotated with old VEP.
    new_vcf : str
        Filename of the sample annotated with new VEP.
    old_version : str
        Old VEP version.
    new_version : str
        New VEP version.
    directory : str
        Path to directory with VCFs of both VEP versions

    Returns
    -------
    annotations_old : pd.DataFrame
        Dataframe of the sample annotated with old VEP.
    annotations_new : pd.DataFrame
        Dataframe of the sample annotated with new VEP.
    '''
    old_vep_file = os.path.join(directory, old_vcf)
    new_vep_file = os.path.join(directory, new_vcf)

    annotations_old = pd.read_csv(
        old_vep_file,
        sep='\t',
        names=["CHROM", "POS", "REF", "ALT", f"Consequence_{old_version}", "Feature", f"HGVSc_{old_version}", f"HGVSp_{old_version}"],
        dtype={'CHROM': 'str'}
    )
    
    annotations_new = pd.read_csv(
        new_vep_file,
        sep='\t',
        names=["CHROM", "POS", "REF", "ALT", f"Consequence_{new_version}", "Feature", f"HGVSc_{new_version}", f"HGVSp_{new_version}"],
        dtype={'CHROM': 'str'}
    )

    return annotations_old, annotations_new


def find_mismatches(annotations_old, annotations_new, old_version, new_version):
    '''
    Merge old and new variant lines on variant and feature.
    Find lines where Consequence, HGVSc, or HGVSp don't match between old and new.
    Remove NR transcripts.

    Parameters
    ----------
    annotations_old : pd.DataFrame
        Dataframe of the sample annotated with old VEP.
    annotations_new : pd.DataFrame
        Dataframe of the sample annotated with new VEP.
    old_version : str
        Old VEP version.
    new_version : str
        New VEP version.


    Returns
    -------
    mismatches_not_nr : pd.DataFrame
        Dataframe of mismatches between old and new VEP annotations. 
    '''
    # merge based on certain headers
    merged = pd.merge(annotations_old, annotations_new, on=["CHROM", "POS", "REF", "ALT", "Feature"], how='inner')
    
    # get any mismatches between old and new Consequence, HGVSc, and HGVSp
    any_mismatches = merged[(merged[f"Consequence_{old_version}"] != merged[f"Consequence_{new_version}"]) 
                            | (merged[f"HGVSc_{old_version}"] != merged[f"HGVSc_{new_version}"]) 
                            | (merged[f"HGVSp_{old_version}"] != merged[f"HGVSp_{new_version}"])]
    
    # drop any NR transcripts
    mismatches_not_nr = any_mismatches[~any_mismatches['Feature'].str.startswith('NR_')]

    return mismatches_not_nr


def to_csv_chunks(df, output_file, sep=",", chunksize=25000):
    '''
    Save pandas dataframe to csv with max 25000 rows. If there are more than 25000 rows,
    create files with similar number of rows.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to write out.
    output_file : str
        Filename to output.
    sep : str
        Separator between each column per row in the dataframe. 

    Returns
    -------
    uniq_vars : CSV file
        CSV file with the unique variants from the mismatches for all samples. 
    uniq_transcripts : CSV file
        CSV file with the unique transcripts from the mismatches for all samples. 
    '''
    total_rows = len(df)
    num_csv = (total_rows // chunksize) + (1 if total_rows % chunksize else 0)

    for i, chunk in enumerate(np.array_split(df, num_csv)):
        chunk.to_csv(f"{output_file}_{i}.txt", sep=sep, index=False, header=False)


def create_variantvalidator_inputs(all_dfs):
    '''
    Get unique variants and transcripts for VariantValidator input. Split to 25000 variant per file.

    Parameters
    ----------
    all_dfs : pd.DataFrame
        Concatenated dataframe of mismatches between new and old VEP for all samples.

    Returns
    -------
    uniq_vars : CSV file
        CSV file with the unique variants from the mismatches for all samples. 
    uniq_transcripts : CSV file
        CSV file with the unique transcripts from the mismatches for all samples. 
    '''
    # get unique variants
    chr_pos_ref_alt = all_dfs[['CHROM','POS', 'REF', 'ALT']]
    uniq_vars = chr_pos_ref_alt.drop_duplicates()
    to_csv_chunks(uniq_vars, "unique_variants", "-")

    # get unique transcripts
    transcript_only = all_dfs[['Feature']]
    uniq_transcripts = transcript_only.drop_duplicates()
    to_csv_chunks(uniq_transcripts, "unique_transcripts")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    ----------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Python script comparing VCFs with different vep versions"
    )

    parser.add_argument(
        '-dir', '--directory', type=str, help="Directory that contains all VCFs"
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

    sample_pairs = get_sample_pairs(args.directory)
    all_dfs=[]

    for old_vcf,new_vcf in sample_pairs:
        annotations_old, annotations_new = read_file(old_vcf, new_vcf, args.old_version, args.new_version, args.directory)
        all_dfs.append(find_mismatches(annotations_old, annotations_new, args.old_version, args.new_version))
    
    all_dfs = pd.concat(all_dfs)
    all_dfs = all_dfs.drop_duplicates()
    all_dfs.to_csv('any_mismatches.tsv', sep='\t', header=True, index=False)
    
    create_variantvalidator_inputs(all_dfs)


if __name__ == "__main__":
    main()
