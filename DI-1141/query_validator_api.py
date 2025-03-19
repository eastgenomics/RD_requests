import argparse
import os
import pandas as pd
import requests
import sys

from ratelimit import limits, sleep_and_retry
from tenacity import retry, stop_after_attempt, wait_exponential
from tqdm import tqdm


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    Returns
    ----------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(description="Query VariantValidator API")

    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        help=(
            "File with any mismatches between two VEP versions for "
            "consequence, HGVSc and HGVSp"
        ),
    )

    parser.add_argument(
        "-g",
        "--genome_build",
        type=str,
        choices=["GRCh37", "GRCh38"],
        help="Genome build the samples were run in",
    )

    parser.add_argument(
        "-a",
        "--assay",
        type=str,
        help="Assay name to be used in the output file",
    )

    parser.add_argument(
        "-c",
        "--chunks",
        type=int,
        help="Number of variants to query in each chunk",
    )

    parser.add_argument(
        "-d",
        "--output_dir",
        type=str,
        help="Directory to save output files",
    )

    args = parser.parse_args()

    return args


def read_in_mismatches(mismatches_file, sep):
    """
    Read in the mismatches file to a dataframe and add a column for the
    variant description

    Parameters
    ----------
    mismatches_file : str
        path to file with mismatches between two VEP versions
    sep : str
        separator used in the file

    Returns
    -------
    mismatches_df : pd.DataFrame
        dataframe of mismatches with additional column for variant description
    """
    mismatches_df = pd.read_csv(mismatches_file, sep=sep)

    mismatches_df["variant_description"] = (
        mismatches_df[["CHROM", "POS", "REF", "ALT"]]
        .astype(str)
        .agg("-".join, axis=1)
    )

    return mismatches_df


@sleep_and_retry
@limits(calls=1, period=1)
@retry(
    stop=stop_after_attempt(5),  # Retry up to 5 times
    wait=wait_exponential(
        multiplier=1, min=1, max=10
    ),  # Exponential backoff (1s, 2s, 4s, etc.)
)
def fetch_endpoint(server, request, content_type):
    """
    Query an endpoint and return the response

    Parameters
    ----------
    server : str
        the base URL of the server
    request : str
        the endpoint to query
    content_type : str
        the content type to be returned

    Returns
    -------
    response: dict
        the response from the server

    Raises
    ------
    Exception
        Raised when rate limit is hit
    """
    response = requests.get(server + request, headers={"Accept": content_type})

    if response.status_code == 429:  # 429 = Too Many Requests
        print("Rate limit hit. Retrying...")
        raise Exception("Rate limit hit")  # Trigger retry

    # There are a couple of variants which trigger an internal error so print
    # and skip these
    if response.status_code == 500:  # Internal Server Error
        print(f"Skipping {request} due to server error (500).")
        return None

    if not response.ok:
        response.raise_for_status()
        sys.exit()

    if content_type == "application/json":
        return response.json()
    else:
        return response.text


def get_df_unique_variants_and_their_transcripts(mismatches_df):
    """
    Group rows for each unique variant, concatenate the transcripts
    into a single string and output new dataframe

    Parameters
    ----------
    mismatches_df : pd.DataFrame
        dataframe of mismatches

    Returns
    -------
    tx_by_variant_df : pd.DataFrame
        dataframe with each unique variant and its transcripts to query
    """
    tx_by_variant_df = (
        mismatches_df.groupby(["variant_description"])
        .agg({"Feature": lambda x: "|".join(x)})
        .reset_index()
    )

    return tx_by_variant_df


def query_variant_with_lovd(
    genome_build,
    variant_description,
    transcript_model,
    transcripts,
    checkonly,
    liftover,
):
    """
    Query the LOVD endpoint of VariantValidator

    Parameters
    ----------
    genome_build : str
        genome build - can be "GRCh37", "GRCh38", "hg19" or "hg38"
    variant_description : str
        description of variant in form chr-pos-ref-alt
    transcript_model : str
        transcript model - can be "refseq", "ensembl" or "all"
    transcripts : str
        can be "None", "all", "raw", "select", "mane", "mane_select", a
        specific transcript(s), separated by '|'
    checkonly : str
        can be "True" (do not return tx and protein descriptions), "False", or
        "tx" (stop at transcript level, exclude protein)
    liftover : str
        can be "True" (liftover to all genomic loci), "primary" (lift to
        primary assembly), or "False" (do not liftover)

    Returns
    -------
    results : dict
        response from endpoint
    """
    results = fetch_endpoint(
        "https://rest.variantvalidator.org/",
        f"LOVD/lovd/{genome_build}/{variant_description}/{transcript_model}/{transcripts}/{checkonly}/{liftover}",
        "application/json",
    )

    return results


def format_results(response):
    """
    Format the results from the API into a list with a dict for each
    annotation of a variant against a transcript

    Parameters
    ----------
    response : dict
        response from endpoint

    Returns
    -------
    results : list
        list of dicts with variant, transcript, VV HGVSc and HGVSp
    """
    results = []
    for variant, variant_info in response.items():
        if variant != "metadata":
            hgvs_info = variant_info.get(variant, {}).get("hgvs_t_and_p", {})
            for transcript, hgvs in hgvs_info.items():
                results.append(
                    {
                        "variant": variant,
                        "Feature": transcript,
                        "HGVSc_validator": hgvs.get("t_hgvs"),
                        "HGVSp_validator": hgvs.get("p_hgvs_tlc"),
                    }
                )

    return results


def chunks(list_to_chunk, n):
    """
    Get n-sized chunks from a list

    Parameters
    ----------
    list_to_chunk : list
        list to be chunked
    n : int
        size of each chunk

    Returns
    -------
    tx_by_variant_df : pd.DataFrame
        dataframe with unique variants and their transcripts
    """
    for i in range(0, len(list_to_chunk), n):
        yield list_to_chunk[i : i + n]


def query_chunks_and_output(chunked_records, output_dir, assay, genome_build):
    """
    Query each chunk of variants and write to file

    Parameters
    ----------
    chunked_records : list
        list of dicts, each with a variant to query
    output_dir : str
        name of output dir to save chunked output to
    assay : str
        assay we're looking at
    genome_build : str
        the genome build
    """
    for index, variant_chunk in enumerate(chunked_records):
        filename = f"{output_dir}/validator_{assay}_chunk{index+1}.tsv"
        if not os.path.exists(filename):
            all_results = []
            for row in tqdm(variant_chunk, desc="Querying variants"):
                results = query_variant_with_lovd(
                    genome_build,
                    row["variant_description"],
                    "refseq",
                    row["Feature"],
                    False,
                    False,
                )
                if results:
                    formatted = format_results(results)
                    all_results.extend(formatted)
            results_df = pd.DataFrame(all_results)
            results_df.to_csv(
                filename,
                index=False,
                sep="\t",
            )


def gather_vv_files_and_merge_with_vep_mismatches(
    output_dir, mismatches_df, assay
):
    """
    _summary_

    Parameters
    ----------
    output_dir : str
        name of output dir to save chunked output to
    mismatches_df : pd.DataFrame
        dataframe of VEP mismatches
    assay : str
        the assay we're looking at
    """
    validator_files = [
        f
        for f in os.listdir(output_dir)
        if f.startswith("validator") and f.endswith(".csv")
    ]
    all_validator_dfs = []
    for vv_file in validator_files:
        df = pd.read_csv(os.path.join(output_dir, vv_file), sep="\t")
        all_validator_dfs.append(df)
    merged_validator_dfs = pd.concat(all_validator_dfs, ignore_index=True)

    # Merge VV results with the any_mismatches
    final_merged = pd.merge(
        mismatches_df,
        merged_validator_dfs,
        left_on=["variant_description", "Feature"],
        right_on=["variant", "Feature"],
        how="left",
    )

    final_merged.to_csv(
        f"{output_dir}/{assay}_validator_all_merged.csv",
        sep="\t",
        index=False,
    )


def main():
    args = parse_args()
    mismatches_df = read_in_mismatches(args.input_file, "\t")
    tx_by_variant = get_df_unique_variants_and_their_transcripts(mismatches_df)
    records = tx_by_variant.to_dict("records")
    chunked = list(chunks(records, args.chunks))

    query_chunks_and_output(
        chunked, args.output_dir, args.assay, args.genome_build
    )

    gather_vv_files_and_merge_with_vep_mismatches(
        args.output_dir, mismatches_df, args.assay
    )


if __name__ == "__main__":
    main()
