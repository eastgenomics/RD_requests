import argparse
import pandas as pd
import requests
import sys

from collections import defaultdict
from openpyxl.utils import get_column_letter


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Information required to query SpliceAI and GTEx"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Excel file of genes to query",
    )

    parser.add_argument(
        "-t",
        "--tissues",
        required=True,
        nargs="+",
        help=(
            "One or more GTEx tissueSiteDetailId values to get median gene"
            " expression for"
        ),
    )

    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Name of output file"
    )

    return parser.parse_args()


def query_api(server: str, endpoint: str, params: dict) -> dict:
    """
    Send a GET request to a REST API endpoint and return the JSON response.

    Parameters
    ----------
    server : str
        Base URL of the API server (e.g., "https://gtexportal.org/api/v2")
    endpoint : str
        API endpoint path to append to the server URL
    params : dict
        Dictionary of query parameters to include in the request

    Returns
    -------
    dict
        Parsed JSON response from the API
    """
    response = requests.get(server + endpoint, params)

    if not response.ok:
        response.raise_for_status()
        sys.exit()
    return response.json()


def get_gtex_data(excel_file: list, tissues: list) -> pd.DataFrame:
    """
    Retrieve median gene expression data from the GTEx API for genes
    listed in the input DataFrame and for certain tissues. This requires
    querying to get the gene ID (and version) used in the relevant GTEx dataset version.

    Parameters
    ----------
    excel_file : pd.DataFrame
        Input DataFrame containing at least a Gene column with
        gene symbols
    tissues: list[str]
        list of tissue IDs to search

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per gene symbol and median expression
        values for the requested tissues
    """
    symbols = excel_file["Gene"].unique()
    info_dict = defaultdict(dict)
    for symbol in symbols:
        symbol = symbol.strip()
        gene_info = query_api(
            server="https://gtexportal.org/api/v2",
            endpoint="/reference/geneSearch",
            params={"geneId": symbol, "gencodeVersion": "v39"},
        )
        gene_id = None
        for gene in gene_info["data"]:
            if gene["geneSymbol"] == symbol:
                gene_id = gene["gencodeId"]
        info_dict[symbol]["gene_id"] = gene_id

    gene_ids = [v["gene_id"] for k, v in info_dict.items()]

    gtex_response = query_api(
        server="https://gtexportal.org/api/v2/",
        endpoint="expression/medianGeneExpression",
        params={
            "datasetId": "gtex_v10",
            "gencodeId": gene_ids,
            "tissueSiteDetailId": tissues,
        },
    )
    df = pd.DataFrame(gtex_response["data"])
    gtex_data = df.pivot(
        index="geneSymbol", columns="tissueSiteDetailId", values="median"
    ).reset_index()
    return gtex_data


def get_spliceai_scores(excel_file):
    """
    Get SpliceAI scores (GRCh38) for variants, using the transcript provided
    in the input file, or the MS transcript as fallback if that does not
    exist in SpliceAI.

    Parameters
    ----------
    excel_file : pd.DataFrame
        Input DataFrame containing at least the following columns:
        - "Location g. (GRCh38)"
        - "Variant c."

    Returns
    -------
    pd.DataFrame
            DataFrame containing SpliceAI scores per variant, including:
            - Variant
            - Input_HGVSc
            - SpliceAI_transcript_set
            - SpliceAI_DS_AL, SpliceAI_DP_AL
            - SpliceAI_DS_DL, SpliceAI_DP_DL
            - SpliceAI_DS_AG, SpliceAI_DP_AG
            - SpliceAI_DS_DG, SpliceAI_DP_DG
    """
    results = []

    for idx, row in excel_file.iterrows():
        variant = row["Location g. (GRCh38)"].replace(":", "-")
        refseq_target_full = row["Variant c."]
        refseq_target_base = refseq_target_full.split(".")[0]

        # Query SpliceAI
        response = query_api(
            server="https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/",
            endpoint="",
            params={"hg": "38", "variant": variant},
        )

        scores_list = response.get("scores", [])
        matched = False

        for score in scores_list:
            refseq_ids = score.get("t_refseq_ids")
            if not refseq_ids:
                continue
            refseq_ids_base = [r.split(".")[0] for r in refseq_ids]

            # Only keep score if RefSeq matches (ignoring version)
            if refseq_target_base in refseq_ids_base:
                results.append(
                    {
                        "Variant": variant.replace("-", ":"),
                        "Input_HGVSc": refseq_target_full,
                        "SpliceAI_transcript_set": ",".join(refseq_ids_base),
                        "SpliceAI_DS_AL": score.get("DS_AL"),
                        "SpliceAI_DP_AL": score.get("DP_AL"),
                        "SpliceAI_DS_DL": score.get("DS_DL"),
                        "SpliceAI_DP_DL": score.get("DP_DL"),
                        "SpliceAI_DS_AG": score.get("DS_AG"),
                        "SpliceAI_DP_AG": score.get("DP_AG"),
                        "SpliceAI_DS_DG": score.get("DS_DG"),
                        "SpliceAI_DP_DG": score.get("DP_DG"),
                    }
                )
                matched = True
                break

        if not matched:
            ms_scores = [s for s in scores_list if s.get("t_priority") == "MS"]
            if ms_scores:
                score = ms_scores[0]
                ms_refseqs = score.get("t_refseq_ids") or ["NA"]
                cleaned = [tx.split(".")[0] for tx in ms_refseqs]
                results.append(
                    {
                        "Variant": variant.replace("-", ":"),
                        "Input_HGVSc": refseq_target_full,
                        "SpliceAI_transcript_set": ",".join(cleaned),
                        "SpliceAI_DS_AL": score.get("DS_AL"),
                        "SpliceAI_DP_AL": score.get("DP_AL"),
                        "SpliceAI_DS_DL": score.get("DS_DL"),
                        "SpliceAI_DP_DL": score.get("DP_DL"),
                        "SpliceAI_DS_AG": score.get("DS_AG"),
                        "SpliceAI_DP_AG": score.get("DP_AG"),
                        "SpliceAI_DS_DG": score.get("DS_DG"),
                        "SpliceAI_DP_DG": score.get("DP_DG"),
                    }
                )
                print(
                    f"Variant {variant} had no match; using MS"
                    f" transcript(s) {ms_refseqs} instead."
                )

    spliceai_df = pd.DataFrame(results)
    spliceai_df = spliceai_df.drop_duplicates(keep="first")

    return spliceai_df


def format_final_file(merged_df: pd.DataFrame, tissues: list) -> pd.DataFrame:
    """
    Format the final Excel file by dropping unnecessary columns and
    re-ordering.

    Parameters
    ----------
    merged_df : pd.DataFrame
        dataframe with both SpliceAI and GTEx data
    tissues: list[str]
        list of tissue IDs which were searched

    Returns
    -------
    pd.DataFrame
        dataframe with final columns in required order
    """
    merged_splice_ai = merged_df.drop(
        columns=[
            "Variant",
            "Input_HGVSc",
            "Δ  score ",
            "position",
            "Δ  score .1",
            "position.1",
            "Δ  score .2",
            "position.2",
            "Δ  score .3",
            "position.3",
            "geneSymbol",
            "Cultured fibroblasts",
            "EBV-transformed lymphocytes",
            "Whole blood",
        ]
    )
    # Add GTEx to beginning of those columns
    gtex_cols = tissues
    gtex_rename = {
        col: f"GTEx_median_gene_expression_{col}" for col in gtex_cols
    }
    merged_splice_ai = merged_splice_ai.rename(columns=gtex_rename)

    # Reorder: original excel cols, then SpliceAI, then GTEx
    original_cols = [
        c
        for c in merged_splice_ai.columns
        if not c.startswith("SpliceAI_") and not c.startswith("GTEx_")
    ]
    spliceai_final_cols = [
        c for c in merged_splice_ai.columns if c.startswith("SpliceAI_")
    ]
    gtex_final_cols = [
        c for c in merged_splice_ai.columns if c.startswith("GTEx_")
    ]
    merged_splice_ai = merged_splice_ai[
        original_cols + spliceai_final_cols + gtex_final_cols
    ]

    return merged_splice_ai


def main():
    args = parse_args()
    all_sheets = pd.read_excel(args.input, sheet_name=None, header=2)
    sheet_names = list(all_sheets.keys())
    first_sheet_name = sheet_names[0]

    excel_file = all_sheets[first_sheet_name].copy()
    excel_file["Variant c."] = excel_file["Variant c."].str.strip()

    gtex_data = get_gtex_data(excel_file, args.tissues)

    merged = pd.merge(
        excel_file,
        gtex_data,
        left_on="Gene",
        right_on="geneSymbol",
        how="left",
    )

    spliceai_df = get_spliceai_scores(excel_file)

    merged_splice_ai = pd.merge(
        merged,
        spliceai_df,
        how="left",
        left_on=["Variant c.", "Location g. (GRCh38)"],
        right_on=["Input_HGVSc", "Variant"],
    )
    merged_splice_ai = format_final_file(merged_splice_ai, args.tissues)
    all_sheets[first_sheet_name] = merged_splice_ai

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        for name, df in all_sheets.items():
            df.to_excel(writer, sheet_name=name, index=False)

            worksheet = writer.sheets[name]

            # Set all column widths to ~5 cm
            for i in range(1, len(df.columns) + 1):
                col_letter = get_column_letter(i)
                worksheet.column_dimensions[col_letter].width = 16.0


if __name__ == "__main__":
    main()
