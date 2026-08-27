import argparse
import pandas as pd
import dxpy
import re
from pathlib import Path
import pysam
import sys
import subprocess


def parse_args():
    """
    Parse command line arguments.

    Returns
    -------
    argparse.Namespace
        parsed arguments: csv, gene, position and output
    """
    parser = argparse.ArgumentParser(
        description=(
            "Search DNAnexus filtered VCFs for a given gene / c. notation"
            " position across a list of specimens."
        )
    )
    parser.add_argument(
        "--csv",
        required=True,
        help=(
            "Path to input CSV containing a 'Specimen Identifier' column,"
            " e.g. cf_cases.csv"
        ),
    )
    parser.add_argument(
        "--gene",
        required=True,
        help="Gene symbol to search for, e.g. CFTR",
    )
    parser.add_argument(
        "--position",
        required=True,
        help="c. notation prefix to match against CSQ_HGVSc, e.g. 'c.1210-'",
    )
    parser.add_argument(
        "--output",
        default="variant_hits.csv",
        help=(
            "Filename for the output CSV of matching variants"
            " (default: variant_hits.csv)"
        ),
    )
    parser.add_argument(
        "--start",
        default="2023-06-01",
        help=(
            "Start date to look for projects from, e.g. 2023-06-01"
            " (default: 2023-06-01)"
        ),
    )
    return parser.parse_args()


def load_sample_ids(csv_path: str) -> pd.DataFrame:
    """
    Load cases from CSV and derive specimen_id from the
    'Specimen Identifier' column.

    Parameters
    ----------
    csv_path : str
        path to the cases CSV file

    Returns
    -------
    pd.DataFrame
        cases dataframe with an added 'specimen_id' column
    """
    cases = pd.read_csv(csv_path)
    cases["specimen_id"] = cases["Specimen Identifier"].str.strip("SP-")
    return cases


def find_projects(name, start=None):
    """
    Find DNANexus projects by name

    Parameters
    ----------
    project_name : str
        project name
    start : str (optional)
        start date to look for projects from
    end: str (optional)
        end date to look for projects until

    Returns
    -------
    projects : list
        list of DNAnexus projects
    """
    projects = list(
        dxpy.find_projects(
            name=name,
            created_after=start,
            name_mode="regexp",
            describe={"fields": {"name": True}},
        )
    )

    return projects


def determine_genome_build(project_name: str) -> str:
    """
    Determine genome build from a DNAnexus project name.

    Projects processed against GRCh38 have names ending in '38_CEN' or
    '38_TWE'; GRCh37 projects have the same suffix without the '38',
    e.g. '..._CEN' or '..._TWE'.

    Parameters
    ----------
    project_name : str
        DNAnexus project name

    Returns
    -------
    str
        'GRCh38', 'GRCh37', or 'unknown' if the name doesn't match
        either expected suffix pattern
    """
    if not isinstance(project_name, str):
        return "unknown"

    if re.search(r"38_(CEN|TWE)$", project_name):
        return "GRCh38"
    elif re.search(r"_(CEN|TWE)$", project_name):
        return "GRCh37"

    return "unknown"


def find_files_in_project(
    project_id: str, name: str, name_mode: str = "regexp"
) -> list:
    """
    Find files in a project by name.

    Parameters
    ----------
    project_id : str
        The project ID.
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
            project=project_id,
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


def find_all_vcfs(
    sample_ids: list, project_pattern: str, start: str
) -> pd.DataFrame:
    """
    Search DNAnexus for filtered VCFs matching the given sample IDs across
    projects matching project_pattern, and return as a flattened df.

    Parameters
    ----------
    sample_ids : list
        list of specimen IDs to build the file search regex from
    project_pattern : str
        regex pattern to match project names, e.g. "^002.*(CEN|TWE)$"
    start : str
        start date to look for projects from, e.g. "2023-06-01"

    Returns
    -------
    pd.DataFrame
        flattened dataframe of matching files, with columns project_id,
        project_name, genome_build, file_id, name, archivalState,
        created and specimen_id
    """
    ids_pattern = "|".join(re.escape(s) for s in sample_ids)
    file_pattern = rf".*(?:{ids_pattern}).*filter\.vcf\.gz$"
    projects = find_projects(name=project_pattern, start=start)
    project_name_map = {p["id"]: p["describe"]["name"] for p in projects}

    all_vcfs = []
    for project in projects:
        vcfs = find_files_in_project(
            project["id"], name=file_pattern, name_mode="regexp"
        )
        all_vcfs.extend(vcfs)

    vcf_df = pd.json_normalize(all_vcfs)
    vcf_df = vcf_df.rename(
        columns={
            "project": "project_id",
            "id": "file_id",
            "describe.name": "name",
            "describe.archivalState": "archivalState",
            "describe.created": "created",
        }
    )
    vcf_df = vcf_df.drop(columns=["describe.id"])
    vcf_df["specimen_id"] = vcf_df["name"].str.split("-").str[1]
    vcf_df["project_name"] = vcf_df["project_id"].map(project_name_map)
    vcf_df["genome_build"] = vcf_df["project_name"].apply(
        determine_genome_build
    )
    return vcf_df


def dedupe_to_earliest(vcf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only the earliest-created file per specimen_id. Note: this is because
    in one case a reanalysis triggered the whole run to be reprocessed
    and in the other case eggd_dias_batch was run due to an issue with
    CNV calling, so the earliest file is the correct/same one.

    Parameters
    ----------
    vcf_df : pd.DataFrame
        dataframe of files, must contain 'specimen_id' and 'created'
        columns

    Returns
    -------
    pd.DataFrame
        deduplicated dataframe with one row per specimen_id
    """
    sorted_df = vcf_df.sort_values("created")

    # groups with more than one file = the ones that need deduping
    dup_specimens = sorted_df["specimen_id"][
        sorted_df["specimen_id"].duplicated(keep=False)
    ].unique()

    if len(dup_specimens) > 0:
        print(
            f"Found {len(dup_specimens)} specimen_id(s) with multiple files:"
        )
        for specimen_id in dup_specimens:
            group = sorted_df[sorted_df["specimen_id"] == specimen_id]
            kept = group.iloc[0]
            dropped = group.iloc[1:]

            print(f"\nspecimen_id: {specimen_id}")
            print(
                f"  KEPT   : {kept['name']} (file_id: {kept['file_id']},"
                f" created: {kept['created']})"
            )
            for _, row in dropped.iterrows():
                print(
                    f"  DROPPED: {row['name']} (file_id: {row['file_id']},"
                    f" created: {row['created']})"
                )
    else:
        print("No duplicate specimen_ids found.")

    return sorted_df.drop_duplicates(subset="specimen_id", keep="first")


def bulk_unarchive_per_project(df: pd.DataFrame):
    """
    Bulk unarchive files per project.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe of files to unarchive, must contain 'project_id' and
        'file_id' columns

    Returns
    -------
    None
        triggers unarchiving via the DNAnexus API; raises RuntimeError
        if any project's unarchive call fails
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


def download_vcfs(df: pd.DataFrame, out_dir: str = "vcfs") -> pd.DataFrame:
    """
    Download each file in df to out_dir, skipping any that already
    exist locally.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe of files to download, must contain 'file_id',
        'project_id' and 'name' columns
    out_dir : str, optional
        directory to download files into, by default "vcfs"

    Returns
    -------
    pd.DataFrame
        copy of df with an added 'local_path' column pointing to each
        file (whether newly downloaded or already present)
    """
    Path(out_dir).mkdir(exist_ok=True)
    local_paths = []
    skipped = 0
    for _, row in df.iterrows():
        out_path = Path(out_dir) / row["name"]
        if out_path.exists():
            skipped += 1
        else:
            dxpy.bindings.dxfile_functions.download_dxfile(
                row["file_id"], str(out_path), project=row["project_id"]
            )
        local_paths.append(str(out_path))

    if skipped:
        print(
            f"Skipped {skipped} files already downloaded, downloaded"
            f" {len(df) - skipped} new."
        )

    df = df.copy()
    df["local_path"] = local_paths
    return df


def unarchive_and_download(
    final_df: pd.DataFrame, out_dir: str = "vcfs"
) -> pd.DataFrame:
    """
    Trigger unarchiving for any non-live files, then download whatever
    is currently live. Files just triggered for unarchive won't be live
    yet, so this run will only pick up files that were already live —
    you'll need to re-run later to catch the rest.

    Parameters
    ----------
    final_df : pd.DataFrame
        dataframe of files to process, must contain 'project_id',
        'file_id', 'name' and 'archivalState' columns
    out_dir : str, optional
        directory to download live files into, by default "vcfs"

    Returns
    -------
    pd.DataFrame
        subset of final_df corresponding to files that were live and
        downloaded, with an added 'local_path' column. Empty (with a
        'local_path' column) if no files were live.
    """
    non_live = final_df[final_df["archivalState"] != "live"]
    if not non_live.empty:
        bulk_unarchive_per_project(non_live)

    live_vcfs = final_df[final_df["archivalState"] == "live"]
    print(
        f"{len(live_vcfs)} files live and ready to download, {len(non_live)}"
        " still unarchiving."
    )

    if live_vcfs.empty:
        return pd.DataFrame(columns=final_df.columns.tolist() + ["local_path"])
    return download_vcfs(live_vcfs, out_dir=out_dir)


def index_vcfs(download_vcf: pd.DataFrame):
    """
    Run bcftools index on each downloaded VCF, removing any existing
    index first so a fresh one is always created with an up-to-date
    timestamp.

    Parameters
    ----------
    download_vcf : pd.DataFrame
        dataframe of downloaded files, must contain 'local_path'

    Returns
    -------
    None
        indexes each VCF in place (creates a .tbi alongside it);
        raises RuntimeError if indexing fails for any file
    """
    for local_path in download_vcf["local_path"]:
        tbi_path = Path(f"{local_path}.tbi")
        csi_path = Path(f"{local_path}.csi")
        for index_path in (tbi_path, csi_path):
            if index_path.exists():
                index_path.unlink()

        try:
            subprocess.run(
                ["bcftools", "index", "--tbi", local_path],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as error:
            print(f"Error indexing {local_path}: {error.stderr}")
            raise RuntimeError(f"Failed to index {local_path}") from error


def find_gene_variants(
    vcf_path: str, gene: str, c_notation_prefix: str
) -> list:
    """
    Parse VCF records and return variants where CSQ_SYMBOL matches the
    given gene and c. notation (from CSQ_HGVSc) starts with the given
    prefix.

    Parameters
    ----------
    vcf_path : str
        path to the VCF file
    gene : str
        gene symbol to match against CSQ_SYMBOL, e.g. "CFTR"
    c_notation_prefix : str
        c. notation prefix to match, e.g. "c.1210-277_1210-276"

    Returns
    -------
    list of dict
        one dict per matching record, with chrom/pos/ref/alt, gene,
        transcript, hgvsc, filter, sample and gt
    """
    matches = []
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            record_gene = record.info.get("CSQ_SYMBOL")
            hgvsc = record.info.get("CSQ_HGVSc")

            if isinstance(record_gene, tuple):
                record_gene = record_gene[0] if record_gene else None
            if isinstance(hgvsc, tuple):
                hgvsc = hgvsc[0] if hgvsc else None

            if record_gene != gene or not hgvsc:
                continue

            transcript, _, c_dot = hgvsc.partition(":")
            if not c_dot.startswith(c_notation_prefix):
                continue

            # FILTER: empty keys() means PASS, otherwise list of failed filter names
            filter_val = (
                ";".join(record.filter.keys())
                if record.filter.keys()
                else "PASS"
            )

            # single-sample VCF: grab the one sample's GT directly
            sample_name = next(iter(record.samples.keys()))
            sample_data = record.samples[sample_name]
            gt = sample_data.get("GT")
            if gt is not None:
                sep = "|" if sample_data.phased else "/"
                gt_value = sep.join(
                    str(a) if a is not None else "." for a in gt
                )
            else:
                gt_value = "./."

            matches.append(
                {
                    "chrom": record.chrom,
                    "pos": record.pos,
                    "ref": record.ref,
                    "alt": ",".join(record.alts) if record.alts else "",
                    "gene": record_gene,
                    "transcript": transcript,
                    "c_dot": c_dot,
                    "filter": filter_val,
                    "sample": (
                        sample_name.split("-")[1]
                        if "-" in sample_name
                        else sample_name
                    ),
                    "gt": gt_value,
                }
            )
    return matches


def run_gene_search(
    download_vcf: pd.DataFrame, gene: str, c_notation_prefix: str
) -> pd.DataFrame:
    """
    Run find_gene_variants across all downloaded VCFs and return the
    combined results.

    Parameters
    ----------
    download_vcf : pd.DataFrame
        dataframe of downloaded files, must contain 'local_path',
        'file_id', 'name', 'specimen_id' and 'genome_build' columns
    gene : str
        gene symbol to search for, e.g. "CFTR"
    c_notation_prefix : str
        c. notation prefix to match, e.g. "c.1210-"

    Returns
    -------
    pd.DataFrame
        one row per matching variant, with file_id, sample_file,
        specimen_id and genome_build added to the fields returned by
        find_gene_variants
    """
    all_matches = []
    for _, row in download_vcf.iterrows():
        matches = find_gene_variants(
            row["local_path"], gene, c_notation_prefix
        )
        for m in matches:
            m["file_id"] = row["file_id"]
            m["sample_file"] = row["name"]
            m["specimen_id"] = row["specimen_id"]
            m["genome_build"] = row["genome_build"]
            all_matches.append(m)
    return pd.DataFrame(all_matches)


def report_missing_samples(
    all_sample_ids: list,
    vcf_df: pd.DataFrame,
    download_vcf: pd.DataFrame,
    results_df: pd.DataFrame,
):
    """
    Print specimen_ids at each stage where they dropped out: not found
    on DNAnexus at all, found but not downloaded (e.g. still archived),
    or downloaded/searched but with no matching variant.

    Parameters
    ----------
    all_sample_ids : list
        full list of specimen_ids expected (e.g. from the input csv)
    vcf_df : pd.DataFrame
        dataframe of all files found on DNAnexus (before dedup/download),
        must contain 'specimen_id'
    download_vcf : pd.DataFrame
        dataframe of files actually downloaded and searched, must
        contain 'specimen_id'
    results_df : pd.DataFrame
        dataframe of matching variants, must contain 'specimen_id' if
        non-empty

    Returns
    -------
    None
        prints a breakdown of sample counts and lists of specimen_ids
        at each drop-off point
    """
    expected = set(all_sample_ids)
    found_on_dx = set(vcf_df["specimen_id"])
    searched = set(download_vcf["specimen_id"])
    hit = set(results_df["specimen_id"]) if not results_df.empty else set()

    not_found = sorted(expected - found_on_dx)
    found_not_searched = sorted(found_on_dx - searched)
    searched_no_hit = sorted(searched - hit)

    print(f"Expected {len(expected)} samples.")
    print(
        f"{len(found_on_dx)} found on DNAnexus, {len(not_found)} NOT found at"
        " all."
    )
    print(
        f"{len(searched)} downloaded and searched, {len(found_not_searched)}"
        " found but not searched (e.g. still archived)."
    )
    print(
        f"{len(hit)} had matching variants, {len(searched_no_hit)}"
        " searched with no hit."
    )

    if not_found:
        print("\nSamples not found on DNAnexus:")
        for s in not_found:
            print(f"  {s}")
    if found_not_searched:
        print("\nSamples found but not searched (check archivalState):")
        for s in found_not_searched:
            print(f"  {s}")
    if searched_no_hit:
        print("\nSamples searched with no matching variant:")
        for s in searched_no_hit:
            print(f"  {s}")


def main():
    args = parse_args()

    cases = load_sample_ids(args.csv)
    sample_ids = cases["specimen_id"].tolist()

    vcf_df = find_all_vcfs(
        sample_ids, project_pattern="^002.*(CEN|TWE)$", start=args.start
    )
    final_df = dedupe_to_earliest(vcf_df)
    final_df.to_csv("vcfs_to_process.csv", index=False)

    download_vcf = unarchive_and_download(final_df, out_dir="vcfs")
    if download_vcf.empty:
        print("No files downloaded — nothing to search. Exiting.")
        sys.exit()

    index_vcfs(download_vcf)

    results_df = run_gene_search(download_vcf, args.gene, args.position)
    results_df.to_csv(args.output, index=False)
    print(f"\nWrote {len(results_df)} variant hits to {args.output}")

    report_missing_samples(
        all_sample_ids=sample_ids,
        vcf_df=vcf_df,
        download_vcf=download_vcf,
        results_df=results_df,
    )


if __name__ == "__main__":
    main()
