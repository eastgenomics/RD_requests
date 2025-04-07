import argparse

import dxpy


def get_vcfs(vcf_folder):
    """Find the vcfs for the given DNAnexus folder path

    Parameters
    ----------
    vcf_folder : str
        String of the DNAnexus path to find the vcfs in

    Returns
    -------
    list
        List of the file ids for the found vcfs
    """

    vcfs = []

    for vcf in dxpy.find_data_objects(
        name="*_markdup_recalibrated_Haplotyper.vcf.gz",
        name_mode="glob",
        folder=vcf_folder,
        recurse=True,
    ):
        vcfs.append(vcf["id"])

    return vcfs


def match_vcfs(truth_vcfs, query_vcfs):
    """Get the sample name from the truth vcf and create a dict containing the
    input names for the happy app and the file ids for the vcfs

    Parameters
    ----------
    truth_vcfs : list
        List of file ids for the truth vcfs
    query_vcfs : list
        List of file ids for the query vcfs

    Returns
    -------
    dict
        Dict containing the sample name, the input names for the happy app and
        the appropriate file ids
    """

    matched_vcfs = {}

    for truth_vcf in truth_vcfs:
        truth_sample_id = dxpy.describe(truth_vcf)["name"].split("_")[0]

        for query_vcf in query_vcfs:
            query_sample_id = dxpy.describe(query_vcf)["name"].split("_")[0]

            if truth_sample_id == query_sample_id:
                matched_vcfs[truth_sample_id] = {
                    "truth_vcf": dxpy.dxlink(truth_vcf),
                    "query_vcf": [dxpy.dxlink(query_vcf)],
                }

    return matched_vcfs


def run_happy(vcfs, folder_path, alias, type_assay):
    """Run the happy app

    Parameters
    ----------
    vcfs : list
        List containing vcf file ids
    folder_path : str
        DNAnexus folder path where you want to store job output, start with /
    alias : str
        Sample name to annotate the job name with
    type_assay: str
        String representing the type of assay which will decide some of the
        inputs used

    Returns
    -------
    str
        String representing the job id
    """

    # Set the job input
    if type_assay == "TWE":
        job_input = vcfs | {
            "panel_bed": dxpy.dxlink("file-G2V8k90433GVQ7v07gfj0ggX"),
            "high_conf_bed": dxpy.dxlink("file-FjkzKPQ4yBGgBgJ6KvyZ563P"),
            "reference_fasta_tar": dxpy.dxlink(
                "file-F3zxG0Q4fXX9YFjP1v5jK9jf"
            ),
            "reference_sdf_tar": dxpy.dxlink("file-Fjp4k3j4yBGpK4x73bzbG2P0"),
            "happy_docker": dxpy.dxlink("file-GFGbK48433GzV4y54b25p43Z"),
            "reppy_docker": dxpy.dxlink("file-GFGbK48433Gk4xYG8KK05QqY"),
        }
    elif type_assay == "CEN":
        job_input = vcfs | {
            "panel_bed": dxpy.dxlink("file-G620390433GYGY34Jq6Zq1Xf"),
            "high_conf_bed": dxpy.dxlink("file-FjkzKPQ4yBGgBgJ6KvyZ563P"),
            "reference_fasta_tar": dxpy.dxlink(
                "file-F3zxG0Q4fXX9YFjP1v5jK9jf"
            ),
            "reference_sdf_tar": dxpy.dxlink("file-Fjp4k3j4yBGpK4x73bzbG2P0"),
            "happy_docker": dxpy.dxlink("file-GFGbK48433GzV4y54b25p43Z"),
            "reppy_docker": dxpy.dxlink("file-GFGbK48433Gk4xYG8KK05QqY"),
        }

    # Run the job
    job_run = dxpy.DXApp(dxid="app-GgXpYg84G2yXVqz0j10Q7qy5").run(
        job_input,
        folder=f"{folder_path}",
        name=f"happy_{alias}",
    )

    # Return the job_id
    return job_run.describe(fields={"id": True})["id"]


def main(truth_vcf_folder, query_vcf_folder, assay, destination):
    # set standard variables
    standard_path = destination

    truth_vcfs = get_vcfs(truth_vcf_folder)
    query_vcfs = get_vcfs(query_vcf_folder)

    happy_jobs = []

    matched_vcfs = match_vcfs(truth_vcfs, query_vcfs)

    for sample_name, vcf_dict in matched_vcfs.items():
        job_id = run_happy(
            vcf_dict,
            alias=sample_name,
            folder_path=standard_path,
            type_assay=assay,
        )
        happy_jobs.append(job_id)

    # Print me out the jobs that I want
    print("Here are the list of jobs for happy:")
    print(happy_jobs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("truth_vcf_folder", help="Truth vcf folder")
    parser.add_argument("query_vcf_folder", help="Query vcf folder")
    parser.add_argument("assay", choices=["TWE", "CEN"], help="Type of assay")
    parser.add_argument(
        "-d",
        "--destination",
        required=True,
        help="Destination of the jobs",
    )

    args = parser.parse_args()
    main(
        args.truth_vcf_folder,
        args.query_vcf_folder,
        args.assay,
        args.destination,
    )
