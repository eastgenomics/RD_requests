import argparse

import dxpy


def parse_file_ids(fastqs):
    """Parse the fastq file ids into a dict with the sample id as key

    Parameters
    ----------
    fastqs : list
        List of file ids for the fastqs

    Returns
    -------
    dict
        Dict of dxlink arranged according to the appropriate sample id and
        input name
    """

    data = {}

    for file_id in fastqs:
        if ":" in file_id:
            file_id = file_id.split(":")[-1]

        file_data = dxpy.describe(file_id)

        sample_name = get_sample_name(file_data["name"])

        data.setdefault(sample_name, {})
        data[sample_name].setdefault("reads_fastqgzs", [])
        data[sample_name].setdefault("reads2_fastqgzs", [])

        if "_R1_" in file_data["name"]:
            # have to do this to order the file ids otherwise the sentieon job
            # fails
            if (
                "_L001_" in sample_name
                and data[sample_name]["reads_fastqgzs"] != {}
            ):
                data[sample_name]["reads_fastqgzs"].insert(
                    0, dxpy.dxlink(file_id)
                )
            else:
                data[sample_name]["reads_fastqgzs"].append(
                    dxpy.dxlink(file_id)
                )

        if "_R2_" in file_data["name"]:
            if (
                "_L001_" in sample_name
                and data[sample_name]["reads2_fastqgzs"] != {}
            ):
                data[sample_name]["reads2_fastqgzs"].insert(
                    0, dxpy.dxlink(file_id)
                )
            else:
                data[sample_name]["reads2_fastqgzs"].append(
                    dxpy.dxlink(file_id)
                )

    return data


def get_sample_name(file_name):
    """Get the sample id from the given file name

    Parameters
    ----------
    file_name : str
        File name

    Returns
    -------
    str
        Sample id extracted from the file name
    """

    return file_name.split("_")[0]


def run_sentieon_fastq_to_vcf(fastqs, folder_path, version, sample_id):
    """Run the sentieon fastq2vcf

    Parameters
    ----------
    fastqs : list
        List containing fastq reads file ids
    folder_path : str
        DNAnexus folder path where you want to store job output, start with /
    version : str
        Version of the sentieon app used
    sample_id: str
        Sample id

    Returns
    -------
    str
        Job id for the sentieon job
    """

    # Set the job input
    job_input = fastqs | {
        "genomebwaindex_targz": dxpy.dxlink("file-F404y604F30fbxQG68KF3GZb"),
        "genome_fastagz": dxpy.dxlink("file-F403K904F30y2vpVFqxB9kz7"),
        "gatk_resource_bundle": dxpy.dxlink("file-F3zx7gj4fXX8QG3Q42BzpyZJ"),
        "ignore_decoy": True,
        "output_metrics": True,
    }

    # Run the job
    job_run = dxpy.DXApp(dxid="app-Gy4j5z00PPyQ5qv5FBXy0ZZp").run(
        job_input,
        folder=f"{folder_path}/{sample_id}",
        name=f"sentieon-germline-fastq2vcf_v{version}_{sample_id}",
    )

    # Return the job_id
    return job_run.describe(fields={"id": True})["id"]


def main(fastqs, destination):
    # set standard variables
    standard_path = destination

    sample_files = parse_file_ids(fastqs)

    # run sentieon 10 times and store jobs in a list
    sentieon_jobs = []

    for sample, job_input in sample_files.items():
        job_id = run_sentieon_fastq_to_vcf(
            job_input,
            folder_path=standard_path,
            version="5.1.0",
            sample_id=sample,
        )
        sentieon_jobs.append(job_id)

    print("All jobs have started running")

    # Print me out the jobs that I want
    print("Here are the list of jobs")
    print(sentieon_jobs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file_id", nargs="+", help="Fastqs DNAnexus file ids for one sample"
    )
    parser.add_argument(
        "-d",
        "--destination",
        help="Destination of the job outputs",
    )

    args = parser.parse_args()
    main(args.file_id, args.destination)
