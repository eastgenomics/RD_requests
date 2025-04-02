import argparse

import dxpy


def parse_file_ids(fastqs):
    data = {}
    data.setdefault("reads_fastqgzs", [])
    data.setdefault("reads2_fastqgzs", [])

    for file_id in fastqs:
        file_data = dxpy.describe(file_id)

        if "_R1_" in file_data["name"]:
            data["reads_fastqgzs"].append(dxpy.dxlink(file_id))

        if "_R2_" in file_data["name"]:
            data["reads2_fastqgzs"].append(dxpy.dxlink(file_id))

    return data


def run_sentieon_fastq_to_vcf(fastqs, folder_path, alias, job_number):
    """
    Run the sentieon_bwa command for me

    Parameters
    ----------
    fastqs : list
        List containing fastq reads file ids
    folder_path : str
        DNAnexus folder path where you want to store job output, start with /
    alias : str
        version number of the sentieon_bwa app

    Returns
    -------
    str
        App will run and the function will return the associated job_id

    """

    # Set the job input
    job_input = fastqs | {
        "genomebwaindex_targz": dxpy.dxlink("file-Gb76f204XGybZ3J6F731xkBp"),
        "genome_fastagz": dxpy.dxlink("file-Gb757784XGyY3FPvkPQ74K9z"),
    }

    # Run the job
    job_run = dxpy.DXApp(dxid="app-Gy4j5z00PPyQ5qv5FBXy0ZZp").run(
        job_input,
        folder=f"{folder_path}/{folder_path}_{job_number}",
        name=f"sentieon-germline-fastq2vcf_v{alias}_{job_number}",
    )

    # Return the job_id
    return job_run.describe(fields={"id": True})["id"]


def main(fastqs, destination):
    # set standard variables
    standard_path = destination

    parsed_fastqs = parse_file_ids(fastqs)

    # run sentieon 10 times and store jobs in a list
    sentieon_jobs = []

    for i in range(1, 11):
        job_id = run_sentieon_fastq_to_vcf(
            parsed_fastqs,
            folder_path=standard_path,
            alias="5.1.0",
            job_number=i,
        )
        sentieon_jobs.append(job_id)

    print("All jobs have started running")

    # Print me out the jobs that I want
    print("Here are the list of jobs:")
    print(sentieon_jobs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file_id", nargs="+", help="Fastqs DNAnexus file ids for one sample"
    )
    parser.add_argument(
        "-d",
        "--destination",
        required=False,
        default="/reproducibility",
        help="Destination of the reproducibility jobs",
    )

    args = parser.parse_args()
    main(args.file_id, args.destination)
