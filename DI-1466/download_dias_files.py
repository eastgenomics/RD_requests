"""
Given an eggd_dias_batch job ID, this script downloads the associated SNV/CNV
  reports if they contain variants and prints a short summary describing the
  number of SNV and CNV reports downloaded. The script also downloads the
  artemis file, QC status file and MultiQC report - if these are able to be
  found.
"""

import concurrent
import argparse
import dxpy


def parse_args():
    """
    Parse arguments given at cmd line.

        Args: None

        Returns:
            - args (Namespace): object containing parsed arguments.
    """

    parser = argparse.ArgumentParser(
        description="Download SNV/CNV reports, artemis file, QC status file\
            and Multi QC file"
        )

    parser.add_argument(
        "-b", "--batch_job_id",
        help="eggd_dias_batch job ID that was used to create the reports you\
            wish to download", required=True
        )

    args = parser.parse_args()

    return args


def read_batch_job_metadata(batch_job_id):
    """
    Using a eggd_dias_batch job ID return a dict containing: the project ID for
     the project in which eggd_dias_batch was run, a list of the launched jobs,
     the QC status file ID and multiqc report file ID (when possible). If the
     QC status and/or MultiQC report file IDs cannot be retrieved from the
     batch job metadata, then a message is printed to say they were not found.
     If no artemis job was launched by dias batch, a message is printed to
     alert the user.

    Parameters
    ----------
    batch_job_id : str
        eggd_dias_batch job ID

    Returns
    -------
    dx_ids : dict
        Dictionary object containing: the project ID for the project in which
         eggd_dias_batch was run, a list of the launched jobs, the QC status
         and multiqc report file IDs (when possible).
    """

    dx_ids = {}

    dx_describe = dxpy.bindings.DXJob(dxid=batch_job_id).describe(
        fields={"output": True, "input": True, "project": True}
    )

    dx_ids["batch_project_id"] = dx_describe["project"]

    dx_ids["launched_jobs_list"] = dx_describe["output"][
        "launched_jobs"].split(',')

    # QC status, multiqc report and artemis are optional input for batch so
    # may not be present
    qc_file_describe = dx_describe["input"].get("qc_file")
    if qc_file_describe is not None:
        dx_ids["qc_file"] = qc_file_describe["$dnanexus_link"]
    else:
        print(
            "QC status file could not be found in eggd_dias_batch job metadata"
        )

    multiqc_report_describe = dx_describe["input"].get("multiqc_report")
    if multiqc_report_describe is not None:
        dx_ids["multiqc_report"] = multiqc_report_describe["$dnanexus_link"]
    else:
        print(
            "MultiQC report could not be found in eggd_dias_batch job metadata"
        )

    artemis_describe = dx_describe["input"].get("artemis")
    if artemis_describe is None or artemis_describe is False:
        print(
            "No artemis job associated with the provided eggd_dias_batch job"
        )

    return dx_ids


def call_in_parallel(func, items, ignore_missing=True,
                     ignore_all_errors=True, **kwargs) -> list:
    """
    Calls the given function in parallel using concurrent.futures on
     the given set of items (i.e for calling dxpy.describe() on multiple
     object IDs).

    Additional arguments specified to kwargs are directly passed to the
     specified function.

    Parameters
    ----------
    func : callable
        function to call on each item
    items : list
        list of items to call function on
    ignore_missing : bool
        controls if to just print a warning instead of raising an
         exception on a dxpy.exceptions.ResourceNotFound being raised.
         This is most likely from a file that has been deleted and we are
         just going to default to ignoring these
    ignore_all_errors : bool
        controls if to just print a warning instead of raising an exception.

    Returns
    -------
    list
        list of responses
    """
    results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        concurrent_jobs = {
            executor.submit(func, item, **kwargs): item for item in items
        }

        for future in concurrent.futures.as_completed(concurrent_jobs):
            # access returned output as each is returned in any order
            try:
                results.append(future.result())
            except Exception as exc:
                if (
                    ignore_missing and
                    isinstance(exc, dxpy.exceptions.ResourceNotFound)
                ):
                    # dx object does not exist and specifying to skip,
                    # just print warning and continue'
                    print(
                        f'WARNING: {concurrent_jobs[future]} could not be '
                        'found, skipping to not raise an exception'
                    )
                    continue
                # catch any other errors that might get raised during querying
                print(
                    f"Warning: Error getting data for {concurrent_jobs[future]}: {exc}"
                )

                if ignore_all_errors is False:
                    raise exc

    return results


def get_file_ids(launched_jobs_list):
    """
    Loop over launched jobs/analysis IDs and retrieve the SNV/CNV .xlsx report
     IDs and artemis file ID (if present).

    Parameters
    ----------
    launched_jobs_list : list
        List of launched jobs/analyses derived from the launched_jobs output
         of dias batch.

    Returns
    -------
    launched_jobs_dict : dict
        Dictionary object containing SNV/CNV .xlsx report IDs and artemis
         file ID (if present).
    """
    launched_jobs_dict = {
        "snv_report_ids": [],
        "cnv_report_ids": [],
    }

    describe_dicts = call_in_parallel(
        dxpy.describe, launched_jobs_list, fields={"executableName": True,
                                                   "output": True}
    )

    failed_report_jobs = []

    for desc in describe_dicts:
        try:
            if desc["executableName"].startswith("dias_reports"):
                launched_jobs_dict["snv_report_ids"].append(
                    desc["output"]["stage-rpt_generate_workbook.xlsx_report"][
                        "$dnanexus_link"]
                )

            elif desc["executableName"].startswith("dias_cnvreports"):
                launched_jobs_dict["cnv_report_ids"].append(
                    desc["output"]["stage-cnv_generate_workbook.xlsx_report"][
                        "$dnanexus_link"]
                )

            elif desc["executableName"] == "eggd_artemis":
                # Should only ever be one or zero artemis job IDs in launched
                # jobs therefore do not need to append to list as done for
                # reports
                launched_jobs_dict["artemis_file_id"] = dxpy.describe(
                    desc["id"], fields={"output": True}
                )["output"]["url_file"]

        except KeyError:
            failed_report_jobs.append(desc["id"])

    if len(failed_report_jobs) > 0:
        print(
            "Warning: output files could not be gathered for the following "
            f"jobs:\n {failed_report_jobs}"
        )

    return launched_jobs_dict


def get_details(file_id, project_id):
    """
    Get details of a file in DNAnexus

    Parameters
    ----------
    file_id : str
        DNAnexus file ID.
    project_id : str
        DNANexus project ID.

    Returns
    -------
    tuple
        Tuple containing (DNAnexus file ID, corresponding file details).
    """
    return file_id, dxpy.DXFile(dxid=file_id, project=project_id).get_details()


def organise_report_files(reports_details, report_type):
    """
    Print the number of SNV/CNV reports which: contain variants, do not contain
     variants or have no associated file details (if such reports exist) and
     return a list of file IDs for reports which are to be downloaded.

    Parameters
    ----------
    reports_details : list
        List of tuples containing DNAnexus file IDs for the SNV/CNV .xlsx
         reports and their corresponding details.
    report_type : str
        String either "SNV" or "CNV", denoting the type of reports.

    Returns
    -------
    reports_for_download : list
        List containing DNAnexus file IDs for the the SNV/CNV .xlsx
         reports which are to be downloaded (i.e. contain variants)
    """

    details_key = {
        "SNV": "included",
        "CNV": "variants"
        }

    no_of_reports_without_details = 0
    no_of_reports_without_vars = 0
    reports_for_download = []

    for report_id, details in reports_details:
        # If old report and no "file details" metadata is available then
        # no_of_vars will = None
        no_of_vars = details.get(details_key[report_type])
        if no_of_vars is None:
            no_of_reports_without_details += 1
        elif no_of_vars == 0:
            no_of_reports_without_vars += 1
        elif no_of_vars > 0:
            reports_for_download.append(report_id)
        else:
            print(
                "Warning: DNAnexus 'file details' metadata "
                f"for {report_id} are not interpretable"
            )

    print(
        f'{len(reports_for_download)} {report_type} report(s) with variants'
    )
    print(
        f'{no_of_reports_without_vars} {report_type} report(s) without'
        ' variants'
    )
    # Unlikely to have reports without details unless using an old batch
    # job, therefore do not need to print the number of reports without details
    # every time
    if no_of_reports_without_details > 0:
        print(
            f'{no_of_reports_without_details} {report_type} report(s) do not '
            'have DNAnexus metadata about number of variants in the report'
        )

    return reports_for_download


def download_single_file(dxid, project):
    """
    Given a single dx file ID, download it with the original filename

    Parameters
    ----------
    dxid : str
        file ID of object to download
    project : str
        project containing the file
    """
    dxpy.bindings.dxfile_functions.download_dxfile(
        dxid,
        dxpy.describe(dxid).get('name'),
        project=project
    )


def main():
    args = parse_args()
    dx_ids = read_batch_job_metadata(args.batch_job_id)
    project_id = dx_ids["batch_project_id"]

    launched_jobs_dict = get_file_ids(dx_ids["launched_jobs_list"])

    snv_reports_details = call_in_parallel(
     get_details, launched_jobs_dict["snv_report_ids"], project_id=project_id
    )
    snv_reports_for_download = organise_report_files(
        snv_reports_details, report_type="SNV"
    )

    # Gather file IDs common to both CEN and TWE
    file_ids_for_download = snv_reports_for_download
    file_ids_for_download.append(dx_ids.get("qc_file"))
    file_ids_for_download.append(dx_ids.get("multiqc_report"))
    file_ids_for_download.append(launched_jobs_dict.get("artemis_file_id"))

    # TWE may have no CNV reports workflows in launched jobs, therefore no
    # CNV report file IDs would be retrieved
    if len(launched_jobs_dict["cnv_report_ids"]) > 0:
        cnv_reports_details = call_in_parallel(
         get_details, launched_jobs_dict["cnv_report_ids"],
         project_id=project_id
        )
        cnv_reports_for_download = organise_report_files(
            cnv_reports_details, report_type="CNV"
        )
        file_ids_for_download += cnv_reports_for_download

    # Remove "None" values from download list where file IDs were not found -
    # user is made  aware of missing files through warnings printed in
    # read_batch_job_metadata()
    file_ids_for_download = [
        file_id for file_id in file_ids_for_download if file_id is not None
    ]

    call_in_parallel(
        download_single_file, file_ids_for_download, project=project_id
    )


if __name__ == "__main__":
    main()
