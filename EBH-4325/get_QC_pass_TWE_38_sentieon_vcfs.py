import dxpy
import pandas as pd


def get_projects(
    search_term,
    number_of_projects=None,
    after_date=None,
    before_date=None,
    search_mode="regexp",
):
    """
    Find projects within specified data range using specified search term and
        mode.

    Parameters
    ----------
    search_term : str
        Search term or regexp pattern used to find projects.
    number_of_projects : int (optional)
        Number of projects to use. If set, takes most recent projects based
        on the date in the project name rather than creation date.
        Defaults to None.
    after_date : str (optional)
        Select projects created after this date.
    before_date : str (optional)
        Select projects created before this date. Defaults to None.
    search_mode : str (optional)
        Type of dxpy search mode to use, acceptable search_modes "regexp",
        "glob" and "exact". Defaults to "regexp".

    Returns
    -------
    projects : list
        list of dictionaries containing information (project ID/name)
        on the selected projects
    """
    projects = list(
        dxpy.bindings.search.find_projects(
            name=search_term,
            name_mode=search_mode,
            created_after=after_date,
            created_before=before_date,
            describe={"fields": {"name": True}},
        )
    )

    projects = sorted(projects, key=lambda x: x["describe"]["name"])

    if number_of_projects is not None:
        projects = projects[-int(number_of_projects):]

    return projects


def find_files(filename_pattern, project_id, name_mode="regexp", folder=None):
    """
    Find files in a project using the specified search term and mode.

    Parameters
    ----------
    filename_pattern : str
        Search term or regexp pattern used to find files.
    project_id : str
        DNAnexus project ID
    name_mode : str
        Type of dxpy search mode to use, acceptable search
        modes "regexp", "glob" and "exact"
    folder : str (optional)
        Folder path to search within the project. Defaults to None.

    Returns
    -------
    files_found : list
        list of dictionaries containing info (file ID/name) about
        the selected files

    Raises
    ------
    AssertionError
        Raised if no files are found for the specified search term
    """
    files_found = list(
        dxpy.bindings.search.find_data_objects(
            classname="file",
            name=filename_pattern,
            name_mode=name_mode,
            project=project_id,
            folder=folder,
            describe={"fields": {"name": True}},
        )
    )

    assert (
        len(files_found) > 0
    ), f"No files found for {filename_pattern} in {project_id}"

    return files_found


def read2df(
    file_id: str,
    project: dict,
    separator,
    mode,
    skiprows=None
):
    """
    Read in file to pandas df

    Parameters
    ----------
    file_id : str
        DNAnexus file ID of tsv file to be read in.
    project : dict
        DNAnexus project ID of project containing file to be read in.
    separator : str
        Separator used within the file
    mode : str
        the mode to open the file with using dxpy
    skiprows : int, optional
        number of rows to skip at the start of the file, default is None

    Returns:
    df : pd.DataFrame
        pd.DataFrame object of file
    """
    file = dxpy.open_dxfile(file_id, project=project, mode=mode)
    df = pd.read_csv(file, sep=separator, skiprows=skiprows)
    return df


DESTINATION_PROJ_ID = "project-J2KFZPj4Vg5Kkj4VFyG55zXK"
DESTINATION_FOLDER_PATH = "/exome_grch38"
N_VCFS = 100


def main():

    # Only collect lastest 38_TWE projects created after 12th August 2025
    twe_38_projects_ids = [
        proj["id"] for proj in get_projects(
            search_term="^002.*38_TWE$",
            after_date="2025-08-12")
    ]

    vcf_file_details = []
    for proj in twe_38_projects_ids:
        # Stop looping through projects once we have N VCF file details
        if len(vcf_file_details) >= N_VCFS:
            break

        # Find and read in manifest file to get list of samples which have
        # passed QC
        manifest = find_files(
            filename_pattern="25-NGSA.*csv$",
            project_id=proj,
            folder="/"
        )

        # Expect only one manifest per project
        assert len(
            manifest) == 1, f"Did not find exactly one manifest in {proj}"

        # Skip the first two lines of the manifest as we do not need them,
        # use ";" as the separator
        manifest_df = read2df(
            file_id=manifest[0]["id"],
            project=proj,
            separator=";",
            mode="r",
            skiprows=2
        )

        # Join the values in the 3rd (specimen ID) and 4th (instrument id)
        # columns with an "-" to get the sample ID
        sample_ids = manifest_df.apply(
            lambda row: f"{row.iloc[3]}-{row.iloc[2]}", axis=1).tolist()

        for sample in sample_ids:

            # Find the Sentieon VCF and index files for the sample
            sentieon_vcf = find_files(
                filename_pattern=f"^{sample}.*Haplotyper.vcf.gz",
                project_id=proj
            )

            # Expect only two files, the VCF and its index
            assert len(sentieon_vcf) == 2, f"Did not find exactly two Sentieon VCF files (vcf and index) for sample {sample} in project {proj}"

            vcf_id, vcf_index_id = None, None
            for file in sentieon_vcf:
                if file["describe"]["name"].endswith(".vcf.gz"):
                    vcf_id = file["id"]

                if file["describe"]["name"].endswith(".vcf.gz.tbi"):
                    vcf_index_id = file["id"]

            if vcf_id is None or vcf_index_id is None:
                raise ValueError(f"Could not determine VCF and index IDs for sample {sample} in project {proj}")

            vcf_file_details.append((sample, proj, vcf_id, vcf_index_id))

            # Stop adding VCF file details once we have N entries
            if len(vcf_file_details) >= N_VCFS:
                break

    # Output summary of files to be cloned in a csv
    output_df = pd.DataFrame(vcf_file_details, columns=["Sample_ID", "Project_ID", "VCF_ID", "VCF_Index_ID"])
    output_df.to_csv("vcfs.csv", index=False)

    # Clone the VCF and index files to the destination project
    for proj_id, vcf_id, vcf_index_id in zip(
        output_df["Project_ID"],
        output_df["VCF_ID"],
        output_df["VCF_Index_ID"]
    ):
        # Clone the VCF
        dxpy.bindings.DXFile(dxid=vcf_id, project=proj_id).clone(
            project=DESTINATION_PROJ_ID,
            folder=DESTINATION_FOLDER_PATH
        )

        # Clone the VCF index
        dxpy.bindings.DXFile(dxid=vcf_index_id, project=proj_id).clone(
            project=DESTINATION_PROJ_ID,
            folder=DESTINATION_FOLDER_PATH
        )


if __name__ == "__main__":
    main()
