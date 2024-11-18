## Find and merge VCFs for creation of a 38 POP AF VCF

### Description
This folder holds code to create an internal GRCh38 POP AF VCF [EBH-435](https://cuhbioinformatics.atlassian.net/browse/DI-435).

### Python script to find VCFs
The Python script `find_vcfs_to_merge.py` takes inputs:
- `-a --assay`: The project prefix to look for in DNAnexus, e.g. `"*CEN38"`
- `-o --outfile_prefix`: What to name the output TSV file which lists VCF files to merge
- `-r --run_mode`: A choice of whether to find QC status files ('find_qc') or not ('no_qc'). no_qc would be used as an input when searching Medicover projects
- `-s --start (optional)`: A date used to find DNAnexus projects created after
- `-e --end (optional)`: A date used to find DNAnexus projects created before
- `-n --number_of_projects (optional)`: The number of most recent projects to keep

Example Python script command for CEN to find QC status files:
`python3 find_vcfs_to_merge.py --assay "*CEN38" --end 2024-05-03 --outfile_prefix CEN38 --run_mode find_qc`

Example Python script command for Medicover:
`python3 find_vcfs_to_merge.py --assay "*TWE38M"--outfile_prefix TWE38M --run_mode no_qc`

How the script works in `find_qc` mode:
1. Finds all DNAnexus projects with suffix `--assay` and between `--start` and `--end` dates (if provided) and then, if `--number_of_projects` is given, subsets this to the last n runs.
2. Finds the related GRCh37 project for each of these projects and reads in all of the QC status files into one merged dataframe. If multiple QC status files exist in a project, the one created last is used.
3. Finds all raw VCFs in each of the DNAnexus projects.
4. Splits these VCFs into a list of validation (including control) samples and non-validation samples based on naming conventions
5. Removes the first instance of a sample which is duplicated or any which failed QC at any time based on information within the QC status files.
6. Creates a final list of all VCFs to merge and writes this out to file with `--outfile_prefix`.

How the script works in `no_qc` mode:
1. Finds all DNAnexus projects with suffix `--assay` and between `--start` and `--end` dates (if provided) and then, if `--number_of_projects` is given, subsets this to the last n runs.
2. Finds all raw VCFs in each of the DNAnexus projects.
3. Removes all samples which are duplicates.
4. Creates a final list of all VCFs to merge and writes this to file with `--outfile_prefix`.

Output:
- A TSV listing the VCF files for all samples to merge (named by `{outfile_prefix}_files_to_merge.txt`)
- If `--run_mode` is set to `find_qc`, a CSV of all validation samples found (`{outfile_prefix}_validation_samples.csv`)
- If `--run_mode` is set to `no_qc`, a TXT file with all duplicate samples which have been removed (`{outfile_prefix}_all_dup_rows.txt`).

### Bash script to merge VCFs
The bash script is run in a DNAnexus cloud workstation and requires positional inputs of:
- The output file generated from the Python script above
- The job ID for the cloud workstation running
- The reference genome for GRCh38

Example bash script command:
`bash merge_VCF_AF.sh CEN38_vcf_to_merge.txt job-Gjb39Z04bxf82XZ12gPJ2bbV GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18_noChr.fasta`

Output:
Merged VCF file and index named from the job ID given, i.e. `final_merged_job-Gjb39Z04bxf82XZ12gPJ2bbV.vcf.gz` and `final_merged_job-Gjb39Z04bxf82XZ12gPJ2bbV.vcf.gz.tbi` which are both uploaded to the DNAnexus project the cloud workstation is running within.