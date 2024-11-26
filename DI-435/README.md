## Find and merge VCFs for creation of a 38 POP AF VCF

### Description
This folder holds code to create an internal GRCh38 POP AF VCF [EBH-435](https://cuhbioinformatics.atlassian.net/browse/DI-435). There are 2 modes implemented:
- `find_qc` searches for GRCh38 projects and finds the related GRCh37 project to retrieve the QC status file (used to remove failed samples) and keep only the first instance of any duplicated samples.
- `no_qc` does not search for QC status files (keeps all files) but removes any duplicated samples completely.

### Python script to find VCFs
#### Usage
Example command to find all CEN GRCh38 VCFs for all runs:
```
python3 find_vcfs_to_merge.py \
    --assay "*CEN38" \
    --outfile_prefix CEN38 \
    --run_mode find_qc
```
Example command to find TWE GRCh38 VCFs where the GRCh37 runs were between two dates:
```
python3 find_vcfs_to_merge.py \
    --assay "*TWE38" \
    --outfile_prefix TWE38 \
    --run_mode find_qc \
    --start 2024-03-15 \
    --end 2024-08-04
```

Example command to find Medicover VCFs with no QC status files:
```
python3 find_vcfs_to_merge.py \
    --assay "*TWE38M" \
    --outfile_prefix TWE38M \
    --run_mode no_qc
```

#### Inputs
- `-a --assay`: The GRCh38 project search term for in DNAnexus, e.g. `"*CEN38"`.
- `-o --outfile_prefix`: Prefix to use to name the output files, e.g. `CEN38`.
- `-r --run_mode`: A choice of whether to find and use QC status files (`find_qc`) or not (`no_qc`).
In `find_qc` mode the following inputs are allowed:
  - `-s --start (optional)`: A date used to find DNAnexus (GRCh37) projects created after, see searching dates section below.
  - `-e --end (optional)`: A date used to find DNAnexus (GRCh37) projects created before, see searching dates section below.

**Searching Dates**
These dates restrict the projects collated to only GRCh38 projects which have corresponding GRCh37 projects which were created within the specified dates.

#### find_qc mode
How the script works in `find_qc` mode:
1. Finds all DNAnexus projects with suffix `--assay`.
2. Finds the related GRCh37 project (and between `--start` and `--end` dates, if provided) for each of these projects and reads in all of the QC status files into one merged dataframe. If multiple QC status files exist in a project, the one created last is used.
3. Finds all raw VCFs in each of the DNAnexus projects.
4. Splits these VCFs into a list of validation (including control) samples and non-validation samples based on naming conventions
5. Removes the first instance of a sample which is duplicated or any which failed QC at any time based on information within the QC status files.
6. Creates a final list of all VCFs to merge and writes this out to file with `--outfile_prefix`.
##### Outputs
- A TSV listing the VCF files for all non-validation samples to merge (named by `{outfile_prefix}_files_to_merge.txt`)
- A CSV of all validation samples found (`{outfile_prefix}_validation_samples.csv`)
- A CSV of all projects within search but missing QC file in DNAnexus and therefore not included (`{outfile_prefix}_projects_missing_QC.csv`).


#### no_qc mode
How the script works in `no_qc` mode:
1. Finds all DNAnexus projects with suffix `--assay`
2. Finds all raw VCFs in each of the DNAnexus projects.
3. Removes all samples which are duplicates.
4. Creates a final list of all VCFs to merge and writes this to file with `--outfile_prefix`.
##### Outputs
- A TSV listing the VCF files to merge (named by `{outfile_prefix}_files_to_merge.txt`)
- A TSV of all duplicate files which have been removed (`{outfile_prefix}_all_dup_rows.txt`)

### Bash script to merge VCFs
The bash script is run in a DNAnexus cloud workstation and requires positional inputs of:
- The output file generated from the Python script above
- The reference genome for GRCh38

Example bash script command:
```
bash merge_VCF_AF.sh \
    CEN38_files_to_merge.txt \
    GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18_noChr.fasta
```

Output:
Merged VCF file and index named from the input files_to_merge txt file, i.e. `final_merged_CEN38_files_to_merge.vcf.gz` and `final_merged_CEN38_files_to_merge.vcf.gz.tbi` which are both uploaded to the DNAnexus project the cloud workstation is running within.