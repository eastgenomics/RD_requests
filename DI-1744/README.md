## Create file with paths to XLSX files on clingen

### Description
This script is used to create a file with paths to complete XLSX workbooks on clingen, which can then be processed by the [clinvar submissions process](https://github.com/eastgenomics/clinvar_submissions)

### Inputs
The script `get_xlsx_pathjs.py` takes inputs:
- `--test_code_mapping (str)`: path to csv file with two columns: R code, and CEN or WES depending on which assay the panel is on
- `--clarity_extract (str)`: path to clarity extract

### Running
Example command:
`python3 RD_requests/DI-1744/get_xlsx_paths.py --test_code_mapping /Downloads/mapping.csv --clarity_extract /Downloads/mock_clarity.csv`

### Output
The output of the script is two files:
- paths.csv, a file of just the paths to the XLSX files on clingen
- all_data.csv, a file with all the data collated during the process (useful for troubleshooting)