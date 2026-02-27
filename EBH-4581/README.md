## Query GTEx and SpliceAI data

### Original request
I am supervising a STP project which is looking at VUS’s predicted to affect splicing to see if they could potentially be re-classified using RNA analysis. We would like help with obtaining SpliceAI scores and GTEX isoform expression data (TPM format) for variants/genes in bulk"

### Inputs
The script `query_spliceai_and_gtex.py` takes inputs:
- `--input (str)`: Input Excel file.
- `--tissues (list)`: List of tissues to be used to query GTEx.
- `--output (str)`: The name of the output Excel file.

### Running
Example command:
```
python query_spliceai_and_gtex.py \
    --input variants_file.xlsx \
    --tissues Whole_Blood Cells_Cultured_fibroblasts Cells_EBV-transformed_lymphocytes \
    --output variants_file_plus_gtex_and_spliceai.xlsx \
```
