## Compare VEP annotations between versions

### query_validator_api.py
#### Inputs
- `--input_file (str)`: Full path to the tab-separated input file with mismatches in Consequence, HGVSc or HGVSp annotations between VEP versions.
- `--genome_build (str)`: Genome build the samples were run in (GRCh37 or GRCh38)
- `--assay (str)`: Specifies the assay which results are being compared.
- `--chunks (int)`: Specifies how many variants should be in each chunk when querying the API.
- `--output_dir (str)`: Full path to the directory where outputs should be saved.

#### Running
Example command:
```
python query_validator_api.py \
    --input_file MYE_any_mismatches.tsv \
    --genome_build GRCh38 \
    --assay MYE \
    --chunks 500 \
    --output_dir /home/vep_113/bulk_comparison/MYE
```

#### Output
- Files with the VariantValidator results for each chunk of the variants, named `validator_{assay}_chunk{n}.tsv`. The number of files depends on how many variants you have inputted and how large the chunks are that you have specified.
- A final file with all rows from the initial mismatches file inputted, plus columns with output from VariantValidator `HGVSc_validator` and `HGVSp_validator`
