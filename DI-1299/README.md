# Compare Somalier outputs in GRCh37 and GRCh38

### Description
The `compare_b37_and_b38.py` script is used to compare Somalier Predicted_Sex values for samples in GRCh37 and GRCh38 and return any mismatches.

### Inputs
- `--config (str)`: Path to a config file with relevant variables for the assay.

Example config:
```json
{
    "assay": "TWE",
    "search_term": "*TWE38",
    "number_of_projects": 10,
    "filename": "Multiqc_somalier.samples.tsv",
    "column_to_compare": "Predicted_Sex",
    "sample_column": "sample_id",
    "variables_to_plot": [
        [
            "X_n", "X_het", "X_hom_ref", "X_hom_alt", "X_depth_mean"
        ],
        [
            "Y_n", "Y_depth_mean"
        ],
        [
            "n_het", "n_hom_ref", "n_unknown", "n_hom_alt"
        ],
        [
            "depth_mean", "depth_sd", "gt_depth_mean", "gt_depth_sd",
            "ab_mean", "ab_std", "p_middling_ab"
        ]
    ]
}
```

### Running
Example command:
```
python compare_b37_and_b38.py --config TWE_config.json
```

### How it works
The script:
- Find all projects using the `search_term` in GRCh38 in DNAnexus.
- Gets the related GRCh37 project for each GRCh38 project with the prefix `002_` and the suffix `_{assay}`.
 - Within the GRCh37 and GRCh38 projects, find all files by `filename` and reads them all in with pandas.
- Identifies any mismatches in `column_to_compare` between the GRCh37 and GRCh38 values for each sample and investigate.
- If mismatches are found, plots the variables in `variables_to_plot` to look at differences in these values for each sample between genome builds.

### Output
- `{assay}_all_results.tsv`: A TSV with Somalier results for all samples found; headers have '_GRCh37' and '_GRCh38' suffixes to show which genome the result is from.
- `{assay}_all_mismatches.tsv`: A TSV with all rows where there is a mismatch between the `column_to_compare` value in GRCh37 and GRCh38.
- If there are samples with mismatches found, one plot is created for all metrics included in each nested list within `variables_to_plot` for each sample in GRCh37 and GRCh38. Each variable in each nested list is plotted as a scatter plot, one per row, and each plot is saved to HTML.