## DI-1480 - Compare variant filtering
This script compares the total and filtered in variants for a sample between two different conditions, e.g. GRCh37 and GRCh38.

#### Inputs
- `--dx_old (str)`: Path in DNAnexus to where the soft-filtered VCFs are with the old filtering. Format `project-XXX:/folder`
- `--dx_new (str)`: Path in DNAnexus to where the soft-filtered VCFs are with the new filtering. Format `project-XXX:/folder`
- `--old_path (str)`: Full path to local directory where soft-filtered VCFs with the old filtering will be downloaded.
- `--new_path (str)`: Full path to local directory where soft-filtered VCFs with the new filtering will be downloaded.
- `--old_name (str)`: Name of old filter, e.g. `'GRCh37'`
- `--new_name (str)`: Name of new filter, e.g. `'GRCh38'`
- `--fields (str)`: Fields for the PASS variants which we want to include in our final output. E.g. `'CSQ_HGVSc,CSQ_HGVSp,CSQ_SYMBOL,CSQ_Consequence'`. Note: requires CSQ_HGVSc to be included, as variants are merged on this field.
- `--outfile (str)`: Name of output .xlsx file.

#### Running
Example command:
```
python compare_variant_filtering.py \
    --dx_old project-GygXY5044Qqx5fb3jV48Pfkg:/output/CEN-250214_1131/dias_reports_v2.2.2_SNV/250217_1255/eggd_optimised_filtering-1.1.0 \
    --dx_new project-Gz801KQ4qyJGGP0k1KqVpjqK:/output/38_CEN-250303_1143/dias_reports_v2.2.2_SNV/250311_1141/eggd_optimised_filtering-1.1.0 \
    --old_path /home/filtering_comparison/CEN/grch37 \
    --new_path /home/filtering_comparison/CEN/grch38 \
    --old_name GRCh37 \
    --new_name GRCh38 \
    --fields CSQ_HGVSc,CSQ_HGVSp,CSQ_SYMBOL,CSQ_Consequence,CSQ_gnomADe_AF,CSQ_gnomADg_AF,CSQ_TWE_AF,CSQ_ClinVar_CLNSIGCONF,CSQ_ClinVar_CLNSIG,CSQ_HGMD_CLASS,MOI \
    --outfile /home/filtering_comparison/CEN/cen_GRCh37_vs_GRCh38_filtering_comparison.xlsx
```

#### Output
The output .xlsx file is named with `--outfile` which includes one sheet `by_sample` with one row per sample, the sample variant totals (total in panel, total PASS variants and the PASS variants themselves) for both the old and new filtering. The second sheet `by_variants` contains the sample name as the index and one row per variant for each sample; PASS variants with the same HGVS for the same sample present with both old and new filtering are shown in the same row.