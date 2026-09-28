#!/bin/bash

# get samples with 1 heterozygous variant
vcf=lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz

# samples with exactly one included het
bcftools query \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\n]' "$vcf" \
  | sort -u > included_hets.tsv

cut -f1 included_hets.tsv | sort | uniq -u > single_het_samples.txt

# rare intronic variants within ±20 of a boundary, any FILTER, in those samples
bcftools query \
    -S single_het_samples.txt \
    -i 'GT="alt"
        & CSQ_HGVSc~"c[.][-*]?[0-9]+[+-]([1-9]|1[0-9]|20)([^0-9]|$)"
        & (CSQ_gnomADg_AF="."      | CSQ_gnomADg_AF<=0.01)
        & (CSQ_gnomADe_AF="."      | CSQ_gnomADe_AF<=0.01)
        & (CSQ_TWE_AF="."          | CSQ_TWE_AF<=0.01)
        & (CSQ_TWE_WES_v1_AF="."   | CSQ_TWE_WES_v1_AF<=0.01)' \
    -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_HGVSc\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\t%CSQ_TWE_WES_v1_AF\n' \
    "$vcf" > single_het_intron_38.tsv


# get heterozygous variants b37

# get samples with 1 heterozygous variant
vcf=lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz

# samples with exactly one included het
bcftools query \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\n]' "$vcf" \
  | sort -u > included_hets.tsv

cut -f1 included_hets.tsv | sort | uniq -u > single_het_samples.txt

# rare intronic variants within ±20 of a boundary, any FILTER, in those samples
bcftools query \
    -S single_het_samples.txt \
    -i 'GT="alt"
        & CSQ_HGVSc~"c[.][-*]?[0-9]+[+-]([1-9]|1[0-9]|20)([^0-9]|$)"
        & (CSQ_gnomADg_AF="."      | CSQ_gnomADg_AF<=0.01)
        & (CSQ_gnomADe_AF="."      | CSQ_gnomADe_AF<=0.01)
        & (CSQ_TWE_AF="."          | CSQ_TWE_AF<=0.01)
        & (CSQ_TWE_WES_v1_AF="."   | CSQ_TWE_WES_v1_AF<=0.01)' \
    -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_HGVSc\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\t%CSQ_TWE_WES_v1_AF\n' \
    "$vcf" > single_het_intron_37.tsv