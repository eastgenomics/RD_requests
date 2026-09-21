#!/bin/bash

# get b38 homozygous variants
bcftools query \
    -i 'GT="AA" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz | uniq >> included_homozygous_38.txt 

# get b37 homozygous variants
bcftools query \
    -i 'GT="AA" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz | uniq >> included_homozygous_37.txt
