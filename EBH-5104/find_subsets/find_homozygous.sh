#!/bin/bash

# b37 files uses the annotation CSQ_TWE_AF
# b38 files uses the annotation CSQ_TWE_WES_v1_AF
# get b38 homozygous variants
bcftools query \
    -i 'GT="AA" & FILTER!="EXCLUDE" & CSQ_Feature="NM_002303.6"' \
    -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_WES_v1_AF\n]' \
    lepr_variants_38_annotated_filtered.vcf.gz | sort -u >> included_homozygous_38.txt 

# get b37 homozygous variants
bcftools query \
    -i 'GT="AA" & FILTER!="EXCLUDE" & CSQ_Feature="NM_002303.6"' \
    -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\n]' \
    lepr_variants_37_annotated_filtered.vcf.gz | sort -u >> included_homozygous_37.txt
