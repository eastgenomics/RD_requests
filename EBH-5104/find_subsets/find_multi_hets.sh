#!/bin/bash

# rather than looping through all hets, find sampls with 2 or more hets and then loop through those samples to get the variant info 

## GRCH38 ## 
# identify samples with 2 or more included hets
bcftools query \
-i 'GT="het" & FILTER!="EXCLUDE" & CSQ_Feature="NM_002303.6"' \
-f '[%SAMPLE\n]' \
lepr_variants_38_annotated_filtered.vcf.gz | sort | uniq -c | awk '$1 >=2 { print $2 }' > multi_hets.txt

# get het variants of samples with 2 or more included hets

bcftools query \
    -S multi_hets.txt \
    -i 'GT="het" & FILTER!="EXCLUDE" & CSQ_Feature="NM_002303.6"' \
    -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_WES_v1_AF\n]' \
    lepr_variants_38_annotated_filtered.vcf.gz > temp_included_multi_het_38.txt

sort temp_included_multi_het_38.txt > included_multi_het_38.txt

## GRCH37 ##
bcftools query \
 -i 'GT="het" & FILTER!="EXCLUDE" & CSQ_Feature="NM_002303.6"' \
 -f '[%SAMPLE\n]' \
 lepr_variants_37_annotated_filtered.vcf.gz | sort | uniq -c | awk '$1 >=2 { print $2 }' > multi_hets.txt


bcftools query \
    -S multi_hets.txt \
    -i 'GT="het" & FILTER!="EXCLUDE" & CSQ_Feature="NM_002303.6"' \
    -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\n]' \
    lepr_variants_37_annotated_filtered.vcf.gz > temp_included_multi_het_37.txt

sort temp_included_multi_het_37.txt > included_multi_het_37.txt