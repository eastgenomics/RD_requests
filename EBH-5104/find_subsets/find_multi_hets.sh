#!/bin/bash

# rather than looping through all hets, find sampls with 2 or more hets and then loop through those samples to get the variant info 

## GRCH38 ## 
multi_hets=$(bcftools query \
-i 'GT="het" & FILTER!="EXCLUDE"' \
-f '[%SAMPLE]' \
lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz | sort | uniq -c | awk '$1 >=2 { print $2 }' )

for sample in $multi_hets; do
    bcftools query \
    -s "$sample" \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz >> included_multi_het_38.txt
done


## GRCH37 ##
multi_hets=$(bcftools query \
-i 'GT="het" & FILTER!="EXCLUDE"' \
-f '[%SAMPLE]' \
lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz | sort | uniq -c | awk '$1 >=2 { print $2 }' )

for sample in $multi_hets; do
    bcftools query \
    -s "$sample" \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz >> included_multi_het_37.txt
done