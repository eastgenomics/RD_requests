#!/bin/bash

# get heterozygous variants b38
for sample in $(bcftools query -l lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz);do
    count=$(bcftools query \
    -s "$sample" \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz | uniq |  wc -l)
    if [ $count -ge 2 ]; then
        bcftools query \
    -s "$sample" \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz | uniq >> included_multi_het_38.txt
    fi
done

# get heterozygous variants b37
for sample in $(bcftools query -l lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz);do
    count=$(bcftools query \
    -s "$sample" \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz | uniq |  wc -l)
    if [ $count -ge 2 ]; then
        bcftools query \
    -s "$sample" \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz | uniq >> included_multi_het_37.txt
    fi
done