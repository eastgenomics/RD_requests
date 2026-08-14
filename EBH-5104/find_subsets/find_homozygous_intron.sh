#!/bin/bash

# get homozygous variants b38
for sample in $(bcftools query -l lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz);do
    old_annotation_vatiants=$(bcftools query \
    -s "$sample" \
    -i 'GT="AA" & CSQ_Consequence="intron_variant"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\n' \
    lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz) 
    if [ "$old_annotation_vatiants" = "" ]; then
        new_annotation_vatiants=$(bcftools query \
        -s "$sample" \
        -i 'GT="AA" & CSQ_Consequence="intron_variant"' \
        -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_WES_v1_AF\n' \
        lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz | uniq)
        if [ "$new_annotation_vatiants" != "" ]; then
        echo "$sample" uses new annotation >> homozygous_intron_38.txt
        echo "$new_annotation_vatiants" >> homozygous_intron_38.txt
        fi

    else
        echo "$sample" uses old annotation >> homozygous_intron_38.txt
        echo "$old_annotation_vatiants" >> homozygous_intron_38.txt
    fi
done


# get homozygous variants b37
for sample in $(bcftools query -l lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz);do
    old_annotation_vatiants=$(bcftools query \
    -s "$sample" \
    -i 'GT="AA" & CSQ_Consequence="intron_variant"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\n' \
    lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz | uniq)
    if [ "$old_annotation_vatiants" = "" ]; then
        new_annotation_vatiants=$(bcftools query \
        -s "$sample" \
        -i 'GT="AA" & CSQ_Consequence="intron_variant"' \
        -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_WES_v1_AF\n' \
        lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz | uniq)
        if [ "$new_annotation_vatiants" != "" ]; then
            echo "$sample" uses new annotation >> homozygous_intron_37.txt
            echo "$new_annotation_vatiants" >> homozygous_intron_37.txt
        fi
    else
        echo "$sample" uses old annotation >> homozygous_intron_37.txt
        echo "$old_annotation_vatiants" >> homozygous_intron_37.txt
    fi
done