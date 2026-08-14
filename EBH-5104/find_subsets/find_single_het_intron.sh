#!/bin/bash

# get heterozygous variants b38
for sample in $(bcftools query -l lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz);do

    # check for 1 heterozygous variant
    hets=$(bcftools query \
    -s "$sample" \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz | uniq)
    count=$(echo "$hets" | wc -l)
    if [ "$hets" = "" ]; then
    # $STRING is empty
        continue
    fi
    if [ $count -eq 1 ]; then
        
        intron_variants_old_annotation=$(bcftools query \
            -s "$sample" \
            -i 'GT!="0/0" &CSQ_Consequence="intron_variant"' \
            -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\n' \
            lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz| uniq)
            if [ "$intron_variants_old_annotation" = "" ]; then
                intron_variants_new_annotation=$(bcftools query \
                -s "$sample" \
                -i 'GT!="0/0" &CSQ_Consequence="intron_variant"' \
                -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_WES_v1_AF\n' \
                lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz | uniq )

                count_intron=$(echo "$intron_variants_new_annotation" | wc -l)
                if [ $count_intron -ge 1 ]; then
                    echo "$hets" >> single_het_intron_38.txt
                    echo "$intron_variants_new_annotation" >> single_het_intron_38.txt
                fi
            else

                count_intron=$(echo "$intron_variants_old_annotation" | wc -l)
                if [ $count_intron -ge 1 ]; then
                    echo "$hets" >> single_het_intron_38.txt
                    echo "$intron_variants_old_annotation" >> single_het_intron_38.txt
                fi
        fi


    # print intronic variants

    fi
done


# get heterozygous variants b37
for sample in $(bcftools query -l lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz);do

    # check for 1 heterozygous variant
    hets=$(bcftools query \
    -s "$sample" \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz | uniq)
    count=$(echo "$hets" | wc -l)
    if [ "$hets" = "" ]; then
    # $STRING is empty
        continue
    fi
    if [ $count -eq 1 ]; then
        
        intron_variants_old_annotation=$(bcftools query \
            -s "$sample" \
            -i 'GT!="0/0" &CSQ_Consequence="intron_variant"' \
            -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\n' \
            lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz| uniq)
            if [ "$intron_variants_old_annotation" = "" ]; then
                intron_variants_new_annotation=$(bcftools query \
                -s "$sample" \
                -i 'GT!="0/0" &CSQ_Consequence="intron_variant"' \
                -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_WES_v1_AF\n' \
                lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz | uniq )

                count_intron=$(echo "$intron_variants_new_annotation" | wc -l)
                if [ $count_intron -ge 1 ]; then
                    echo "$hets" >> single_het_intron_37.txt
                    echo "$intron_variants_new_annotation" >> single_het_intron_37.txt
                fi
            else

                count_intron=$(echo "$intron_variants_old_annotation" | wc -l)
                if [ $count_intron -ge 1 ]; then
                    echo "$hets" >> single_het_intron_37.txt
                    echo "$intron_variants_old_annotation" >> single_het_intron_37.txt
                fi
        fi


    fi
done