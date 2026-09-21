#!/bin/bash

# get samples with 1 heterozygous variant
hets=$(bcftools query \
-s "$sample" \
-i 'GT="het" & FILTER!="EXCLUDE"' \
-f '[%SAMPLE]' \
lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz | sort | uniq -c | awk '$1 ==1 { print $2 }' )

# search for intronic variants in those samples
# some samples have old annotation for TWE frequency (CSQ_TWE_AF)
# some have new annotation (CSQ_TWE_WES_v1_AF), so need to check both
# want AF for indication of rare intron variants
for sample in $hets; do
    intron_variants_old_annotation=$(bcftools query \
        -s "$sample" \
        -i 'GT!="./." &CSQ_Consequence~"intron*"' \
        -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\n' \
        lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz)
        if [ "$intron_variants_old_annotation" = "" ]; then
            intron_variants_new_annotation=$(bcftools query \
            -s "$sample" \
            -i 'GT!="./." &CSQ_Consequence~"intron*"' \
            -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_WES_v1_AF\n' \
            lepr_variants_TWE_38_sorted_annotated_filtered.vcf.gz)

            count_intron=$(echo "$intron_variants_new_annotation" | grep -c .)
            if [ $count_intron -ge 1 ]; then
                echo "$hets" >> single_het_intron_38.txt
                echo "$intron_variants_new_annotation" >> single_het_intron_38.txt
            fi
        else

            count_intron=$(echo "$intron_variants_old_annotation" | grep -c .)
            if [ $count_intron -ge 1 ]; then
                echo "$hets" >> single_het_intron_38.txt
                echo "$intron_variants_old_annotation" >> single_het_intron_38.txt
            fi
    fi

done


# get heterozygous variants b37


    # check for 1 heterozygous variant
hets=$(bcftools query \
    -s "$sample" \
    -i 'GT="het" & FILTER!="EXCLUDE"' \
    -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\n' \
    lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz | sort | uniq -c | awk '$1 ==1 { print $2 }' )

for sample in $hets; do
    intron_variants_old_annotation=$(bcftools query \
        -s "$sample" \
        -i 'GT!="0/0" &CSQ_Consequence~"intron_variant"' \
        -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\n' \
        lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz)
        if [ "$intron_variants_old_annotation" = "" ]; then
            intron_variants_new_annotation=$(bcftools query \
            -s "$sample" \
            -i 'GT!="0/0" &CSQ_Consequence~"intron_variant"' \
            -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_WES_v1_AF\n' \
            lepr_variants_TWE_37_sorted_annotated_filtered.vcf.gz)

            count_intron=$(echo "$intron_variants_new_annotation" | grep -c .)
            if [ $count_intron -ge 1 ]; then
                echo "$hets" >> single_het_intron_37.txt
                echo "$intron_variants_new_annotation" >> single_het_intron_37.txt
            fi
        else

            count_intron=$(echo "$intron_variants_old_annotation" | grep -c .)
            if [ $count_intron -ge 1 ]; then
                echo "$hets" >> single_het_intron_37.txt
                echo "$intron_variants_old_annotation" >> single_het_intron_37.txt
            fi
        fi
done