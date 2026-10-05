#!/bin/bash

# get homozygous intron variants b38
# b37 files uses the annotation CSQ_TWE_AF
# b38 files uses the annotation CSQ_TWE_WES_v1_AF
# want AF for indication of rare intron variants
# adding  2>/dev/null as the regex CSQ_Consequence~"intron_variant" returns an evaluation for each variant
# ie  pass=1 [intron_variant] or pass=0 [missense_variant]
bcftools query \
 -i 'GT="AA" & CSQ_Consequence~"intron_variant"' \
 -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%FILTER\t%CSQ_Feature\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_WES_v1_AF\n]' \
 lepr_variants_38_annotated_filtered.vcf.gz 2>/dev/null > homozygous_intron_38.txt


# get homozygous variants b37
bcftools query \
 -i 'GT="AA" & CSQ_Consequence~"intron_variant"' \
 -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%FILTER\t%CSQ_Feature\t%CSQ_Consequence\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_TWE_AF\n]' \
    lepr_variants_37_annotated_filtered.vcf.gz 2>/dev/null > homozygous_intron_37.txt

