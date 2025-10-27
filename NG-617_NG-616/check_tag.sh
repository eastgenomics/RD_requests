#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 <tag> <clinvar_dev_file> <clinvar_prod_file>"
    echo "Example: bash $0 'Pathogenic,_low_penetrance' clinvar_20251019_GRCh38.vcf.gz clinvar_20250907_GRCh38.vcf.gz"
    exit 1
}

# Check if no arguments provided
if [ $# -eq 0 ]; then
    usage
fi

# Declare inputs
tag="$1"
clinvar_dev="$2"
clinvar_prod="$3"

printf "Processing tag: %s\n" "$tag"
# This bcftools command will generate a tab separated table of 3 rows, CHROM:POS, CLNSIG and CLNHGVS
bcftools query -f "%CHROM:%POS\t%CLNSIG\t%CLNHGVS" $clinvar_dev | grep $'\t'${tag}$'\t' | \
while IFS=$'\t' read -r chrom_pos clnsig clnhgvs; do
    echo "Processing variant: $chrom_pos"
    echo "CLNSIG: $clnsig"
    echo "CLNHGVS: $clnhgvs"
    echo "---"

    results=$(bcftools query -r "${chrom_pos}" \
        --regions-overlap 0 \
        -f "%CHROM\t%POS\t%REF\t%ALT\t%CLNSIG\t%CLNHGVS" \
        $clinvar_prod | grep "$clnhgvs")

    # Check if the result is empty
    if [ -n "$results" ]; then
        echo "Found variant:"
        echo "$results"
    else
        echo "Variant in $chrom_pos could not be found"
    fi
    echo "================================"
done

