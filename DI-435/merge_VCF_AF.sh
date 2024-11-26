#!/bin/bash
set -exo pipefail
# Script for merging VCF files and adding AF tag
#
# Inputs:
#   $1 -> input file with VCF files to merge
#   $2 -> Reference genome file

input_file=$1
genome=$2

job="${DX_JOB_ID}"
num_proc="$(grep -c ^processor /proc/cpuinfo)"

# Check if the input file exists and is not empty
if [ ! -s "$input_file" ]; then
    echo "ERROR: Input file '$input_file' is empty or does not exist!"
    exit 1
fi

# Validate if the job ID is set and length is not zero
if [ -z "$job" ]; then
    echo "ERROR: DNAnexus job ID is not set!"
    exit 1
fi

# Validate genome file
if [ ! -f "$genome" ]; then
    echo "ERROR: Reference genome file '$genome' does not exist!"
    exit 1
fi

#blank array to store the project files
project_files=()

# Read the input file line by line and construct the project_file array
# keep track of line number for error reporting
line_number=0
while IFS=$'\t' read -r _ _ field3 field4 _; do
        ((line_number++))

        # Validate fields are not empty
        if [ -z "$field3" ] || [ -z "$field4" ]; then
            echo "ERROR: Missing project or file ID at line $line_number"
            exit 1
        fi

        # Validate field format (assuming they should match dx:// format)
        if [[ ! "$field3" =~ ^project-.* ]] || [[ ! "$field4" =~ ^file-.* ]]; then
            echo "ERROR: Invalid project or file ID format at line $line_number"
            exit 1
        fi

        project_files+=("${field3}:${field4}")
done < "$input_file"

echo "Downloading files"
echo "${project_files[@]}" | tr ' ' '\n' | xargs -n 1 -P "${num_proc}" -I {} dx download --no-progress "{}"

# Index VCFs
echo "Indexing VCFs"
echo *vcf.gz | tr ' ' '\n' | xargs -n 1 -P "${num_proc}" -I {} bcftools index "{}"

# Normalising VCFs
mkdir norm
echo "Normalising VCFs"
for vcf in *vcf.gz; do
    bcftools norm -m -any -f "${genome}" -Oz "${vcf}" > "norm/${vcf}"
done

# Indexing normalised VCFs
echo "Indexing normalised VCFs"
cd norm || exit
echo *vcf.gz | tr ' ' '\n' | xargs -n 1 -P "${num_proc}" -I {} bcftools index -f "{}"

# Merging normalised VCFs
echo "Collating normalised VCFs"
vcf_files=()
for vcf in *vcf.gz; do
    vcf_files+=("$vcf")
done
echo "Merging collated VCFs into ../merged.vcf"
bcftools merge --output-type v -m none --missing-to-ref "${vcf_files[@]}" > ../merged.vcf

# Bgzip and index merged VCF file
echo "Bgzip and indexing merged file"
cd ..
bgzip merged.vcf
bcftools index merged.vcf.gz

# Final merging and processing
echo "Final processing"

echo "Normalising merged VCF"
bcftools norm -m -any -f "${genome}" -Ou merged.vcf.gz \
    | bcftools +fill-tags --output-type v -o merge_tag.vcf -- -t AN,AC,NS,AF,MAF,AC_Hom,AC_Het,AC_Hemi

echo "Sorting and compressing final VCF"
bcftools sort merge_tag.vcf -Oz -o "final_merged_${job}.vcf.gz"

echo "Indexing final VCF"
tabix -p vcf "final_merged_${job}.vcf.gz"

dx upload "final_merged_${job}.vcf.gz"
dx upload "final_merged_${job}.vcf.gz.tbi"
dx terminate "${job}"