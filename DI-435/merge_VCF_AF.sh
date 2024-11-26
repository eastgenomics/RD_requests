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

exec 1> >(tee -a -i "/home/dnanexus/dx_stdout")
exec 2> >(tee -a -i "/home/dnanexus/dx_stderr")

# Prefixes all lines of commands written to stdout with datetime
PS4='\000[$(date)]\011'
export TZ=Europe/London

# set frequency of instance usage in logs to 10 seconds
sudo kill $(ps aux | grep pcp-dstat | head -n1 | awk '{print $2}')
/usr/bin/dx-dstat 10

_get_peak_usage() {
    : '''
    Reports the peak memory and storage usage from dstat, to be called at end of script
    '''
    dx watch $DX_JOB_ID --no-follow --quiet > job.log

    peak_mem=$(grep 'INFO CPU' job.log | cut -d':' -f5 | cut -d'/' -f1 | sort -n | tail -n1)
    total_mem="$(($(grep MemTotal /proc/meminfo | grep --only-matching '[0-9]*')/1024))"

    peak_storage=$(grep 'INFO CPU' job.log | cut -d':' -f6 | cut -d'/' -f1 | sort -n | tail -n1)
    total_storage=$(df -Pk / | awk 'NR == 2' | awk '{printf("%.0f", $2/1024/1024)}')

    echo "Memory usage peaked at ${peak_mem}/${total_mem}MB"
    echo "Storage usage peaked at ${peak_storage}/${total_storage}GB"
}

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

# Blank array to store the project files
project_files=()

# Read the input file line by line and construct the project_file array
# keep track of line number for error reporting
line_number=0
while IFS=$'\t' read -r _ field2 field3; do
    ((line_number++))

    # Validate fields are not empty
    if [ -z "$field2" ] || [ -z "$field3" ]; then
        echo "ERROR: Missing project or file ID at line $line_number"
        exit 1
    fi

    # Validate field format (assuming they should match dx:// format)
    if [[ ! "$field2" =~ ^project-.* ]] || [[ ! "$field3" =~ ^file-.* ]]; then
        echo "ERROR: Invalid project or file ID format at line $line_number"
        exit 1
    fi
    project_files+=("${field2}:${field3}")

done < "$input_file"

# Download VCFs in parallel
echo "Downloading files"
echo "${project_files[@]}" | tr ' ' '\n' | xargs -n 1 -P "${num_proc}" -I {} dx download --no-progress "{}"

# Index VCFs
echo "Indexing VCFs"
find . -maxdepth 1 -name '*.vcf.gz' -print0 | xargs -n 1 -P "${num_proc}" -I {} bcftools index "{}"

# Normalise multiallelics and left align VCFs in parallel
echo "Normalising VCFs"
mkdir -p norm
find . -maxdepth 1 -name '*.vcf.gz' -print0 | xargs -n 1 -P "${num_proc}" -I{} bcftools norm -m -any -f "${genome}" -Oz "{}" -o norm/"$(basename {})"

# Indexing normalised VCFs
echo "Indexing normalised VCFs"
cd norm || exit
find . -maxdepth 1 -name '*.vcf.gz' -print0 | xargs -0 -P8 -I{} bcftools index -f "{}"

# Merging normalised VCFs
echo "Collating VCFs"
file_prefix="${input_file%.txt}"
vcf_files=()
for vcf in *vcf.gz; do
    vcf_files+=("$vcf")
done

echo "Merging collated VCFs into ../merged.vcf"
bcftools merge --output-type v -m none --missing-to-ref "${vcf_files[@]}" > "../${file_prefix}_merged.vcf"

# Bgzip and index merged VCF file
echo "Bgzipping and indexing merged file"
cd ..
bgzip "${file_prefix}_merged.vcf"
bcftools index "${file_prefix}_merged.vcf.gz"

# Final merging and processing
echo "Final processing"

echo "Normalising merged VCF"
bcftools norm -m -any -f "${genome}" -Ou "${file_prefix}_merged.vcf.gz" \
    | bcftools +fill-tags --output-type v -o "${file_prefix}_merge_tag.vcf" \
    -- -t AN,AC,NS,AF,MAF,AC_Hom,AC_Het,AC_Hemi

echo "Sorting and compressing final VCF"
bcftools sort "${file_prefix}_merge_tag.vcf" -Oz -o "final_merged_${file_prefix}.vcf.gz"

echo "Indexing final VCF"
tabix -p vcf "final_merged_${file_prefix}.vcf.gz"

dx upload "final_merged_${file_prefix}.vcf.gz"
dx upload "final_merged_${file_prefix}.vcf.gz.tbi"
_get_peak_usage
dx terminate "${job}"
