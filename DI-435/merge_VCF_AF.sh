#!/bin/bash
# Script for merging VCF files and adding AF tag
#
# Inputs:
#   $1 -> input file with VCF files to merge
#   $2 -> Reference genome file

input_file=$1
genome=$2

set -exo pipefail

exec 1> >(tee -a -i "/home/dnanexus/dx_stdout")
exec 2> >(tee -a -i "/home/dnanexus/dx_stderr")

# prefixes all lines of commands written to stdout with datetime
PS4='\000[$(date)]\011'
export TZ=Europe/London

# set frequency of instance usage in logs to 10 seconds
sudo kill $(ps aux | grep pcp-dstat | head -n1 | awk '{print $2}')
/usr/bin/dx-dstat 10

_get_peak_usage() {
    : '''
    Reports the peak memory and storage usage from dstat, to be called at end of app
    '''
    dx watch $DX_JOB_ID --no-follow --quiet > job.log

    peak_mem=$(grep 'INFO CPU' job.log | cut -d':' -f5 | cut -d'/' -f1 | sort -n | tail -n1)
    total_mem="$(($(grep MemTotal /proc/meminfo | grep --only-matching '[0-9]*')/1024))"

    peak_storage=$(grep 'INFO CPU' job.log | cut -d':' -f6 | cut -d'/' -f1 | sort -n | tail -n1)
    total_storage=$(df -Pk / | awk 'NR == 2' | awk '{printf("%.0f", $2/1024/1024)}')

    echo "Memory usage peaked at ${peak_mem}/${total_mem}MB"
    echo "Storage usage peaked at ${peak_storage}/${total_storage}GB"
}


# Download all VCFs in parallel
echo "Downloading VCFs"
to_download=$(awk -F ' ' '{print $2":"$3}' "$input_file")
echo "$to_download" | tr ' ' '\n' | xargs -P8 -n1 -I{} dx download {} --no-progress

# Index VCFs
echo "Indexing VCFs"
find . -maxdepth 1 -name '*.vcf.gz' -print0 | xargs -0 -P8 -I{} bcftools index -f "{}"

# Normalise multiallelics and left align VCFs in parallel
echo "Normalising VCFs"
mkdir -p norm
find . -maxdepth 1 -name '*.vcf.gz' -print0 | xargs -0 -P8 -I{} bcftools norm -m -any -f "${genome}" -Oz "{}" -o norm/"$(basename {})"

# Indexing normalised VCFs
echo "Indexing normalised VCFs"
cd norm || exit
find . -maxdepth 1 -name '*.vcf.gz' -print0 | xargs -0 -P8 -I{} bcftools index -f "{}"

# Merging normalised VCFs
echo "Merging normalised VCFs"
file_prefix="${input_file%.txt}"
command="bcftools merge --output-type v -m none --missing-to-ref"
# Add the VCF files names to the command
for vcf in ./*vcf.gz; do
    command="${command} $vcf";
done
command="${command} > ../${file_prefix}_merged.vcf"
echo "${command}"
eval "${command}"

# Bgzip and index merged VCF file
echo "Bgzipping and indexing merged file"
cd ..
bgzip "${file_prefix}_merged.vcf"
bcftools index "${file_prefix}_merged.vcf.gz"

# Final merging and processing
echo "Final processing"
command="bcftools norm -m -any -f ${genome} -Ou ${file_prefix}_merged.vcf.gz"
command="${command} | bcftools +fill-tags --output-type v -o ${file_prefix}_merge_tag.vcf -- -t AN,AC,NS,AF,MAF,AC_Hom,AC_Het,AC_Hemi"
command="${command} ; bcftools sort ${file_prefix}_merge_tag.vcf -Oz > final_merged_${file_prefix}.vcf.gz"
command="${command} ; tabix -p vcf final_merged_${file_prefix}.vcf.gz"
eval "$command"
dx upload "final_merged_${file_prefix}.vcf.gz"
dx upload "final_merged_${file_prefix}.vcf.gz.tbi"

_get_peak_usage
dx terminate "${DX_JOB_ID}"
