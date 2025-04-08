# Sentieon germline fastq to vcf update v4.2.2 to v5.1.0

```bash
TWE_samples=(135168069-25028R0075 135222258-25030R0051 135359325-25036R0111 135377736-25037R0041 135806665-25057Q0138)
CEN_samples=(135648585-25049R0117 135648637-25050R0031 135648671-25050R0028 135648731-25049R0118 135744929-25052R0109)

TWE_project=002_250305_A01303_0525_BH3L35DMX2_TWE
CEN_project=002_250306_A01295_0497_AHVJJTDRX5_CEN

for sample in ${TWE_samples[@]}; do file=$(dx find data --name ${sample}*.bam.bai --path ${TWE_project}:/output/TWE-250307_0030/sentieon-dnaseq-4.2.2/ --brief); dx cp $file 004_250327_Sentieon_RD_v5.1.0:/TWE_expected/; done
for sample in ${CEN_samples[@]}; do file=$(dx find data --name ${sample}*.bam.bai --path ${CEN_project}:/output/CEN-250307_1320/sentieon-dnaseq-4.2.2/ --brief); dx cp $file 004_250327_Sentieon_RD_v5.1.0:/CEN_expected/; done

# fastqs
TWE_fastq_ids=$(for sample in ${TWE_samples[@]}; do file=$(dx find data --name ${sample}*fastq.gz --path 001_Staging_Area52:/250305_A01303_0525_BH3L35DMX2/ --brief); echo $file; done)
CEN_fastq_ids=$(for sample in ${CEN_samples[@]}; do file=$(dx find data --name ${sample}*fastq.gz --path 001_Staging_Area52:/250306_A01295_0497_AHVJJTDRX5/ --brief); echo $file; done)
```

## Reproducibility

```bash
for bam in $(dx ls /TWE_expected/*bam); do id=$(echo ${bam} | cut -d"-" -f2); fastqs=$(dx find data --name *$id* --path 001_Staging_Area52:/250305_A01303_0525_BH3L35DMX2/Data/Intensities/BaseCalls/ --brief); dx cp $fastqs 004_250327_Sentieon_RD_v5.1.0:/TWE_inputs/; done

for bam in $(dx ls /CEN_expected/*bam); do id=$(echo ${bam} | cut -d"-" -f2); fastqs=$(dx find data --name *$id* --path 001_Staging_Area52:/250306_A01295_0497_AHVJJTDRX5/Data/Intensities/BaseCalls/ --brief); dx cp $fastqs 004_250327_Sentieon_RD_v5.1.0:/CEN_inputs/; done

python3 DI-1526/reproducibility.py file-Gz557F84qFpJ5GQX95B4FzFq file-Gz557GQ4qFp4Y3K7z93JX17V file-Gz55JVQ4qFp274ZBF7Ky6k4K file-Gz55JY84qFp6byYg842989x4
All bwa-mem jobs have started running
Here are the list of jobs for bwa-mem:
['job-GzXZ7y04JXpgfjf5075QVkbg', 'job-GzXZ7y84JXpVB7pG6qBQpgBp', 'job-GzXZ7y84JXpv7783kq499YkB', 'job-GzXZ7y84JXpVB7pG6qBQpgBv', 'job-GzXZ7yQ4JXpXZvpY4j6x5jZX', 'job-GzXZ7yQ4JXpk0YJYVjBPgbpF', 'job-GzXZ7yQ4JXpk0YJYVjBPgbpJ', 'job-GzXZ7yj4JXpXZvpY4j6x5jZZ', 'job-GzXZ7yj4JXpXZbkKjJYX7y6P', 'job-GzXZ7yj4JXpV2x2g12x6PJpG']

dx run app-cloud_workstation --ssh -y
mkdir reproducibility
sudo apt-get install tabix
i=1; for bam in $(dx find data --name *.bam --path 004_250327_Sentieon_RD_v5.1.0:/reproducibility/ --brief); do name_bam=$(dx describe ${bam} --json | jq -r .name); dx download ${bam} -o reproducibility/${i}_${name_bam}; i=$(($i+1)); done

i=1; for vcf in $(dx find data --name *Haplotyper.vcf.gz --path 004_250327_Sentieon_RD_v5.1.0:/reproducibility/ --brief); do name_vcf=$(dx describe ${vcf} --json | jq -r .name); dx download ${vcf} -o reproducibility/${i}_${name_vcf}; i=$(($i+1)); done

diff --from-file 1_135168069-25028R0075-25NGWES8-9526-F-103698_S13_L001_markdup.bam *bam > output_diff_bam
dx upload output_diff_bam --path 004_250327_Sentieon_RD_v5.1.0:/reproducibility/
```

## Samtools flagstat

```bash
dx mkdir -p /flagstat_comparison/TWE_v4.2.2/
dx mkdir -p /flagstat_comparison/CEN_v4.2.2/
# get the flagstat for TWE
for bam in $(dx ls /TWE_expected/*bam); do id=$(echo ${bam} | cut -d"-" -f2); flagstat=$(dx find data --name *${id}* --path 002_250305_A01303_0525_BH3L35DMX2_TWE:/output/TWE-250307_0030/eggd_samtools_flagstat-1.1.0/ --brief); dx cp $flagstat 004_250327_Sentieon_RD_v5.1.0:/flagstat_comparison/TWE_v4.2.2/; done

# get the flagstat for CEN
for bam in $(dx ls /CEN_expected/*bam); do id=$(echo ${bam} | cut -d"-" -f2); flagstat=$(dx find data --name *${id}* --path 002_250306_A01295_0497_AHVJJTDRX5_CEN:/output/CEN-250307_1320/eggd_samtools_flagstat-1.1.0/ --brief); dx cp $flagstat 004_250327_Sentieon_RD_v5.1.0:/flagstat_comparison/CEN_v4.2.2/; done

TWE_fastq_ids=$(for sample in ${TWE_samples[@]}; do file=$(dx find data --name ${sample}*fastq.gz --path 001_Staging_Area52:/250305_A01303_0525_BH3L35DMX2/ --brief); echo $file; done)
CEN_fastq_ids=$(for sample in ${CEN_samples[@]}; do file=$(dx find data --name ${sample}*fastq.gz --path 001_Staging_Area52:/250306_A01295_0497_AHVJJTDRX5/ --brief); echo $file; done)
python3 DI-1526/run_sentieon.py $TWE_fastq_ids -d /TWE_sentieon_v5.1.0_output/
python3 DI-1526/run_sentieon.py $CEN_fastq_ids -d /CEN_sentieon_v5.1.0_output/

# generate batch inputs for TWE
dx generate_batch_inputs -iinput_bam="(.*)_(.*)_(.*)_(.*).bam$" -iinput_bam_index="(.*)_(.*)_(.*)_(.*).bam.bai$" --path /TWE_sentieon_v5.1.0_output/ -o TWE_flagstat
dx generate_batch_inputs -iinput_bam="(.*)_(.*)_(.*)_(.*).bam$" -iinput_bam_index="(.*)_(.*)_(.*)_(.*).bam.bai$" --path /CEN_sentieon_v5.1.0_output/ -o CEN_flagstat

dx run app-GPJ8Zqj4F534YK4FFZ94QVqj --batch-tsv TWE_flagstat.0000.tsv --destination /flagstat_comparison/TWE_v5.1.0 -y
dx run app-GPJ8Zqj4F534YK4FFZ94QVqj --batch-tsv CEN_flagstat.0000.tsv --destination /flagstat_comparison/CEN_v5.1.0 -y

python3 DI-1526/compare_flagstat.py 135168069-25028R0075 135222258-25030R0051 135359325-25036R0111 135377736-25037R0041 135806665-25057Q0138 -p /flagstat_comparison/
135168069-25028R0075 OK
135222258-25030R0051 OK
135359325-25036R0111 OK
135377736-25037R0041 OK
135806665-25057Q0138 OK

python3 DI-1526/compare_flagstat.py 135648585-25049R0117 135648637-25050R0031 135648671-25050R0028 135648731-25049R0118 135744929-25052R0109 -p /flagstat_comparison/
135648585-25049R0117 OK
135648637-25050R0031 OK
135648671-25050R0028 OK
135648731-25049R0118 OK
135744929-25052R0109 OK
```

## Hap.py

```bash
python3 DI-1526/run_happy.py /TWE_expected /TWE_sentieon_v5.1.0_output TWE -d /happy_TWE/
Here are the list of jobs for happy:
['job-GzfP4504JXpzXYGp626p0J70', 'job-GzfP4504JXpkkgPBgV8QPyy0', 'job-GzfP4584JXpx347VGQfGJ71Z', 'job-GzfP4584JXpkkgPBgV8QPyy2', 'job-GzfP4584JXpq4vg6fjF03xqf']

python3 DI-1526/run_happy.py /CEN_expected /CEN_sentieon_v5.1.0_output CEN -d /happy_CEN/
Here are the list of jobs for happy:
['job-GzfP48Q4JXpjY5QV0jz0yQKp', 'job-GzfP48Q4JXpV0QybjpV1FP9x', 'job-GzfP48Q4JXpQ7Jpzf7X75KFg', 'job-GzfP48j4JXpQ7Jpzf7X75KFk', 'job-GzfP48j4JXpz7pGz58xQxXYY']

for file in $(dx find data --name *.csv --path /happy_CEN --brief); do dx cat $file | cut -d"," -f10,11; done
```
