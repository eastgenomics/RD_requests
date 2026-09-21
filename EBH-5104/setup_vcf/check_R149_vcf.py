import argparse
import fnmatch
import subprocess
import pandas as pd
import io as io
# want to readvcf files
def parse_arguments():
        """
        Parse command-line arguments.

        Returns
        -------
        args: argparse.Namespace
        Parsed arguments
        """
        parser = argparse.ArgumentParser(
                description="read VCFs and build table of LEPR variants"
        )
        parser.add_argument(
        "--input_file", type=str, required=True,
        help="Input file containing VCF file information."
        )
        return parser.parse_args()

# in included find LEPR variants, in excluded find lepr variants 
args = parse_arguments()
R149_samples = pd.read_csv(args.input_file, sep="\t")

print(R149_samples)
lepr_samples_TWE_38 = pd.DataFrame()
lepr_samples_TWE_37 = pd.DataFrame()
vcf_38 = pd.DataFrame()
vcf_37 = pd.DataFrame()
build38_samples = R149_samples[R149_samples["project_name"].str.contains("_38_TWE")]
build37_samples = R149_samples[~R149_samples["project_name"].str.contains("_38_TWE")]

# check bed files exist 
try:
    # File path
        bed_file = "GRCh38_LEPR_exons_plus25.bed"
    
        # Try to open the file
        with open(bed_file, 'r') as file:
                print("GRCh38_LEPR_exons_plus25.bed exists and is ready to read")
except FileNotFoundError:
        print("GRCh38_LEPR_exons_plus25.bed does not exist in the current folder.")

try:
    # File path
        bed_file = "GRCh37_LEPR_exons_plus25.bed"
    
        # Try to open the file
        with open(bed_file, 'r') as file:
                print("GRCh37_LEPR_exons_plus25.bed exists and is ready to read")
except FileNotFoundError:
        print("GRCh37_LEPR_exons_plus25.bed does not exist in the current folder.")

for file_id in build38_samples["project_file"]:
        # check file is live and not archived, if archived, skip and print message
        is_live_cmd = (
            f"dx describe {file_id} --json | jq -r '.archivalState'"
            )
        is_live_output = subprocess.run(is_live_cmd, shell=True,
                             capture_output=True, check=True)
        archival_state = is_live_output.stdout.decode().strip()
        if archival_state != "live":
            print(f"File {file_id} is archived. Skipping.")
            continue
        else:
              print(f"File {file_id} is live. Proceeding with processing.")
        try :
            print(file_id)
            cmd = (
            f"dx cat {file_id} | bcftools view -T GRCh38_LEPR_exons_plus25.bed|  grep -v '^##'"
            )
            try:        
                output = subprocess.run(cmd, shell=True,
                             capture_output=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command for file {file_id}: {e}")
                print(f"Command output: {subprocess.run(cmd, shell=True,
                             capture_output=True).stderr.decode()}")
                continue
            df = pd.read_csv(io.BytesIO(output.stdout), sep="\t")
            print(df)
        except pd.errors.EmptyDataError as error:
            print(f"Error reading file {file_id}: {error}")
            continue
        file_info = R149_samples[R149_samples["project_file"]==file_id].reset_index()
        sample_name = file_info.iloc[0]['sample']
        project_name = file_info.iloc[0]['project_name']
        vcf_38= pd.concat([vcf_38, df])
        df[["project"]] = project_name
        df[["sample"]] = sample_name
        df.rename(columns={df.columns[9]: "sample_info"}, inplace=True)

        print(df)
        lepr_samples_TWE_38 = pd.concat([lepr_samples_TWE_38, df], axis=0, ignore_index= True)

vcf_38.fillna("0/0:.:.:.:.", inplace=True)
print(vcf_38)
vcf_38.to_csv(
        f"vcf_format_lepr_variants_TWE_38.tsv",
        sep="\t",
        index=False,
)


for file_id in build37_samples["project_file"]:
        # print(file_id)
                # check file is live and not archived, if archived, skip and print message
        is_live_cmd = (
            f"dx describe {file_id} --json | jq -r '.archivalState'"
            )
        is_live_output = subprocess.run(is_live_cmd, shell=True,
                             capture_output=True, check=True)
        archival_state = is_live_output.stdout.decode().strip()
        if archival_state != "live":
            print(f"File {file_id} is archived. Skipping.")
            continue
        else:
              print(f"File {file_id} is live. Proceeding with processing.")
        try :
            print(file_id)
            cmd = (
            f"dx cat {file_id} | bcftools view -T GRCh37_LEPR_exons_plus25.bed|  grep -v '^##'"
            )
            try:        
                output = subprocess.run(cmd, shell=True,
                             capture_output=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command for file {file_id}: {e}")
                print(f"Command output: {subprocess.run(cmd, shell=True,
                             capture_output=True).stderr.decode()}")
                continue
            df = pd.read_csv(io.BytesIO(output.stdout), sep="\t")
            print(df)

        except pd.errors.EmptyDataError as error:
                print(f"Error reading file {file_id}: {error}")
                print(f"Check archival status of {file_id}")
                continue
    
        file_info = R149_samples[R149_samples["project_file"]==file_id].reset_index()
        vcf_37 = pd.concat([vcf_37, df])

        sample_name = file_info.iloc[0]['sample']
        project_name = file_info.iloc[0]['project_name']
        df[["project"]] = project_name
        df[["sample"]] = sample_name

        df.rename(columns={df.columns[9]: "sample_info"}, inplace=True)
        lepr_samples_TWE_37 = pd.concat([lepr_samples_TWE_37, df], axis=0, ignore_index= True)

vcf_37.fillna("0/0:.:.:.:.", inplace=True)
print(vcf_37)
vcf_37.to_csv(
        f"vcf_format_lepr_variants_TWE_37.tsv",
        sep="\t",
        index=False,
)