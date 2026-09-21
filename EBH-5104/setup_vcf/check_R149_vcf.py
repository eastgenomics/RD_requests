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
vcf_38_list = []
vcf_37_list = []
build38_samples = R149_samples[R149_samples["project_name"].str.contains("_38_TWE")]
build37_samples = R149_samples[~R149_samples["project_name"].str.contains("_38_TWE")]

# check if bed files exist in current folder, if not, print message and exit
try:
        bed_file = "GRCh38_LEPR_exons_plus25.bed"
        with open(bed_file, 'r') as file:
                print("GRCh38_LEPR_exons_plus25.bed exists and is ready to read")
except FileNotFoundError:
        print("GRCh38_LEPR_exons_plus25.bed does not exist in the current folder.")

try:
        bed_file = "GRCh37_LEPR_exons_plus25.bed"
        with open(bed_file, 'r') as file:
                print("GRCh37_LEPR_exons_plus25.bed exists and is ready to read")
except FileNotFoundError:
        print("GRCh37_LEPR_exons_plus25.bed does not exist in the current folder.")

# run if there are build 38 samples, otherwise skip to build 37 samples
if len(build38_samples) > 0:
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
                except pd.errors.EmptyDataError as error:
                        print(f"Error reading file {file_id}: {error}")
                        continue
                file_info = R149_samples[R149_samples["project_file"]==file_id].reset_index()
                sample_name = file_info.iloc[0]['sample']
                project_name = file_info.iloc[0]['project_name']
                df[["project"]] = project_name
                df[["sample"]] = sample_name
                df.rename(columns={df.columns[9]: "sample_info"}, inplace=True)
                vcf_38_list.append(df)

        lepr_samples_TWE_38 = pd.concat(vcf_38_list, axis=0, ignore_index= True)
        print(lepr_samples_TWE_38)
        # one line per variant/sample combo
        lepr_samples_TWE_38.to_csv(
                f"lepr_variants_TWE_38.tsv",
                sep="\t",
                index=False,
        )

        # convert to multi-sample vcf format, one line per variant, with sample columns
        variant_cols = ["#CHROM","POS","ID", "REF", "ALT"]
        vcf_38_lepr = (
        lepr_samples_TWE_38.pivot_table(
                index=variant_cols,
                columns="sample",
                values="sample_info",
                aggfunc="first"
        )
        .reset_index()
        )
        metadata = (
        lepr_samples_TWE_38.groupby(variant_cols)
        .agg({
          "QUAL": "first",   
          "FILTER": "first",
          "INFO": "first",
          "FORMAT": "first"
         })
         .reset_index()
        )

        vcf_38_lepr = metadata.merge(
        vcf_38_lepr,
        on=variant_cols,
        how="left"
        )

        # fill in missing sample values with "./."
        sample_cols = [c for c in vcf_38_lepr.columns if c not in variant_cols]
        vcf_38_lepr[sample_cols] = vcf_38_lepr[sample_cols].fillna("./.")

        # save to vcf file
        with open("lepr_variants_38.vcf", "w") as f:
                f.write("##fileformat=VCFv4.2\n")
                vcf_38_lepr.to_csv(f, sep="\t", index=False)


# repeat for build 37 samples
if len(build37_samples) > 0:
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

        except pd.errors.EmptyDataError as error:
                print(f"Error reading file {file_id}: {error}")
                print(f"Check archival status of {file_id}")
                continue
    
        file_info = R149_samples[R149_samples["project_file"]==file_id].reset_index()
        

        sample_name = file_info.iloc[0]['sample']
        project_name = file_info.iloc[0]['project_name']
        df[["project"]] = project_name
        df[["sample"]] = sample_name

        df.rename(columns={df.columns[9]: "sample_info"}, inplace=True)
        vcf_37_list.append(df)

        lepr_samples_TWE_37 = pd.concat(vcf_37_list, axis=0, ignore_index= True)

        print(lepr_samples_TWE_37)
        lepr_samples_TWE_37.to_csv(
                f"lepr_variants_TWE_37.tsv",
                sep="\t",
                index=False,
        )

        # convert to multi-sample vcf format, one line per variant, with sample columns
        variant_cols = ["#CHROM","POS","ID", "REF", "ALT"]
        vcf_37_lepr = (
        lepr_samples_TWE_37.pivot_table(
                index=variant_cols,
                columns="sample",
                values="sample_info",
                aggfunc="first"
        )
        .reset_index()
        )
        # gather metadata for each variant, since pivot_table will drop these columns
        metadata = (
        lepr_samples_TWE_37.groupby(variant_cols)
        .agg({
          "QUAL": "first",   
          "FILTER": "first",
          "INFO": "first",
          "FORMAT": "first"
         })
         .reset_index()
        )
        # add variant metadata back to vcf_37_lepr
        vcf_37_lepr = metadata.merge(
        vcf_37_lepr,
        on=variant_cols,
        how="left"
        )

        # fill in missing sample values with "./."
        # add co;lumns for each sample, with missing values filled in with "./."
        sample_cols = [c for c in vcf_37_lepr.columns if c not in variant_cols]
        vcf_37_lepr[sample_cols] = vcf_37_lepr[sample_cols].fillna("./.")

        # save to vcf file
        with open("lepr_variants_37.vcf", "w") as f:
                f.write("##fileformat=VCFv4.2\n")
                vcf_37_lepr.to_csv(f, sep="\t", index=False)
