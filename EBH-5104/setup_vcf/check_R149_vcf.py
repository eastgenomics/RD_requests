import argparse
import fnmatch
import os
import subprocess
import sys
import pandas as pd
import io as io
# want to read xlsx files
# in included find LEPR variants, in excluded find lepr variants 
R149_samples = pd.read_csv("R149_vcfs_VCFs_no_dups_no_control.tsv", sep="\t")
print(R149_samples)
lepr_samples_TWE_38 = pd.DataFrame()
lepr_samples_TWE_37 = pd.DataFrame()
vcf_38 = pd.DataFrame()
vcf_37 = pd.DataFrame()
build38_samples = R149_samples[R149_samples["project_name"].str.contains("_38_TWE")]
build37_samples = R149_samples[~R149_samples["project_name"].str.contains("_38_TWE")]
for file_id in build38_samples["project_file"]:
        # print(file_id)
    try :
            print(file_id)
            cmd = (
            f"dx cat {file_id} | bcftools view -T GRCh38_LEPR_exons_plus25.bed|  grep -v '^##'"
            )
            output = subprocess.run(cmd, shell=True,
                             capture_output=True, check=False)
            df = pd.read_csv(io.BytesIO(output.stdout), sep="\t")
            print(df)

            file_info = R149_samples[R149_samples["project_file"]==file_id].reset_index()
            sample_name = file_info._get_value(0,'sample')
            project_name = file_info._get_value(0,'project_name')            
            vcf_38= pd.concat([vcf_38, df])
            df[["project"]] = project_name
            df[["sample"]] = sample_name
            df.rename(columns={df.columns[9]: "sample_info"}, inplace=True)

            print(df)
            lepr_samples_TWE_38 = pd.concat([lepr_samples_TWE_38, df], axis=0, ignore_index= True)


    except ValueError as error:
            print(f"Error reading file {file_id}: {error}")
            print(f"Check archival status of {file_id}")
            pass
    
for file_id in build37_samples["project_file"]:
        # print(file_id)
    try :
            print(file_id)
            cmd = (
            f"dx cat {file_id} | bcftools view -T GRCh37_LEPR_exons_plus25.bed|  grep -v '^##'"
            )
            output = subprocess.run(cmd, shell=True,
                             capture_output=True, check=False)
            df = pd.read_csv(io.BytesIO(output.stdout), sep="\t")
            print(df)

            file_info = R149_samples[R149_samples["project_file"]==file_id].reset_index()
            vcf_37 = pd.concat([vcf_37, df])

            sample_name = file_info._get_value(0,'sample')
   
            project_name = file_info._get_value(0,'project_name')
            df[["project"]] = project_name
            df[["sample"]] = sample_name

            df.rename(columns={df.columns[9]: "sample_info"}, inplace=True)
            lepr_samples_TWE_37 = pd.concat([lepr_samples_TWE_37, df], axis=0, ignore_index= True)


    except ValueError as error:
            print(f"Error reading file {file_id}: {error}")
            print(f"Check archival status of {file_id}")
            pass
    
print(lepr_samples_TWE_38)
lepr_samples_TWE_38.to_csv(
        f"vcf_lepr_variants_TWE_38.tsv",
        sep="\t",
        index=False,
)
print(lepr_samples_TWE_37)
lepr_samples_TWE_37.to_csv(
        f"vcf_lepr_variants_TWE_37.tsv",
        sep="\t",
        index=False,
)

vcf_38.fillna("0/0:.:.:.:.", inplace=True)
print(vcf_38)
vcf_38.to_csv(
        f"lepr_variants_TWE_38.tsv",
        sep="\t",
        index=False,
)

vcf_37.fillna("0/0:.:.:.:.", inplace=True)
print(vcf_37)
vcf_37.to_csv(
        f"lepr_variants_TWE_37.tsv",
        sep="\t",
        index=False,
)