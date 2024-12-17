#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 11:50:04 2024

@author: arun

Code logic:
    1. get a list of bam files from the project of Staging area 055?
    2. get the list of jsons from the same project.
    3. For each JSON, determine the sample name using regex.
    4. For each BAM, determine the sample name using regex
    5. Merge JSON and BAM files into a dataframe where the common samplename

"""
import re
import pandas as pd
import dxpy

# generate the list of bam files and json files rom the Medicover project
list_of_bam_ids = dxpy.bindings.search.find_data_objects(classname="file",
                                                         name="*.bam",
                                                         name_mode="glob",
                                                         project="project-GXx8Gg8468bZ072Py1Bb1p0F",
                                                         describe=True)

list_of_json_ids = dxpy.bindings.search.find_data_objects(classname="file",
                                                          name="*.json",
                                                          name_mode="glob",
                                                          project="project-GXx8Gg8468bZ072Py1Bb1p0F",
                                                          describe=True)

list_of_bams = [file['describe']['name'] for file in list_of_bam_ids]
list_of_jsons = [file['describe']['name'] for file in list_of_json_ids]


# Generate the list of sample names from the bam files
samplenames_bam = []
for bam in list_of_bams:
    pattern = r'SP[a-zA-Z0-9]+'
    match = re.search(pattern, bam)
    # If match not found, then try another pattern
    if match is None:
        pattern = r'[gG][mM][0-9]+_?[0-9]+'
        match = re.search(pattern, bam)
    
    samplenames_bam.append(match.group())
    # If match still not found, raise error
    if match is None:
        raise TypeError(f"The {bam} has no matched sample name")

# create a dataframe using bam file lists
data_bam = {"Sample_name": samplenames_bam, "Bamfile": list_of_bams}
# Drop duplicates
dataframe_bam = pd.DataFrame(data=data_bam).drop_duplicates(subset='Sample_name',
                                                            keep='first')


# Generate the list of sample names from the Json files
samplenames_json = []
for json in list_of_jsons:
    pattern = r'SP[a-zA-Z0-9]+'
    match = re.search(pattern, json)
    # If match not found, try another pattern
    if match is None:
        pattern = r'[gG][mM][0-9]+_?[0-9]+'
        match = re.search(pattern, json)
    # If match not found, append a None instead
    if match is None:
        print(f"The {json} has no matched sample name")
        """
        The following lines will be printed:
        The 20240809_LH00625_0016_B22NY3JLT3.json has no matched sample name
        The reports_20240726_LH00625_0014_A22HNW3LT3.json has no matched sample name
        The report_20240719_LH00625_0013_A22K2KNLT3.json has no matched sample name
        These files are not jsons belonging to a sample
        """
        samplenames_json.append(None)
        continue

    samplenames_json.append(match.group())

# Create a dataframe using json file lists
data_json = {"Sample_name": samplenames_json, "Json": list_of_jsons}
dataframe_json = pd.DataFrame(data=data_json).dropna()

# Merge the dataframe bam and dataframe json using the Sample_name as key
merged_df = pd.merge(dataframe_bam, dataframe_json,
                     on="Sample_name", how="left")

# Store the dataframe in separate tsv file
merged_df_with_json = merged_df.dropna()
merged_df_with_json.to_csv('Medicover_samples_with_json.tsv',
                           sep='\t', index=False)

merged_df_without_json = merged_df[merged_df['Json'].isna()]
merged_df_without_json.to_csv('Medicover_samples_without_json.tsv',
                              sep='\t', index=False)
