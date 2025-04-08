import argparse
import dxpy
import pandas as pd
import swifter

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--test_code_mapping",
        type=str,
        help=("Path to file mapping R code to NGS assay (CEN or TWE)"),
    )
    parser.add_argument(
        "--clarity_extract",
        type=str,
        help=("Path to Clarity extract file"),
    )

    return parser.parse_args()

def open_files(clarity, assay_mapping):
    '''
    Open files and read in file contents to DataFrames
    '''
    with open(clarity) as f:
        clarity_df = pd.read_csv(f, delimiter='\t')

    with open(assay_mapping) as f:
        assay_df = pd.read_csv(f, names=['Test_Code', 'Assay'])

    return clarity_df, assay_df

def create_path(filename, assay, run):
    '''
    Create a path to a specific file on clingen
    Inputs:
        assay (str): either CEN or WES
        filename (str): filename for the xlsx report
        run (str): sequencing run name
    Outputs:
        path (str): path to the given filename on clingen
    '''
    base_path = "/appdata/clingen/cg/Regional Genetics Laboratories/Molecular Genetics/Data archive/Sequencing HT/"
    if assay == "CEN":
        path = base_path + f"{assay}/Run folders/{run}_{assay}/{filename}"
    elif assay == "WES":
        path = base_path + f"{assay}/{run}_TWE/{filename}"
    return path

def find_file_name(search_query):
    '''
    Find the file name for a given search query on DNAnexus
    Inputs:
        search_query (str): search query for file name when searching DNAnexus
    Outputs:
        filename (str): a file name, or None if no or multiple matches found
    '''
    files = list(dxpy.find_data_objects(
        name=search_query,
        name_mode="glob",
        describe=True
    ))
    for file in files:
        filenames = []
        filenames.append(file.get('describe').get('name'))
        if len(filenames) != 1:
            print(f"{search_query} returned multiple/no files")
            return None
        else:
            return filenames[0]


def main():
    '''
    Script entry point
    '''
    # Parse args
    args = parse_args()

    # Read data into dataframes
    clarity_df, assay_df = open_files(args.clarity_extract, args.test_code_mapping)

    # Process data to construct a path for each specimen
    clarity_df = clarity_df.astype(str)
    df = pd.merge(clarity_df, assay_df, on='Test_Code', how='left')
    df['search'] = (
        df['Instrument_ID'] + '-' + df['Specimen_ID'] + '*' + df['Test_Code']
        + '*SNV*.xlsx'
    )
    df['filename'] = df['search'].swifter.apply(find_file_name)
    df['path'] = df.apply(
        lambda x: create_path(x['filename'], x['Assay'], x['Run_name']),
        axis=1
    )

    # Create output files
    df.to_csv("all_data.csv")
    df.path.to_csv("paths.csv", header=False, index=False)

if __name__ == "__main__":
    main()
