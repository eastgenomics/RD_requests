'''
This script strips the PID fields (C3 and D3) from any variant workbooks in
a given file path
'''
import argparse
import glob
from openpyxl import load_workbook


def parse_args():
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--path",
        type=str,
        help=("Path to workbooks"),
    )

    return parser.parse_args()


def strip_pid(file):
    """
    Strip PID from cells C3, and D3 (which are the first name and last name
    fields) for the workbook in the given filepath

    Args
    ------
    file (str): path to an xlsx file

    Returns
    -------
    None, changes input file
    """
    workbook = load_workbook(file)
    summary = workbook["summary"]
    summary['C3'] = None
    summary['D3'] = None
    workbook.save(file)


def main():
    args = parse_args()

    if not args.path.endswith('/'):
        args.path += '/'

    files = glob.glob(args.path + "*.xlsx")

    for file in files:
        strip_pid(file)


if __name__ == "__main__":
    main()
