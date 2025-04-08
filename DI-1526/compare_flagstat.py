import argparse

import dxpy


def main(samples, path):
    for sample in samples:
        found_files = dxpy.find_data_objects(
            name=f"{sample}*", name_mode="glob", folder=path
        )

        files = [file for file in found_files]

        assert len(files) == 2, "Expected 2 flagstat files to be found"

        file1_data, file2_data = files

        file1_content = dxpy.open_dxfile(file1_data["id"]).read()
        file2_content = dxpy.open_dxfile(file2_data["id"]).read()

        if file1_content == file2_content:
            print(f"{sample} OK")
        else:
            print(f"{sample} not OK")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "samples", nargs="+", help="Sample ids to find flagstat outputs"
    )
    parser.add_argument(
        "-p",
        "--path",
        help="Path in which to find data. There should be 2 flagstat files in this folder or the subfolders",
    )

    args = parser.parse_args()
    main(args.samples, args.path)
