#!/usr/bin/env python

import argparse
import os


def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--dir", "-d", type=str, required=True, help="Directory with tral output files."
    )

    return parser.parse_args()


def concatenate_files(file_list, output_file):
    # get header
    with open(file_list[0], "r") as f:
        header_line = f.readline().strip()
        if not header_line.startswith("begin"):
            raise ValueError("File layout different than expected!")
        header_line = "ID\t{}\n".format(header_line)

    with open(output_file, "w") as o:
        o.write(header_line)
        for file_name in file_list:
            handle = file_name.split("/")[-1]
            prot_id = handle.split(".")[0]
            with open(file_name, "r") as f:
                for line in f:
                    if not line.startswith("begin"):
                        o.write("{}\t{}".format(prot_id, line))


if __name__ == "__main__":
    args = parser()
    filenames = [os.path.join(args.dir, file) for file in os.listdir(args.dir) if file.endswith(".tsv")]
    output_file = "{}/merged.tsv".format(args.dir)
    concatenate_files(filenames, output_file)
