#!/usr/bin/env python3
"""
Sometimes, TRAL returns duplicate TRs (not sure why).
This script will take the output from run_tral.py that has been merged using merge_tral_results.py and remove the
duplicates. This is done by checked if (ID, begin) is unique (pretty naive method but should work)

Author: Max Verbiest
Contact: max.verbiest@zhaw.ch
"""

import argparse


def tral_file_filter(file):
    """In file, check whether there are TRs that belong to the same protein and have the same starting amino acid.
    This should never occur and these instances are therefore considered duplicates.

    Paramters
    file (str): Merged TRAL output file that will be checked for duplicates

    Returns
    line_list (list[str]):
                A list containing lines describing tandem repeats, with the duplicate lines removed
    """
    
    filtered_dict = {}
    line_list = []
    dupe_count = 0
    with open(file, "r") as f:
        for line in f:
            if line.startswith("ID"):
                line_list.append(line)
                continue
            prot_id, begin = line.split("\t")[0], line.split("\t")[1]
            try:
                if begin in filtered_dict[prot_id]:
                    dupe_count += 1
                    continue
                filtered_dict[prot_id].append(begin)
                line_list.append(line)
            except KeyError:
                filtered_dict[prot_id] = [begin]
                line_list.append(line)
    print("Filter found {} cases where TRs in the same protein had the same start position.".format(dupe_count))
    return line_list


def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--file", "-f", type=str, required=True, help="File with merged TRAL results to be filtered."
    )
    parser.add_argument(
        "--overwrite", "-o", type=bool, default=False, help="Overwrite existing file with filtered file?"
    )

    return parser.parse_args()


def main():
    args = parser()
    filtered_lines = tral_file_filter(args.file)
    # overwrite input file?
    if args.overwrite:
        print("Overwriting existing file with filtered file")
        with open(args.file, "w") as o:
            for line in filtered_lines:
                o.write(line)


if __name__ == "__main__":
    main()
