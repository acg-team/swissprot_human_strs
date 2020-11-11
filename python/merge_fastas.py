#!/usr/bin/env python
import argparse

from src.fasta_utils import FastaMerger


def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--file1", "-f1", type=str, required=True, help="First input fasta file"
    )
    parser.add_argument(
        "--file2", "-f2", type=str, required=True, help="Second input fasta file"
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Output destination for merged proteins to be written to"
    )
    parser.add_argument(
        "--nometa", "-n", action="store_true", help="Flag to prevent writing of metadata (default: False)"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parser()

    merger = FastaMerger(
        args.file1,
        args.file2,
        args.output
    )

    print()
    print("Merging the following files:\n{} ({} proteins) \n{} ({} proteins)\n".format(
        args.file1,
        len(merger.file1),
        args.file2,
        len(merger.file2))
    )
    merger.write_merged()
    print("Merged file written to:\n{} ({} proteins total)\n".format(args.output, len(merger.merged)))

    if not args.nometa:
        merger.write_meta_file()
        print("Meta information written to:\n{}".format(merger.out_meta))
