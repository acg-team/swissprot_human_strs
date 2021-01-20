#!/usr/bin/env python
import argparse

from src.fasta_utils import ProteinRetriever


def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--genefile", "-g", required=True, type=str, help="Path to file containing one gene id per line."
    )
    parser.add_argument(
        "--outfile", "-o", required=True, type=str, help="File where protein fasta sequences will be written."
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parser()

    retriever = ProteinRetriever(args.genefile, args.outfile)
    retriever.retrieve_proteins()
