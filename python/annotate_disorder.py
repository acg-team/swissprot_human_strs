#!/usr/bin/env/python
import argparse

from src.disorder_annotation import APIDisorderAnnotator
from src.disorder_annotation import LocalDisorderAnnotator


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fasta", type=str, required=True, help="Fasta file with proteins to be annotated"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Path to output file where annotations will be written"
    )
    parser.add_argument(
        "-m", "--mode", type=str, default="api", choices=["local", "api"],
        help="Choose whether to run MobiDB-Lite locally (local) or retrieve annotations from MobiDB web API (api)"
    )
    parser.add_argument(
        "-t", "--threads", type=int, choices=range(0, 8),
        help="(only relevant in local mode) How many threads should MobiDB-Lite use? (between 1(default) and 7)"
    )
    return parser.parse_args()


def main():
    args = parser()

    if args.mode == "local":
        annotator = LocalDisorderAnnotator(args.fasta, args.output)
        annotator.get_disorder_annotations(threads=args.threads)
    elif args.mode == "api":
        annotator = APIDisorderAnnotator(args.fasta, args.output)
        annotator.get_disorder_annotations()


if __name__ == "__main__":
    main()
