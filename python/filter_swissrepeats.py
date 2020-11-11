#!/usr/bin/env python
import argparse

def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--swissrepeats", "-s", type=str, required=True, help="Data from swissrepeats study"
    )
    parser.add_argument(
        "--pathway", "-p", type=str, required=True, help="File containing list of pathway proteins."
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Output destination for filtered data to be written to"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parser()

    with open(args.pathway, "r") as f:
        proteins = {line.split("\t")[0] for line in f if not line.startswith("ID")}

    with open(args.output, "w") as o:
        with open(args.swissrepeats) as s:
            found = 0
            for line in s:
                if line.startswith("ID"):
                    o.write(line)
                    continue
                id = line.split(",")[0]
                if id in proteins:
                    o.write(line)
                    found += 1

    print("Number of proteins in pathway: {}".format(len(proteins)))
    print("Number of TRs found: {}".format(found))
