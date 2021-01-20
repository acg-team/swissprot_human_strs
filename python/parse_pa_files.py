"""
Parsing Pathology Atlas favourable & unfavourable files, make combined file containing UniProt accession and
whether it is favourable/ unfavourable (csv format).
"""

import argparse


def parse_pa_file(pa_file, correlation):
    if not correlation in {"pa_fav", "pa_unfav"}:
        raise ValueError("ERROR: correlation type must be element of {pa_fav, pa_unfav}")

    accession_dict = {}
    f = open(pa_file, "r")
    # file has header, skip first line
    next(f)
    for line in f:
        uniprot = line.split("\t")[4]
        if not uniprot:
            continue
        # process string of accession(s), add to dict
        uniprot = uniprot.replace(" ", "")
        if "\"" in uniprot:
            uniprot = uniprot.replace("\"", "")
        uniprot_split = uniprot.split(",")
        for accession in uniprot_split:
            accession_dict[accession] = correlation
    f.close()
    return accession_dict


def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--unfav", "-u", type=str, required=True, help="Pathology atlas file with unfavourable correlating genes"
    )
    parser.add_argument(
        "--fav", "-f", type=str, required=True, help="Pathology atlas file with favourable correlating genes"
    )
    parser.add_argument(
        "--site", "-s", type=str, required=True, help="Primary site of cancer"
    )
    parser.add_argument(
        "--outdir", "-d", type=str, required=False, default="./", help="Directory to deposit output file"
    )

    return parser.parse_args()


def main():
    args = cla_parser()

    fav_uniprot = parse_pa_file(args.fav, "pa_fav")
    unfav_uniprot = parse_pa_file(args.unfav, "pa_unfav")

    if not args.outdir.endswith("/"):
        args.outdir = args.outdir + "/"
    output_file = "{}/{}_uniprot_only.csv".format(args.outdir, args.site)
    with open(output_file, "w") as o:
        o.write("ID,pa_{}".format(args.site))
        # concatenate dicts, write to file (assume no accessions are in both dictionaries)
        for key, value in {**unfav_uniprot, **fav_uniprot}.items():
            o.write("\n{},{}".format(key, value))


if __name__ == "__main__":
    main()
