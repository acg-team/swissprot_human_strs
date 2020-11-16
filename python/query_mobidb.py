#!/usr/bin/env python3
import argparse
import urllib.request, urllib.error
import json
from Bio import SeqIO


def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--file", "-f", type=str, required=True, help="Path to fasta file containing proteins"
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Path where output file will be created"
    )

    return parser.parse_args()


def query_mobidb(sp_id, include_curated):
    # Check https://dev.mobidb.org/help/apidoc for information on MobiDB API
    filter = "&projection=prediction-disorder-mobidb_lite"
    if include_curated:
        filter += ",curated-disorder-priority"
    base_url = "https://mobidb.bio.unipd.it/api/download?acc={}{}&format=json"
    prot_url = base_url.format(sp_id, filter)
    try:
        with urllib.request.urlopen(prot_url) as response:
            disorder = response.read().decode("utf-8")
            disorder_dict = json.loads(disorder)
    except urllib.error.URLError:
        raise urllib.error.URLError("SSL: CERTIFICATE_VERIFY_FAILED")

    print("full-disorder-priority" in disorder_dict.keys())
    return disorder_dict


def disorder_from_full(dis_dict):
    region_lists = dis_dict["mobidb_consensus"]["disorder"]["full"][0]["regions"]
    disorder_regions = [(i[0], i[1]) for i in region_lists if i[2].lower() == "d"]
    return disorder_regions


def disorder_from_predictions(dis_dict):
    all_consensus = dis_dict["mobidb_consensus"]["disorder"]["predictors"]
    for consensus in all_consensus:
        if not consensus["method"] == "mobidb_lite":
            continue
        disorder_regions = [(i[0], i[1]) for i in consensus["regions"]]
        return disorder_regions


def main():
    # Uncomment for debugging
    # records = SeqIO.parse("../data/proteins/merged_wnt.fasta", "fasta")
    #
    # # Check https://mobidb.bio.unipd.it/about/mobidb for information on MobiDB API
    # for record in records:
    #     seq_id = record.id.split("|")[1]
    #     try:
    #         disorder_dict = query_mobidb(seq_id, include_curated=False)
    #     except urllib.error.URLError:
    #         print("{} was not found in MobiDB".format(seq_id))
    #         continue
    #
    #     try:
    #         disorder_regions = [(i[0], i[1]) for i in disorder_dict["prediction-disorder-mobidb_lite"]["regions"]]
    #     except KeyError:
    #         print("No disorder predicted for {}".format(seq_id))
    #         continue
    #     print(seq_id, disorder_regions)
    # exit()

    args = parser()
    records = SeqIO.parse(args.file, "fasta")

    # Check https://mobidb.bio.unipd.it/about/mobidb for information on MobiDB API
    with open(args.output, "w") as o:
        header = "ID\tbegin\tend"
        o.write(header)
        for record in records:
            seq_id = record.id.split("|")[1]
            print("Started work on sequence '{}'...".format(seq_id))
            try:
                disorder_dict = query_mobidb(seq_id, include_curated=False)
            except urllib.error.URLError:
                print("'Not found in MobiDB".format(seq_id))
                continue

            try:
                disorder_regions = [(i[0], i[1]) for i in disorder_dict["prediction-disorder-mobidb_lite"]["regions"]]
            except KeyError:
                print("No disorder predicted".format(seq_id))
                continue
            print("Predicted disordered regions: {}".format(disorder_regions))
            for coords in disorder_regions:
                o.write("\n{}\t{}\t{}".format(seq_id, coords[0], coords[1]))


if __name__ == "__main__":
    main()
