#!/usr/bin/env python3

import argparse
import itertools

from Bio import SeqIO
import pandas as pd

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--proteins", "-p", type=str, required=True, help="Path to where SwissProt proteins are stored"
    )
    parser.add_argument(
        "--str_df", "-s", type=str, required=True, help="Path to where dataframe with STR information is stored"
    )

    return parser.parse_args()


def main():
    proteins = "/Users/maxverbiest/PhD/projects/SP_CRC_Pathway_TRs/data/uniprot_human_reviewed_21122020.fasta"
    strs = "/Users/maxverbiest/PhD/projects/SP_CRC_Pathway_TRs/data_for_sub/str_sp_final.tsv"
    duos = {
        ("G", "S"),
        ("A", "P"),
        ("G", "R"),
        ("R", "S"),
        ("G", "P"),
        ("E", "D"),
        ("G", "A"),
        ("P", "S"),
        ("E", "R"),
        ("K", "E")
    }
    count_dict = dict()
    for duo in duos:
        count_dict[f"{duo[0]}{duo[1]}"] = 0
        count_dict[f"{duo[1]}{duo[0]}"] = 0

    str_df = pd.read_csv(strs, sep="\t")[["ID", "begin", "end"]]
    prot_dict = dict()
    for i, row in str_df.iterrows():
        try:
            prot_dict[row["ID"]].append((row["begin"], row["end"]))
        except KeyError:
            prot_dict[row["ID"]] = [(row["begin"], row["end"])]

    prot_parser = SeqIO.parse(proteins, "fasta")
    for record in prot_parser:
        uniprot_id = record.id.split("|")[1]
        for i in range(0, len(record.seq) + 1):
            try:
                dipeptide = record.seq[i:i + 2]
                skip = False
                # try to access count dictionary to see if dipeptide is among those to count
                count_dict[dipeptide]
                try:
                    str_coords = prot_dict[uniprot_id]
                    # protein contains STR, check if current dipeptide overlaps STR
                    for coords in str_coords:
                        # is first dipeptide position in STR?
                        if coords[0] <= i+1 <= coords[1]:
                            skip = True
                            break
                        # is second dipeptide position in STR?
                        elif coords[0] <= i+2 <= coords[1]:
                            skip = True
                            break
                    # dipeptide not in any STR, increment dipeptide count
                    if not skip:
                        count_dict[dipeptide] += 1
                except KeyError:
                    # protein does not have STR(s), increment dipeptide count
                    count_dict[dipeptide] += 1
            except KeyError:
                # dipeptide is not among those to count, continue
                continue
    print(count_dict)

if __name__ == "__main__":
    main()
