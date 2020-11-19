#!/usr/bin/env python3
"""
Script to detect Tandem Repeats in a specified fasta file containing protein sequences using TRAL.
Output will be generated in a specified output directory (which will be made if it does not exist).
Output will consist of one .tsv file and one .pkl file for each protein in which a TR is detected. These files can be
merged into one file containing all Tandem Repeats by running the separate 'merge_tral_results.py' script.

NOTE: this script is heavily based on 'TR_in_multiple_Protein.py' by Matteo Delucchi (https://github.com/matteodelucchi/CRC_TRs)
"""

import sys
import os
import pickle
import argparse
import logging
import logging.config

from Bio import SeqIO

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence
from tral.repeat_list import repeat_list
from tral.hmm import hmm


def find_protein_repeats(sequences_file, result_dir):
    # TODO implement checking of potential results dir for proteins to skip
    logging.config.fileConfig(config_file("logging.ini"))
    log = logging.getLogger('root')

    CONFIG_GENERAL = configuration.Configuration.instance().config
    CONFIG = CONFIG_GENERAL["repeat_list"]
    score = CONFIG["model"]

    proteins = SeqIO.parse(sequences_file, "fasta")

    all_denovo_repeats = 0
    all_filtered_repeats = 0

    counter = 1
    for record in proteins:
        print("Started work on protein number {}".format(counter))
        counter += 1
        seq_name = record.id.split("|")[1]

        # name is protein identifier
        seq = sequence.Sequence(seq=str(record.seq), name=seq_name)

        denovo_list = seq.detect(denovo=True, repeat={"calc_pvalue": True})

        if not denovo_list and len(denovo_list) == 0:
            continue

        for TR in denovo_list.repeats:
            TR.model = None

        seq.set_repeatlist(denovo_list, "denovo_all")

        # add number of denovo found repeats
        all_denovo_repeats += len(seq.get_repeatlist("denovo_all").repeats)

        ##########################################################################
        # Filtering TRs
        denovo_list_filtered = filter_repeatlist(seq.get_repeatlist("denovo_all"))

        if not denovo_list_filtered or len(denovo_list_filtered.repeats) == 0:
            continue

        ##########################################################################
        # Clustering
        # De novo TRs are clustered for overlap (common ancestry). Only best =
        # lowest p-Value and lowest divergence are retained.
        denovo_list_filtered = denovo_list_filtered.filter(
            "none_overlapping", ["common_ancestry"], [("pvalue", "phylo_gap01"), ("divergence", "phylo_gap01")])
        seq.set_repeatlist(denovo_list_filtered, "denovo_filtered")

        ##########################################################################
        # De novo TRs are refined with HMMs
        final_list = refine_repeatlist(seq, "denovo_filtered")
        seq.set_repeatlist(final_list, "denovo_final")

        ##########################################################################
        # Save Tandem Repeats TODO: implement pickle binary dump
        # Create output directory if not already exists.
        try:
            if not os.path.isdir(result_dir):
                os.makedirs(result_dir)
        except:
            raise Exception(
                "Could not create path to result directory: {}".format(
                    os.path.dirname(result_dir)))

        # create filename
        output_tsv_file = os.path.join(result_dir, seq_name + ".tsv")

        # save TR-file in tsv format TODO: prevent writing if no TRs remain after filter and cluster
        write_file(seq.get_repeatlist("denovo_final"), output_tsv_file)

        all_filtered_repeats += len(seq.get_repeatlist("denovo_final"))
        print("\n***", seq_name, "***")
        print("denovo repeats:", len(denovo_list.repeats))
        print("repeats after filtering and clustering:",
              len(seq.get_repeatlist("denovo_final")))

        for i in range(len(seq.get_repeatlist("denovo_final"))):
            print(seq.get_repeatlist("denovo_final")[i])

    return print("\nThere where {} repeats found de novo.".format(all_denovo_repeats),
                 "After filtering and clustering there where only {} repeats left.\n".format(
                     all_filtered_repeats))


def filter_repeatlist(tr_list, pvalue_threshold=0.05, divergence_threshold=0.1, n_threshold=2.5):
    # filtering for pvalue
    tr_list_filtered = tr_list.filter(
        "pvalue",
        "phylo_gap01",
        pvalue_threshold)

    # filtering for divergence
    tr_list_filtered = tr_list_filtered.filter(
        "divergence",
        "phylo_gap01",
        divergence_threshold)

    # filtering for number of repeat units
    tr_list_filtered = tr_list_filtered.filter(
        "attribute",
        "n_effective",
        "min",
        n_threshold)

    return tr_list_filtered


def refine_repeatlist(seq, tag):
    refined_list = []
    for TR in seq.get_repeatlist(tag).repeats:
        use_refined = False
        denovo_hmm = hmm.HMM.create(input_format='repeat', repeat=TR)
        # Run HMM on sequence
        denovo_refined_list = seq.detect(lHMM=[denovo_hmm])
        if denovo_refined_list and denovo_refined_list.repeats:
            TR_refined = denovo_refined_list.repeats[0]
            TR_refined.TRD = TR.TRD
            TR_refined.model = "cpHMM"
            # Check whether new and old TR overlap. Check whether new TR is
            # significant. If not both, put unrefined TR into final.
            if repeat_list.two_repeats_overlap(
                    "shared_char",
                    TR,
                    TR_refined):
                tmp_repeatlist = repeat_list.RepeatList([TR_refined])
                tmp_repeatlist_filtered = filter_repeatlist(tmp_repeatlist)
                if tmp_repeatlist_filtered.repeats:
                    use_refined = True
        if use_refined:
            refined_list.append(TR_refined)
        else:
            refined_list.append(TR)
    return refined_list


def write_file(tr_list, destination):
    header = "\t".join(["begin",
                        "msa_original",
                        "l_effective",
                        "n_effective",
                        "repeat_region_length",
                        "divergence",
                        "pvalue"])
    with open(destination, "w") as f:
        f.write(header)
        for tr in tr_list:
            line = [
                str(i) for i in [
                    tr.begin,
                    ",".join(tr.msa),
                    tr.l_effective,
                    tr.n_effective,
                    tr.repeat_region_length,
                    tr.divergence("phylo_gap01"),
                    tr.pvalue("phylo_gap01")
                ]
            ]
            f.write("\n" + "\t".join(line))



def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta", "-f", type=str, required=True, help="Path to file containing protein sequence(s)"
    )
    parser.add_argument(
        "--outdir", "-o", type=str, required=True, help="Path to output directory"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parser()

    find_protein_repeats(args.fasta, args.outdir)
