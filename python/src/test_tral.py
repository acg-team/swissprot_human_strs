#!/usr/bin/env python3
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
from tral.hmm import hmm


def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta", "-f", type=str, required=True, help="Path to file containing protein sequence(s)"
    )
    parser.add_argument(
        "--outdir", "-o", type=str, required=True, help="Path to output directory"
    )

    return parser.parse_args()


def find_protein_repeats(sequences_file, result_dir, pvalue_threshold=0.05, divergence_threshold=0.1, n_threshold=2.5, l_threshold=3):
    logging.config.fileConfig(config_file("logging.ini"))
    log = logging.getLogger('root')

    CONFIG_GENERAL = configuration.Configuration.instance().config
    CONFIG = CONFIG_GENERAL["repeat_list"]
    score = CONFIG["model"]

    proteins = SeqIO.parse(sequences_file, "fasta")

    all_denovo_repeats = 0
    all_filtered_repeats = 0

    for record in proteins:
        seq_name = record.id.split("|")[1]

        # name is protein identifier
        seq = sequence.Sequence(seq=str(record.seq), name=seq_name)

        denovo_list = seq.detect(denovo=True)

        if not denovo_list:
            continue

        for TR in denovo_list.repeats:
            TR.calculate_pvalues()

        ##########################################################################
        # Filtering TRs
        # add number of denovo found repeats
        all_denovo_repeats += len(denovo_list.repeats)

        # filtering for pvalue
        denovo_list = denovo_list.filter(
            "pvalue",
            score,
            pvalue_threshold)

        # filtering for divergence
        denovo_list = denovo_list.filter(
            "divergence",
            score,
            divergence_threshold)

        # filtering for number of repeat units
        denovo_list = denovo_list.filter(
            "attribute",
            "n_effective",
            "min",
            n_threshold)

        # filtering for length of repeat units
        # denovo_list = denovo_list_remastered.filter(
        #     "attribute",
        #     "l_effective",
        #     "max",
        #     l_threshold)

        ##########################################################################
        # Building HMM with hmmbuild
        # # De novo TRs are remastered with HMM
        denovo_hmm = [hmm.HMM.create(input_format='repeat', repeat=iTR)
                      for iTR in denovo_list.repeats]  # only possible with hmmbuild
        denovo_list_remastered = seq.detect(lHMM=denovo_hmm)

        if not denovo_list_remastered:
            continue

        # filtering for pvalue
        denovo_list_remastered = denovo_list_remastered.filter(
            "pvalue",
            score,
            pvalue_threshold)


        ##########################################################################
        # Clustering
        # De novo TRs were clustered for overlap (common ancestry). Only best =
        # lowest p-Value and lowest divergence were retained.
        denovo_list_remastered = denovo_list_remastered.filter(
            "none_overlapping", ["common_ancestry"], [("pvalue", score), ("divergence", score)])


        ##########################################################################
        # Save Tandem Repeats
        # Create output directory if not already exists.
        try:
            if not os.path.isdir(result_dir):
                os.makedirs(result_dir)
        except:
            raise Exception(
                "Could not create path to result directory: {}".format(
                    os.path.dirname(result_dir)))

        # create filename
        output_pickle_file = os.path.join(result_dir, seq_name + ".pkl")
        output_tsv_file = os.path.join(result_dir, seq_name + ".tsv")

        # save TR-file
        denovo_list_remastered.write(output_format="pickle", file=output_pickle_file)
        denovo_list_remastered.write(output_format="tsv", file=output_tsv_file)

        all_filtered_repeats += len(denovo_list_remastered.repeats)
        print("\n***", seq_name, "***")
        print("denovo repeats:", len(denovo_list.repeats))
        print("repeats after filtering and clustering:",
              len(denovo_list_remastered.repeats))

        for i in range(len(denovo_list_remastered.repeats)):
            print(denovo_list_remastered.repeats[i])

    return print("\nThere where {} repeats found de novo.".format(all_denovo_repeats),
                 "After filtering and clustering there where only {} repeats left.\n".format(
                     all_filtered_repeats))


if __name__ == "__main__":
    args = parser()

    find_protein_repeats(args.fasta, args.outdir)
