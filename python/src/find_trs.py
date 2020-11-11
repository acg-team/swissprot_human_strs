#!/usr/bin/env python
import os
import argparse

from Bio import SeqIO

from tral.sequence import sequence
from tral.hmm import hmm

__all__ = [
    "TRFinder",
]


class TRFinder(object):
    def __init__(self, fasta_file, output_dir, score="phylo_gap01"):
        if not fasta_file.endswith(".fasta"):
            raise ValueError("Input file needs to have '.fasta' extension.")
        self.fasta_handle = fasta_file
        self.sequences = SeqIO.parse(self.fasta_handle, "fasta")
        self.output_dir = self.generate_output_dir(output_dir)
        self.score = score

    def generate_output_dir(self, output_dir):
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        return output_dir

    def merge_repeat_files(self):
        file_list = [os.path.join(self.output_dir, f) for f in os.listdir(self.output_dir) if f.endswith(".tsv")]
        # get header
        with open(file_list[0], "r") as f:
            header_line = f.readline().strip()
            if not header_line.startswith("begin"):
                raise ValueError("File layout different than expected!")
            header_line = "ID\t{}\n".format(header_line)
        output_file = "{}/merged.tsv".format(self.output_dir)
        with open(output_file, "w") as o:
            o.write(header_line)
            for file_name in file_list:
                handle = file_name.split("/")[-1]
                prot_id = handle.split(".")[0]
                with open(file_name, "r") as f:
                    for line in f:
                        if not line.startswith("begin"):
                            o.write("{}\t{}".format(prot_id, line))

    def detect_in_sequence(self, record, remaster=True):
        seq_name = record.id.split("|")[1]
        # name is protein identifier
        seq = sequence.Sequence(seq=str(record.seq), name=seq_name)
        denovo_list = seq.detect(denovo=True)
        ##########################################################################
        # Building HMM with hmmbuild
        # De novo TRs are remastered with HMM
        if remaster:
            denovo_hmm = [hmm.HMM.create(input_format='repeat', repeat=iTR)
                          for iTR in denovo_list.repeats]  # only possible with hmmbuild
            denovo_list = seq.detect(lHMM=denovo_hmm)
        if not denovo_list:
            return
        for TR in denovo_list.repeats:
            TR.calculate_pvalues()
        return denovo_list

    def filter(self, repeat_list, criterion, threshold, **kwargs):
        try:
            if criterion in {"pvalue", "divergence"}:
                return repeat_list.filter(criterion, self.score, threshold)
            elif criterion in {"l_effective", "n_effective"}:
                return repeat_list.filter("attribute", criterion, kwargs["filter_type"], threshold)
        except KeyError:
            raise KeyError("Filter for repeat list called with wrong kwargs")
        raise ValueError("Filter for repeat list called without proper filtering criterion")

    def cluster(self, repeat_list):
        ##########################################################################
        # Clustering
        # De novo TRs are clustered for overlap (common ancestry). Only best =
        # lowest p-Value and lowest divergence are retained.
        repeat_list = repeat_list.filter(
            "none_overlapping", ["common_ancestry"], [("pvalue", self.score), ("divergence", self.score)])
        return repeat_list


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
    finder = TRFinder(fasta_file=args.fasta, output_dir=args.outdir)
    for record in finder.sequences:
        repeat_list = finder.detect_in_sequence(record)
        if not repeat_list:
            continue

        repeat_list = finder.filter(repeat_list=repeat_list, criterion="pvalue", threshold=0.05)
        repeat_list = finder.filter(repeat_list=repeat_list, criterion="divergence", threshold=0.1)
        repeat_list = finder.filter(repeat_list=repeat_list, criterion="n_effective", threshold=2.5, filter_type="min")
        repeat_list = finder.filter(repeat_list=repeat_list, criterion="l_effective", threshold=3, filter_type="max")
        if len(repeat_list.repeats) == 0:
            continue
        clustered_list = finder.cluster(repeat_list)

        # create filename
        seq_name = record.id.split("|")[1]
        output_pickle_file = os.path.join(finder.output_dir, seq_name + ".pkl")
        output_tsv_file = os.path.join(finder.output_dir, seq_name + ".tsv")

        # save TR-file
        clustered_list.write(output_format="pickle", file=output_pickle_file)
        clustered_list.write(output_format="tsv", file=output_tsv_file)

        print("Found {} TR(s) in protein {} (after filtering and clustering)".format(len(clustered_list.repeats), seq_name))
    finder.merge_repeat_files()
