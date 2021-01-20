#!/usr/bin/env python
import os
import urllib.request
from Bio import SeqIO

__all__ = [
    "ProteinRetriever",
    "FastaMerger",
]

class ProteinRetriever(object):
    def __init__(self, file_path, output_file):
        self.file_path = file_path
        self.output_file = self.check_output_dir(output_file)
        self.gene_ids = self.parse_gene_ids()

    def parse_gene_ids(self):
        try:
            with open(self.file_path, "r") as gene_file:
                gene_ids = {line.strip() for line in gene_file}
        except FileNotFoundError:
            exit("ERROR: file '{}' does not exist".format(self.file_path))
        return gene_ids

    def check_output_dir(self, output_file):
        if output_file and not "/" in output_file:
            # file is specified in cwd
            return output_file
        # check if output directory exists
        output_dir = output_file.split("/")[:-1]
        output_dir = "/".join(output_dir)
        if not os.path.isdir(output_dir):
            exit("Specified output directory for file {} does not exist".format(output_file))
        return output_file

    def proteins_from_gene(self, gene_id):
        url = "https://www.uniprot.org/uniprot/?query=gene:{}&format=fasta&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes".format(gene_id)
        try:
            with urllib.request.urlopen(url) as response:
                    proteins = response.read().decode('utf-8')
                    return proteins
        except:
            # print("Couldn't retrieve proteins of gene: {} \n Connection issue or no proteins in Swissprot for this gene.".format(gene_id))
            pass

    def retrieve_proteins(self):
        print("Started retrieving proteins for file {}".format(self.file_path))
        with open(self.output_file, "w") as out_file:
            # [out_file.write(self.proteins_from_gene(gene_id)) for gene_id in self.gene_ids]
            for gene_id in self.gene_ids:
                proteins = self.proteins_from_gene(gene_id)
                if not proteins:
                    print("WARNING: COULD NOT RETRIEVE ANY PROTEINS FOR GENE '{}'".format(gene_id))
                    continue
                out_file.write(proteins)
        print("Proteins retrieved and written to {}".format(self.output_file))

    def print_gene_ids(self):
        for gene in self.gene_ids:
            print(gene)

    def __str__(self):
        return "ProteinRetriever for file: {}".format(self.file_path)


class FastaMerger(object):
    def __init__(self, file1, file2, out_file):
        self.out_file = out_file
        self.out_meta = self.make_meta_file_name()
        self.file1 = SeqIO.to_dict(SeqIO.parse(file1, "fasta"))
        self.file2 = SeqIO.to_dict(SeqIO.parse(file2, "fasta"))
        self.merged = self.merge_fastas()

    def merge_fastas(self):
        return {**self.file1, **self.file2}

    def write_merged(self):
        records_list = [record for record in self.merged.values()]
        SeqIO.write(records_list, self.out_file, "fasta")

    def make_meta_file_name(self):
        if not self.out_file.endswith(".fasta"):
            raise ValueError("Please specify an output .fasta file for proteins to be written to.")
        return self.out_file.replace(".fasta", "_meta.tsv")

    def write_meta_file(self):
        header = "ID\tFile1\tFile2"
        base_string = "\n{}\t{}\t{}"
        all_ids = sorted([prot_id for prot_id in self.merged.keys()])

        with open(self.out_meta, "w") as f:
            f.write(header)
            for prot_id in all_ids:
                try:
                    self.file1[prot_id]
                    try:
                        self.file2[prot_id]
                        filled_string = base_string.format(prot_id, "yes", "yes")
                    except KeyError:
                        filled_string = base_string.format(prot_id, "yes", "no")
                except KeyError:
                    filled_string = base_string.format(prot_id, "no", "yes")
                f.write(filled_string)

    def __str__(self):
        return "\
        File1 proteins before merge: {}\n\
        File2 proteins before merge: {}\n\
        Number proteins after merge: {}".format(len(self.file1), len(self.file2), len(self.merged))
