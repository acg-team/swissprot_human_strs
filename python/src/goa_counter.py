#!/usr/bin/env python3
import os


class GoaCounter(object):
    """
    Class designed to parse a Gene Ontology Annotation (GOA) file and then determine how many proteins
    from the Protein Atlas Pathology Atlas map to the different GO terms (split per anatomical site and
    correlation type)
    """

    def __init__(self, goa_file_handle, pa_file_handle):
        self.goa_file_handle = self.check_file(goa_file_handle)
        self.pa_file_handle = self.check_file(pa_file_handle)

    def check_file(self, file_handle):
        if not isinstance(file_handle, str):
            raise ValueError("Specify file path as string, not '{}'".format(type(file_handle)))
        if not os.path.exists(file_handle):
            raise ValueError("File '{}' does not exist".format(file_handle))
        return file_handle

    def parse_goa_file(self, skip_electronic=True, domains=None):
        if not domains:
            domains = {"C", "F", "P"}
        if not isinstance(domains, (list, set, tuple)) or not (all(i in {"C", "F", "P"} for i in domains)):
            raise ValueError("Specify desired GO domains in a list, set, or tuple with element(s) of {C, F, P}")

        goa_dict = dict()
        with open(self.goa_file_handle, "r") as f:
            for line in f:
                if not line.startswith("!"):
                    line_split = line.strip().split("\t")

                    # Check if filter criteria are met
                    evidence, domain = line_split[6], line_split[8]
                    if evidence == "IEA" and skip_electronic:
                        continue
                    if domain not in domains:
                        continue

                    # Format GO term as '{domain}_{GO_term}'
                    go_term = "{}_{}".format(domain, line_split[4])
                    uniprot_accession = line_split[1]

                    # Add UniProt accession to GO term in dict. If GO term not in dict: make new entry
                    try:
                        goa_dict[go_term].add(uniprot_accession)
                    except KeyError:
                        goa_dict[go_term] = {uniprot_accession}
        return goa_dict

    def parse_pa_file(self):
        f = open(self.pa_file_handle, "r")
        # Get header, extract list of PA groups and initialize pa_dict
        header = next(f)
        pa_groups = header.strip().split("\t")[1:]
        pa_dict = {group: set() for group in pa_groups}

        # iterate over rest of lines in file, add prot ids to pa_dict
        for line in f:
            line_split = line.strip().split("\t")
            current_id = line_split.pop(0)
            for i in range(0, len(line_split)):
                if int(line_split[i]) == 1:
                    pa_dict[pa_groups[i]].add(current_id)

        f.close()
        return pa_dict

    def write_goa_counts_per_group(self, output_file, goa_dict=None, pa_dict=None, scale=True):
        if not goa_dict:
            goa_dict = self.parse_goa_file()
        if not pa_dict:
            pa_dict = self.parse_pa_file()

        go_terms = sorted(goa_dict.keys())
        pa_groups = sorted(pa_dict.keys())
        with open(output_file, "w") as o:
            # write header
            o.write("pa_group\t{}\n".format("\t".join(go_terms)))
            for group in pa_groups:
                # counts_list = [str(len(pa_dict[group].intersection(goa_dict[term]))) for term in go_terms]
                counts_list = []
                for term in go_terms:
                    count = len(pa_dict[group].intersection(goa_dict[term]))
                    if scale:
                        # divide count by total number of proteins found for group, try to offset effect of group size
                        count = str(count / len(pa_dict[group]))
                    counts_list.append(count)
                o.write("{}\t{}\n".format(group, "\t".join(counts_list)))

    def __str__(self):
        return "GoaCounter using GOA file at:\n'{}'\nand PA file at:\n'{}'".format(self.goa_file_handle,
                                                                                   self.pa_file_handle)

    def __repr__(self):
        return "GoaCounter(goa_file_handle={}, pa_file_handle={})".format(self.goa_file_handle, self.pa_file_handle)


def main():
    goa_handle = "/Users/maxverbiest/PhD/projects/SP_CRC_Pathway_TRs/data/GO_PA_counts/goa_human.gaf"
    pa_handle = "/Users/maxverbiest/PhD/projects/SP_CRC_Pathway_TRs/data/pathology_atlas/pa_str_proteins_wide.tsv"
    counter = GoaCounter(goa_file_handle=goa_handle, pa_file_handle=pa_handle)
    # goa_dict = counter.parse_goa_file()
    # pa_dict = counter.parse_pa_file()

    goa_dict = counter.parse_goa_file(domains={"P"})
    counter.write_goa_counts_per_group(
        output_file="/Users/maxverbiest/PhD/projects/SP_CRC_Pathway_TRs/data/GO_PA_counts/goa_counts_norm_P.tsv",
        goa_dict=goa_dict
    )


if __name__ == "__main__":
    main()
