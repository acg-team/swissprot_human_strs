#!/usr/bin/env python3
"""
Script to retrieve information on signal peptides form SwissProt for all UniProt accession in a file (newline separated)
More specifically, it will retrieve the first and last amino acid position of a signal peptide if a protein contains
one, and write it to a file

Author: Max Verbiest
Contact: max.verbiest@zhaw.ch
"""

import requests
import xmltodict
import argparse


def get_signal_locations(accessions, destination):
    """Retrieve signal peptide begin and end sites from UniProt for all accessions and write them to a file.
    Output file will have format ID\tbegin\tend\n

    Paramters
    accessions (list, set, tuple):
                Iterable containing accessions to be checked for signal peptides
    destination (str):
                Output file where information on signal peptides will be written
    """
    
    with open(destination, "w") as o:
        counter = 1
        for accession in accessions:
            signal_peptide = None
            query = "https://www.uniprot.org/uniprot/{}.xml".format(accession)

            response_text = requests.get(query).content
            xml_dict = xmltodict.parse(response_text)
            try:
                features = xml_dict["uniprot"]["entry"]["feature"]
            except KeyError:
                print("No features were available for accession '{}'".format(accession))
                continue

            for i in features:
                try:
                    if i["@type"] == "signal peptide":
                        begin = i["location"]["begin"]["@position"]
                        end = i["location"]["end"]["@position"]
                        try:
                            if i["@description"] == "Not cleaved":
                                signal_peptide = (begin, end, "retained")
                        except KeyError:
                            signal_peptide = (begin, end, "cleaved")
                except:
                    # print("Problem parsing features for accession '{}'".format(accession))
                    continue
            if signal_peptide:
                # print(accession, signal_peptide)
                o.write("{}\t{}\t{}\t{}\n".format(accession, signal_peptide[0], signal_peptide[1], signal_peptide[2]))
            # Added counter to keep track of how far the script has come, makes it easier to restart after potential
            # crashes, which does happen sometimes
            if counter % 20 == 0:
                print("finished query for accession: {} (number {})".format(accession, counter))
            counter += 1


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--protein_ids", "-p", required=True, type=str, help="Path to file containing one UniProt accession per line."
    )
    parser.add_argument(
        "--outfile", "-o", required=True, type=str, help="File where signal peptide annotation will be written."
    )
    return parser.parse_args()


def main():
    debug = False
    if debug:
        print("RUNNING IN DEBUG MODE")
        accession_list = ["O95445", "P36894"]  # [SP retained, SP cleaved]
        for accession in accession_list:
            query = "https://www.uniprot.org/uniprot/{}.xml".format(accession)

            response_text = requests.get(query).content
            xml_dict = xmltodict.parse(response_text)
            try:
                features = xml_dict["uniprot"]["entry"]["feature"]
            except KeyError:
                print("No features were available for accession '{}'".format(accession))
                continue
            print("")
        exit()

    args = parser()
    with open(args.protein_ids, "r") as f:
        accession_list = [line.strip().split("\t")[0] for line in f]

    get_signal_locations(accession_list, args.outfile)


if __name__ == "__main__":
    main()
