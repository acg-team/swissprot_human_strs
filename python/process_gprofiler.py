"""
Process g:Profiler output to include log2(fold change) and number/fraction of IDR containing proteins per enriched term.
IMPORTANT:
Table downloaded from g:Profiler website comes with quotation marks around each element and is comma separated.
Handling this would be very tedious to implement in this script, so I suggest preprocessing the table using the
following R commands:
    tab <- read.table("gProfiler_file.csv", sep=",", header=TRUE)
    write.table(tab, "reformat.tsv", sep="\t", row.names=FALSE, quote=FALSE)

Author: Max Verbiest
Contact: max.verbiest@zhaw.ch
"""
import argparse
import numpy as np


def process_gprofiler(gprofiler, outfile, idrs=False, signals=False, idr_list=None, signal_file=None):
    """Take a prepocessed, tsv gProfiler output file and add columns for specified measures:
        default: log2(fold change) optional: idr fraction, signal peptide fraction
    Measures will be calculated, added and the processed gProfiler file is written to outfile

    Parameters
    gprofiler (str):    Path to preprocessed gProfiler output (check top of file)
    outfile (str):      Path specifying where processed file will be deposited
    idrs (Bool):        Whether fraction of IDR containing proteins in term should be calculated (requires idr_list)
    signals (str):      Whether fraction of signal peptide containing proteins in term should be calculated (requires signal_file)
    idr_list (str):     Path to file of IDR containing accessions (one per line)
    signal_file (str):  Path to file of signal peptide containing accessions (as produced by signal_peptide_annotation.py)
    """

    if idrs:
        if not idr_list:
            raise (Exception("Need a list of IDR proteins to calculate fraction IDR containing per term"))
        with open(idr_list, "r") as f:
            idr_set = {line.strip() for line in f}
    if signals:
        if not signal_file:
            raise (Exception("Need a file of signal peptide proteins to calculate fraction signal peptide containing per term"))
        with open(signal_file, "r") as f:
            signal_set = {line.split("\t")[0] for line in f}

    f = open(gprofiler, "r")
    header = next(f)
    header_list = ["log2fc"]
    if idrs:
        header_list.append("idr_fraction")
    if signals:
        header_list.append("signal_fraction")
    new_header = extend_line(header, *header_list)
    with open(outfile, "w") as o:
        o.write(new_header)
        for line in f:
            item_list = [get_log2fc(line)]
            if idrs:
                item_list.append(get_feature_fraction(line, idr_set))
            if signals:
                item_list.append(get_feature_fraction(line, signal_set))
            new_line = extend_line(line, *item_list)
            o.write(new_line)
    f.close()


def get_log2fc(line):
    """Calculate the log2(fold change) as a measure of how overrepresented the term is in the gene list

    Parameters
    line (str):     Line from a gProfiler output file (tab separated and without quotation marks, see top of file)

    Returns
    log2fc (float):
                    Log2(fold change) calculated for term
    """

    line_split = line.split("\t")

    observed = int(line_split[7])
    expect = int(line_split[5]) / int(line_split[8]) * int(line_split[6])

    fold_change = observed/expect
    return np.log2(fold_change)


def get_feature_fraction(line, feature_accessions):
    """Calculate the fraction of terms from the intersection of over-represented term (10th column in gProfiler output)
    are also in the set of accessions that represent proteins containing a certain feature (e.g. IDR, signal peptide)

    Parameters
    line (str):     Line from a gProfiler output file (tab separated and without quotation marks, see top of file)
    feature_accessions (set):
                    Set containing accessions that possess the feature of interest

    Returns
    fraction (int):
                    Fraction of accessions in the overrepresented term that contain the feature
    """

    term_accessions = line.split("\t")[9]
    term_accessions = {accession for accession in term_accessions.split(",")}

    with_feature = len(term_accessions.intersection(feature_accessions))
    return with_feature / len(term_accessions)


def extend_line(line, *args):
    """Add elements to a line, separated by tabs

    Parameters:
    line (str):     Line that *args will be added to for tsv format
    args (str):     Elements to add to line

    Returns:
    new_line (str):
                    Input line extended with *args
    """

    new_line = line.strip()
    for i in args:
        new_line = new_line + "\t" + str(i)
    return new_line + "\n"


def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--file", "-f", required=True, type=str, help="Output of g:Profiler web tool to be processed"
    )
    parser.add_argument(
        "--outfile", "-o", required=True, type=str, help="File where output will be deposited"
    )
    parser.add_argument(
        "--with_idr", "-w", action="store_true", help="Determine fraction of IDR containing proteins per term"
    )
    parser.add_argument(
        "--with_signals", "-s", action="store_true", help="Determine number of signal peptide containing proteins per term"
    )
    parser.add_argument(
        "--idr_list", "-i", type=str,
        default="/Users/maxverbiest/PhD/projects/SP_CRC_Pathway_TRs/results/disorder/local_annotations/unique_idr_containing.tsv",
        help="File containing UniProt accessions of IDR containing proteins (one per line) -- has default"
    )
    parser.add_argument(
        "--signal_peptides", "-p", type=str,
        default=None,
        help="File containing locations of signal peptides in UniProt accessions (as created by signal_peptide_annotation.py)"
    )

    return parser.parse_args()


def main():
    args = cla_parser()
    process_gprofiler(gprofiler=args.file,
                      outfile=args.outfile,
                      idrs=args.with_idr,
                      signals=args.with_signals,
                      idr_list=args.idr_list,
                      signal_file=args.signal_peptides)


if __name__ == "__main__":
    main()
